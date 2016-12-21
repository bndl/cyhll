# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from math import ceil, log

import array
from cpython cimport array

from .consts import thresholds, estimates, bias

from .murmer3 cimport hash_x64_128
cdef inline long hash64(value):
    return hash_x64_128(value)[0]


cdef unsigned long MAXINT = -(2u ** (64u - 1u))


cpdef double get_alpha(int p) except -1:
    if not (4 <= p <= 16):
        raise ValueError("p=%d should be in range [4 : 16]" % p)
    if p == 4:
        return 0.673
    if p == 5:
        return 0.697
    if p == 6:
        return 0.709
    return 0.7213 / (1.0 + 1.079 / (1 << p))


cpdef inline long get_rho(long w, long max_width) except -1:
    cdef long rho = max_width - w.bit_length() + 1
    if rho <= 0:
        raise ValueError('w overflow')
    return rho


cdef estimate_bias(E, p):
    nearest_neighbors = get_nearest_neighbors(E, estimates[p - 4])
    bias_vector = bias[p - 4]
    return sum(bias_vector[i] for i in nearest_neighbors) / len(nearest_neighbors)


cdef get_nearest_neighbors(E, estimate_vector):
    distance_map = sorted(
        ((E - val) ** 2, idx)
        for idx, val in enumerate(estimate_vector)
    )
    return [idx for dist, idx in distance_map[:6]]


cdef inline _combine(HyperLogLog a, HyperLogLog b):
    cdef int i
    for i in range(a.m):
        a.Mv[i] = max(a.Mv[i], b.Mv[i])


cdef class HyperLogLog:
    cdef public float error_rate
    cdef public double alpha
    cdef public int p
    cdef public int m
    cdef public array.array M
    cdef long* Mv
    cdef public long mask
    cdef public long max_width

    def __init__(self, error_rate, M=None):
        assert 0 < error_rate < 1
        self.error_rate = error_rate
        self.p = ceil(log((1.04 / error_rate) ** 2, 2))
        self.alpha = get_alpha(self.p)
        self.m = 1 << self.p
        self.M = M or array.array('l', [0] * self.m)
        self.Mv = self.M.data.as_longs
        self.mask = self.m - 1
        self.max_width = 64 - self.p


    def add(self, value):
        cdef unsigned long x = hash64(value) + MAXINT
        cdef long j = x & self.mask
        cdef long w = x >> self.p
        self.Mv[j] = max(self.Mv[j], get_rho(w, self.max_width))
        return self


    def add_all(self, values):
        cdef unsigned long x
        cdef long j
        cdef long w

        cdef int mask = self.mask
        cdef int max_width = self.max_width

        for value in values:
            x = hash64(value) + MAXINT
            j = x & mask
            w = x >> self.p
            self.Mv[j] = max(self.Mv[j], get_rho(w, max_width))

        return self


    def merge(self, *others):
        for other in others:
            if self.m != other.m:
                raise ValueError('Counters precisions should be equal')

        for other in others:
            _combine(self, other)

        return self


    def _Ep(self):
        cdef Py_ssize_t i
        E = self.alpha * self.m ** 2 / sum((2.0 ** -self.Mv[i]) for i in range(self.m))
        return (E - estimate_bias(E, self.p)) if E <= 5 * self.m else E


    def card(self):
        V = self.M.count(0)
        if V > 0:
            H = self.m * log(self.m / V)
            return H if H <= thresholds[self.p - 4] else self._Ep()
        else:
            return self._Ep()
        
        
    def __richcmp__(self, o, op):
        assert isinstance(o, HyperLogLog)
        if op == 0:
            return self.card() < o.card()
        elif op == 1:
            return self.card() <= o.card()
        elif op == 2:
            return self.M == o.M
        elif op == 3:
            return self.M != o.M
        elif op == 4:
            return self.card() > o.card()
        elif op == 5:
            return self.card() >= o.card()
        else:
            assert False


    def __reduce__(self):
        return HyperLogLog, (self.error_rate, self.M)


    def __repr__(self):
        return '<HyperLogLog alpha=%s, p=%s, m=%s>' % (self.alpha, self.p, self.m)
