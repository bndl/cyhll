'''
This module implements the MurmurHash3 algorithm in Cython.

The original source file states:

    MurmurHash3 was written by Austin Appleby, and is placed in the public
    domain. The author hereby disclaims copyright to this source code.
'''


from libc.math cimport ceil
from libc.stdlib cimport malloc, calloc, free

from cython.view cimport memoryview
from builtins import memoryview as py_memoryview

from cpython cimport array, bytearray
import array

from cpython cimport array, bytes, buffer, PyBytes_AsString, PyBytes_AsStringAndSize, \
    PyObject_GetBuffer, PyBytes_FromFormat, PyInt_AsLong, PyInt_GetMax, PyFloat_AsDouble

cdef extern from *:
    char* PyUnicode_AsUTF8AndSize(object o, Py_ssize_t* size)


MIN_LONG = -2 ** 63
MAX_LONG = 2 ** 63 - 1


cdef inline ulong rotl64 (ulong x, char r):
    return (x << r) | (x >> (64 - r))


cdef inline ulong fmix(ulong k):
    k ^= k >> 33
    k *= <ulong> 0xff51afd7ed558ccd
    k ^= k >> 33
    k *= <ulong> 0xc4ceb9fe1a85ec53
    k ^= k >> 33
    return k


cdef inline ulong get_block(char* data, Py_ssize_t i):
    return <ulong> data[i]


cpdef tuple hash_x64_128(obj):
    '''
    Takes a 128 bit murmer3 hash of obj.

    Args:
        obj : bytes, bytearray, memoryview, str, array.array, int, float, tuple or list
            The object to take the hash of. If obj is a tuple or a list, it's elements are hashed
            recursively, XORing the intermediate hashes.

    Returns:
        A 2-tuple of the two 64 bit parts of the 128 bit hash.
    '''
    cdef ulong* target = <ulong*> malloc(2 * sizeof(ulong))
    if not target:
        raise MemoryError()

    try:
        _hash_x64_128(obj, target)
        h1 = <long> target[0]
        h2 = <long> target[1]
        return h1, h2
    finally:
        free(target)


cdef void _hash_x64_128(obj, ulong* target):
    cdef char* data
    cdef bint do_free = False
    cdef int seed
    cdef Py_ssize_t size
    cdef Py_buffer buf

    cdef type t = type(obj)

    if t == tuple or t == list:
        target[0] = 0
        target[1] = 0
        if t == tuple:
            _hash_x64_128_iterable(obj, 1, target)
        elif t == list:
            _hash_x64_128_iterable(obj, 2, target)
        return

    try:
        if t == bytes:
            seed = 1
            PyBytes_AsStringAndSize(obj, &data, &size)
        elif t == str:
            seed = 2
            data = PyUnicode_AsUTF8AndSize(obj, &size)
        elif t == array.array:
            seed = 3
            data = (<array.array> obj).data.as_chars
            size = len(obj)
        elif t == bytearray:
            seed = 4
            data = <bytearray> obj
            size = len(obj)
        elif t == py_memoryview:
            seed = 5
            PyObject_GetBuffer(obj, &buf, buffer.PyBUF_SIMPLE)
            data = <char*> buf.buf
            size = buf.len
        elif t == int:
            seed = 6
            if MIN_LONG <= obj <= MAX_LONG:
                size = <Py_ssize_t> sizeof(long)
                data = <char*> malloc(size)
                if not data:
                    raise MemoryError()
                do_free = True
                ((<long*> data)[0]) = PyInt_AsLong(obj)
            else:
                obj = '%x' % obj
                data = PyUnicode_AsUTF8AndSize(obj, &size)
        elif t == float:
            seed = 7
            size = <Py_ssize_t> sizeof(double)
            data = <char*> malloc(size)
            if not data:
                raise MemoryError()
            do_free = True
            ((<double*> data)[0]) = PyFloat_AsDouble(obj)
        else:
            raise Exception('Unsupported type %s' % t)

        hash_x64_128_raw(data, size, seed, target)
    finally:
        if do_free:
            free(data)


cdef tuple _hash_x64_128_iterable(iterable, int seed, ulong* target):
    cdef ulong *tmp = <ulong*> malloc(2 * sizeof(ulong))
    if not tmp:
        raise MemoryError()

    cdef ulong h1 = seed
    cdef ulong h2 = seed

    try:
        for element in iterable:
            _hash_x64_128(element, tmp)
            h1 ^= tmp[0]
            h2 ^= tmp[1]

            h1 += h2
            h2 += h1

            h1 = fmix(h1)
            h2 = fmix(h2)

            h1 += h2
            h2 += h1

        target[0] = h1
        target[1] = h2
    finally:
        free(tmp)


cdef inline void hash_x64_128_raw(char* data, Py_ssize_t size, int seed, ulong* target):
    '''
    Computes a 128 bit murmer3 hash.

    Args:
        data: A pointer to the byte string to compute the hash on.
        size: The size of the byte string.

    Returns:
        A 2-tuple of the two 64 bit parts of the 128 bit hash.
    '''
    cdef int tail_size = size % 16
    cdef int body_size = size - tail_size
    cdef int nblocks = size // 16

    cdef ulong h1 = seed
    cdef ulong h2 = seed

    cdef ulong c1 = 0x87c37b91114253d5
    cdef ulong c2 = 0x4cf5ad432745937f

    cdef ulong* blocks = <ulong*> data

    cdef int i
    cdef ulong k1, k2

    for i in range(nblocks):
        k1 = blocks[i * 2 + 0]
        k2 = blocks[i * 2 + 1]

        k1 *= c1
        k1 = rotl64(k1, 31)
        k1 *= c2
        h1 ^= k1

        h1 = rotl64(h1, 27)
        h1 += h2
        h1 = h1 * 5 + 0x52dce729

        k2 *= c2
        k2 = rotl64(k2, 33)
        k2 *= c1
        h2 ^= k2

        h2 = rotl64(h2, 31)
        h2 += h1
        h2 = h2 * 5 + 0x38495ab5

    k1 = 0
    k2 = 0

    if tail_size > 8:
        if tail_size >= 15: k2 ^= get_block(data, body_size + 14) << 48
        if tail_size >= 14: k2 ^= get_block(data, body_size + 13) << 40
        if tail_size >= 13: k2 ^= get_block(data, body_size + 12) << 32
        if tail_size >= 12: k2 ^= get_block(data, body_size + 11) << 24
        if tail_size >= 11: k2 ^= get_block(data, body_size + 10) << 16
        if tail_size >= 10: k2 ^= get_block(data, body_size +  9) << 8
        if tail_size >=  9: k2 ^= get_block(data, body_size +  8)

        k2 *= c2
        k2 = rotl64(k2, 33)
        k2 *= c1
        h2 ^= k2

    if tail_size > 0:
        if tail_size >=  8: k1 ^= get_block(data, body_size +  7) << 56
        if tail_size >=  7: k1 ^= get_block(data, body_size +  6) << 48
        if tail_size >=  6: k1 ^= get_block(data, body_size +  5) << 40
        if tail_size >=  5: k1 ^= get_block(data, body_size +  4) << 32
        if tail_size >=  4: k1 ^= get_block(data, body_size +  3) << 24
        if tail_size >=  3: k1 ^= get_block(data, body_size +  2) << 16
        if tail_size >=  2: k1 ^= get_block(data, body_size +  1) << 8
        if tail_size >=  1: k1 ^= get_block(data, body_size     ) << 0

        k1 *= c1
        k1 = rotl64(k1, 31)
        k1 *= c2
        h1 ^= k1

    # finalization
    h1 ^= size
    h2 ^= size

    h1 += h2
    h2 += h1

    h1 = fmix(h1)
    h2 = fmix(h2)

    h1 += h2
    h2 += h1

    # return <long> h1, <long> h2
    target[0] = h1
    target[1] = h2
