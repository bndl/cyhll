from unittest import TestCase
# from hyperloglog.hll import HyperLogLog, get_alpha, get_rho
from cyhll.hyperloglog import HyperLogLog, get_alpha, get_rho
import math
import os
import pickle


class HyperLogLogTestCase(TestCase):
    def test_alpha(self):
        alpha = [get_alpha(b) for b in range(4, 10)]
        self.assertEqual(alpha, [0.673, 0.697, 0.709, 0.7152704932638152, 0.7182725932495458, 0.7197831133217303])

    def test_alpha_bad(self):
        self.assertRaises(ValueError, get_alpha, 1)
        self.assertRaises(ValueError, get_alpha, 17)

    def test_rho(self):
        self.assertEqual(get_rho(0, 32), 33)
        self.assertEqual(get_rho(1, 32), 32)
        self.assertEqual(get_rho(2, 32), 31)
        self.assertEqual(get_rho(3, 32), 31)
        self.assertEqual(get_rho(4, 32), 30)
        self.assertEqual(get_rho(5, 32), 30)
        self.assertEqual(get_rho(6, 32), 30)
        self.assertEqual(get_rho(7, 32), 30)
        self.assertEqual(get_rho(1 << 31, 32), 1)
        self.assertRaises(ValueError, get_rho, 1 << 32, 32)

    def test_init(self):
        s = HyperLogLog(0.05)
        self.assertEqual(s.p, 9)
        self.assertEqual(s.alpha, 0.7197831133217303)
        self.assertEqual(s.m, 512)
        self.assertEqual(len(s.M), 512)

    def test_calc_cardinality(self):
        clist = [1, 5, 10, 30, 60, 200, 1000, 10000, 60000]
        n = 30
        rel_err = 0.05

        for card in clist:
            s = 0.0
            for c in range(n):
                a = HyperLogLog(rel_err)

                for i in range(card):
                    a.add(os.urandom(20))

                s += a.card()

            z = (float(s) / n - card) / (rel_err * card / math.sqrt(n))
            self.assertLess(-3, z)
            self.assertGreater(3, z)

    def test_merge(self):
        a = HyperLogLog(0.05)
        b = HyperLogLog(0.05)
        c = HyperLogLog(0.05)

        for i in range(2):
            a.add(str(i))
            c.add(str(i))

        for i in range(2, 4):
            b.add(str(i))
            c.add(str(i))

        a.merge(b)

        self.assertNotEqual(a, b)
        self.assertNotEqual(b, c)
        self.assertEqual(a, c)


    def test_merge_err(self):
        a = HyperLogLog(0.05)
        b = HyperLogLog(0.01)

        self.assertRaises(ValueError, a.merge, b)

    def test_pickle(self):
        a = HyperLogLog(0.05)
        for x in range(100):
            a.add(str(x))
        b = pickle.loads(pickle.dumps(a))
        self.assertEqual(a.M, b.M)
        self.assertEqual(a.alpha, b.alpha)
        self.assertEqual(a.p, b.p)
        self.assertEqual(a.m, b.m)
