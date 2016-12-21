from unittest import TestCase
import string
from array import array
from cyhll import murmer3


class MurmerHash3TestCase(TestCase):
    def test_types(self):
        s = string.ascii_lowercase
        b = s.encode()

        primitives = [
            (71840294.89317356, 667535430529636161, 1509664973462398741),
            (0.9206813904848233, -7208461354467138051, 3769734253642353695),
            (7, 2410945784122383567, 8796919646904760335),
            (2 ** 63 - 1, 8149519964364859869, -8422929453225419506),
            (2 ** 64, -1682715467720934315, 4680060760960518367),
            (s, 2023876412346819495, 3209753511981878745),
            (b, 1660024533806982568, -4059853047773441421),
            (bytearray(b), -5723744278208638090, -3453631051498498437),
            (memoryview(b), 6857659530477715540, -7412019202647509479),
            (array('u', s), -2280497468652441540, -7946040245602166900),
            (array('b', b), -5366871570902108727, -3361522844783288491),
            (array('l', range(100)), -6924894080844287909, -4009853411618444015),
        ]

        composites = [
            (primitives, 8091227537723635227, -5183301716046233955),
            (tuple(primitives), 4770561775264928040, 5739399004458794637),
            ((primitives[:len(primitives) // 2], primitives[len(primitives) // 2:]),
             6869818053394741145, -8548565857465016348)
        ]

        for value, h1, h2 in primitives + composites:
            h = murmer3.hash_x64_128(value)
            self.assertEqual(h1, h[0])
            self.assertEqual(h2, h[1])
