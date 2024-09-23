class Logical:
    @staticmethod
    def _right_rot(x: int, n: int) -> int:
        return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF

    def Ch(self, x: int, y: int, z: int) -> int:
        return (x & y) ^ (~x & z)

    def Maj(self, x: int, y: int, z: int) -> int:
        return (x & y) ^ (x & z) ^ (y & z)

    def sigma0(self, x: int) -> int:
        return self._right_rot(x, 2) ^ self._right_rot(x, 13) ^ self._right_rot(x, 22)

    def sigma1(self, x: int) -> int:
        return self._right_rot(x, 6) ^ self._right_rot(x, 11) ^ self._right_rot(x, 25)

    def sig0(self, x: int) -> int:
        return self._right_rot(x, 7) ^ self._right_rot(x, 18) ^ (x >> 3)

    def sig1(self, x: int) -> int:
        return self._right_rot(x, 17) ^ self._right_rot(x, 19) ^ (x >> 10)
