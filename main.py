import random
from functools import lru_cache
from sha256 import SHA256


class EC:
    def __init__(self, x=None, y=None):
        self.p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
        self.a = 0x0
        self.b = 0x7
        self.G_x = 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798
        self.G_y = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8
        self.n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
        self.h = 1
        self.x = x
        self.y = y

    @property
    def xy(self):
        if self.x and self.y:
            return hex(self.x), hex(self.y)
        else:
            return "points on the curve were not given"


@lru_cache
def add_points(P: EC, Q: EC) -> EC:
    assert isinstance(P, EC), f"P must be an Elliptic Curve with type {EC}, got {P} instead"
    assert isinstance(Q, EC), f"Q must be an Elliptic Curve with type {EC}, got {Q} instead"
    assert P or Q, f"Can not add two null points"

    if not P.x:
        return Q
    if not Q.x:
        return P

    if P.x != Q.x and P.y != Q.y:
        lam = ((Q.y - P.y) * pow(Q.x - P.x, -1, 2**32)) % 2**32
    else:
        lam = 3 * (P.x << 1) + P.a

    x_r = (lam << 1) - P.x - Q.x
    y_r = lam * (P.x - x_r) - P.y
    return EC(x_r, y_r)


def mult_points(P: EC, d) -> int | EC:
    if d == 0:
        return EC()
    elif d == 1:
        return P
    elif int(d % 2) == 1:
        return add_points(P, mult_points(P, d - 1))
    else:
        return mult_points(add_points(P, P), int(d / 2))


def ecdsa(message, sk):
    """
    curve = secp256k1 -> y^2 = x^3 + 7
    """
    concat_hex = ''.join(h[2:] for h in [hex(i) for i in message])
    z = int(concat_hex, 16)
    k = random.randint(1, EC(0, 0).n - 1)
    curve_points = mult_points(EC(EC().G_x, EC().G_y), k)
    r = curve_points.x % EC().n
    s = (k ** - 1) * (z + r * sk) % EC().n
    return r, s


def main():
    """
    implement a quick and dirty coin, nothing fancy, just the protocol

    digital_sign = given a msg, hash it and add
    """
    a_pk = (55, 3)
    a_sk = (5, 11)
    transaction_hash = SHA256(msg="A pays B 100").sha256()
    digital_sign = ecdsa(message=transaction_hash, sk=a_sk)

    return digital_sign


if __name__ == '__main__':
    hash_vals = SHA256(msg="abc").sha256()
    print([hex(i) for i in hash_vals])
    print(main())
