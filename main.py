import random
from functools import lru_cache
from logical import Logical
from decimal import Decimal


def pre_processing(msg: str) -> list[int]:
    msg_ascii = msg.encode(encoding='ascii', errors='ignore')

    zero_bits = (448 - (len(msg_ascii) * 8) - 1) % 512
    bin_msg = int.from_bytes(msg_ascii, 'big')
    padded_msg = ((((bin_msg << 1) + 1) << zero_bits) << 64) + len(msg_ascii) * 8

    blocks = []
    i = 0
    while i <= padded_msg.bit_length():
        n1 = padded_msg >> i & (2 ** 512 - 1)
        i += 512
        blocks.append(n1)

    return blocks


@lru_cache
def expanded_msg_blocks(block, index_j):
    if index_j <= 15:
        if index_j == 0 and (block.bit_length() < 512):
            return block >> (512 - 32)
        try:
            return block >> (block.bit_length() - 32 * (index_j + 1)) & (2 ** 32 - 1)
        except ValueError:
            return block & ((1 << 32) - 1)
    else:
        return ((Logical().sig1(expanded_msg_blocks(block, index_j - 2)) +
                 Logical().sig0(expanded_msg_blocks(block, index_j - 15)))
                + expanded_msg_blocks(block, index_j - 16)) % 2 ** 32


def sha256(text):
    """
    """
    blocks = pre_processing(text)

    hash_values = [[0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]]
    constant_words = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    ]

    lg = Logical()

    for i in range(1, len(blocks) + 1):
        a, b, c, d, e, f, g, h = hash_values[i - 1]

        for j in range(63 + 1):
            W = expanded_msg_blocks(blocks[i - 1], j)
            T1 = (h + lg.sigma1(e) + lg.Ch(e, f, g) + constant_words[j] + W) % (2 ** 32)
            T2 = (lg.sigma0(a) + lg.Maj(a, b, c)) % (2 ** 32)
            h = g
            g = f
            f = e
            e = d + T1 % (2 ** 32)
            d = c
            c = b
            b = a
            a = T1 + T2 % (2 ** 32)

        hash_values_i = [0] * 8

        for j, var in enumerate([a, b, c, d, e, f, g, h]):
            hash_values_i[j] = var + hash_values[i - 1][j]
        hash_values.append(hash_values_i)

    return hash_values[-1]


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
        return hex(self.x), hex(self.y)


@lru_cache
def add_points(P: EC, Q: EC) -> EC:
    assert isinstance(P, EC), f"P must be an Elliptic Curve with type {EC}, got {P} instead"
    assert isinstance(Q, EC), f"Q must be an Elliptic Curve with type {EC}, got {Q} instead"
    assert P or Q, f"Can not add two null points"

    if not P:
        return Q
    if not Q:
        return P

    if P.x != Q.x and P.y != Q.y:
        lam = (Q.y - P.y) / (Q.x - P.x)
    elif P.x == Q.x:
        return EC()
    else:
        lam = 3 * (P.x << 1) + P.a

    x_r = (lam << 1) - P.x - Q.x
    y_r = lam * (P.x - x_r) - P.y
    return EC(x_r, y_r)


def mult_points(P: EC, d):
    if d == 0:
        return 0
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
    s = k ** - 1 * (z + r * sk) % EC().n
    return r, s


def main():
    """
    implement a quick and dirty coin, nothing fancy, just the protocol

    digital_sign = given a msg, hash it and add
    """
    a_pk = (55, 3)
    a_sk = (5, 11)
    transaction_hash = sha256(text="A pays B 100")
    digital_sign = ecdsa(message=transaction_hash, sk=a_sk)

    return digital_sign


if __name__ == '__main__':
    hash_vals = sha256(text="A pays B 100")
    print([hex(i) for i in hash_vals])
    print(main())
