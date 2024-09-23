from logical import Logical


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


def main(text):
    """
    j = 0;
        T1  = 0xbdbb52b0
            = 3183170224 = 0x5be0cd19 ( h ) + 0x1f5055e7 (sigma1(e) + Ch(e, f, g)) + 24 ( W ) + 0x428a2f98 (K)
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

    return [hex(i) for i in hash_values[-1]]


if __name__ == '__main__':
    hash_vals = main(text="abc")
    print(hash_vals)
