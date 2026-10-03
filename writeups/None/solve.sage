#!/usr/bin/env sage
# 47CON CTF - "none"
# ECDSA P-256, el hint es el nibble bajo del nonce
# Uso: sage solve.sage challenge.json

import json, hashlib, sys

p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
n = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
b = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
E = EllipticCurve(GF(p), [-3, b])
G = E(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,
      0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)

path = sys.argv[1] if len(sys.argv) > 1 else 'challenge.json'
data = json.load(open(path))
Qpub = E(data['public_Q'][0], data['public_Q'][1])
sigs = data['signatures']

# k = t*d + u (mod n), con k = 16*e + hint -> e = (t*d + u - hint)/16
inv16 = int(pow(16, -1, n))
B = n // 16 + 1

raw, T, U = [], [], []
for sg in sigs:
    h = int.from_bytes(hashlib.sha256(bytes.fromhex(sg['m'])).digest(), 'big') % n
    si = int(pow(int(sg['s']), -1, n))
    t, u = si * int(sg['r']) % n, si * h % n
    raw.append((t, u, int(sg['hint'])))
    T.append(t * inv16 % n)
    U.append((u - int(sg['hint'])) * inv16 % n)

# retículo HNP
m = len(T)
W = 2**320
c = W // B

M = Matrix(ZZ, m + 2, m + 2)
for i in range(m):
    M[i, i] = n * c
    M[m, i] = T[i] * c
    M[m + 1, i] = ((U[i] - B // 2) % n) * c
M[m, m] = W // n
M[m + 1, m + 1] = W

inv = [int(pow(t, -1, n)) for t in T]

def hunt(R):
    # cualquier coordenada de cualquier fila puede ser un e_i valido
    seen = set()
    for row in R.rows():
        for i in range(m):
            for sign in (1, -1):
                e = sign * int(row[i]) // c + B // 2
                d = (e - U[i]) * inv[i] % n
                if d in seen or not 0 < d < n:
                    continue
                seen.add(d)
                # filtro rapido antes de multiplicar en la curva
                if all(((tt * d + uu) % n) % 16 == hh for tt, uu, hh in raw[:4]):
                    if d * G == Qpub:
                        return d
    return None

R = M.LLL()
d = hunt(R)
for blk in (30, 45, 55, 60, 62, 65, 68, 70):
    if d:
        break
    R = R.BKZ(block_size=blk, proof=False)
    d = hunt(R)
    print('[+] BKZ-%d -> %s' % (blk, d))

if d is None:
    sys.exit('no recuperada: prueba a subir el tamaño de bloque')

print('d =', d)
print('FLAG: 47CON{%s}' % format(d, 'x'))