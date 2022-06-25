def check(A, B, x, l):
    return all([((a * x) % q) >> l == b for a, b in zip(A, B)])


# module
q = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123

# bit length of module
m = 256

# number of equations
n = 80

# A (Coeff)
A = [
   
]

# B known MSB length
s = 6

# B known MSB (KnownNonce)
KnownNonce = [
    33, 32, 21, 36, 26, 35, 4, 41, 44, 16, 45, 3, 43, 36, 31, 15, 56, 20, 26, 10, 55, 49, 52,
    34, 39, 24, 15, 7, 21, 52, 56, 47, 51, 31, 5, 8, 32, 3, 23, 9, 34, 51, 49, 47, 0, 4, 11, 
    6, 50, 4, 47, 49, 18, 13, 13, 23, 37, 13, 58, 32, 11, 14, 13, 50, 57, 31, 2, 20, 12, 29, 
    52, 62, 21, 42, 48, 17, 55, 53, 28, 56 
]


# if 'DEBUG':
#     x = 10277847721795914431756115951552197892660363513639868708768066915171823347352452322827321894508504645532536737811666
#     KnownNonce = [((a * x) % q) >> (m - s) for a in A]

assert all([_.nbits() <= s for _ in KnownNonce])

# B unknown LSB length
l = m - s

B = [_ * 2^l for _ in KnownNonce]

assert len(A) == n
assert len(B) == n


w = 2^(l - 1)


A_ = [(a * inverse_mod(A[0], q)) % q for a in A]
B_ = [(b - a * B[0] + w) % q for a, b in zip(A_, B)]
#  + w * (a)
L = block_matrix(ZZ, [
    q * identity_matrix(ZZ, n - 1),
    Matrix(ZZ, A_[1:]),
    Matrix(ZZ, B_[1:])
],nrows=3, subdivide=False)

R = Matrix(ZZ, [
    *([[0, 0]] * (n - 1)),
    [    1,    0],
    [    0,    w]
])

M = block_matrix(ZZ, [L, R], ncols=2, subdivide=False)

MB = M.BKZ(block_size=15)

for vec in MB:
    if abs(vec[-1]) == w:

        print(vec.change_ring(RR) / w)
        
        if vec[-1] > 0:
            ks = [- vec[-2]] + [-(_ - w) for _ in vec[:-2]]
        else:
            ks = [vec[-2]] + [_ + w for _ in vec[:-2]]

        res = set()
        
        for a, b, k in zip(A, B, ks):
            if gcd(a, q) == 1:
                res.add((b + k) * inverse_mod(a, q) % q)
        
        assert len(res) == 1

        x = list(res)[0]

        if check(A, KnownNonce, x, l):
            print('found')
            break

print(x)
