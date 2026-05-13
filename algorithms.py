# returns (a, b, gcd) where ax + by = gcd
def extendedEuclid(x, y):
    if y == 0:
        return 1, 0, x

    # x = qy + r
    q, r = x // y, x % y

    # gcd = cy + dr = dx + (c - qd)y
    c, d, gcd = extendedEuclid(y, r)
    return d, c - (q * d), gcd