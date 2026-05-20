from QuadForm import QuadForm as Quad

def factor(target):

    # choose the appropriate discriminant
    disc = -target if -target % 4 == 1 else -4 * target

    for _ in range(100):
        
        # sample a random ideal
        sample = Quad.randomComposable(disc)

        order = sample.order()
        if (order % 2 == 1):
            continue
        ambiguous = sample.exp(order // 2)

        # extract factorisation information
        a, b, c = ambiguous.form()
        
        if b == 0:
            factors = [a, c]
        
        elif a == b:
            factors = [b//2, (2*c)-(b//2)] if b % 2 == 0 else [b, (4*c)-b]

        elif a == c:
            factors = [a-(b//2), a+(b//2)] if b % 2 == 0 else [(2*a)-b, (2*a)+b]

        if 1 in factors:
            continue

        # recursively factor - return concatenated list
        return factor(factors[0]) + factor(factors[1])
    
    # if no factorisation is found, assume with high probability that the target is prime
    return [target]

print(factor(2167583))