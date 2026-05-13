import numpy as np
import bisect
from QuadForm import QuadForm as Quad


# returns a power of the element with order 2 given an upper bound on the element's order
def BSGS(element, orderBound):
    smallBound = int(np.ceil(np.sqrt(orderBound)))
    
    babySteps = [ Quad.identity(element.disc()) ]
    for i in range(1, smallBound):
        babySteps.append(babySteps[-1] * element)

        if babySteps[-1] == Quad.identity(element.disc()):
            if len(babySteps) % 2 == 0:
                return babySteps[len(babySteps) // 2]
            raise Exception("Element order not divisible by 2")
    
    sortedBabySteps = [(i, babySteps[i]) for i in range(len(babySteps))]
    sortedBabySteps.sort(key = lambda x: x[1])

    giantSteps = [ Quad.inverse(babySteps[-1] * element) ]
    while len(giantSteps) <= smallBound:
        sortedIndex = bisect.bisect_left(sortedBabySteps, giantSteps[-1], key = lambda x : x[1])
        
        if sortedIndex < len(sortedBabySteps) and sortedBabySteps[sortedIndex][1] == giantSteps[-1]:
            elementOrder = len(giantSteps) * smallBound + sortedBabySteps[sortedIndex][0]
            if elementOrder % 2 == 0:
                return babySteps[(elementOrder // 2) % smallBound] * giantSteps[(elementOrder // 2) // smallBound]
            raise Exception("Element order not divisible by 2")
        
        giantSteps.append(giantSteps[-1] * giantSteps[0])     
    raise Exception("Bad value of orderBound")


def factor(target):

    # choose the appropriate discriminant
    disc = -target if -target % 4 == 1 else -4 * target

    for i in range(100):
        
        # sample a random ideal
        sample = Quad.random(disc)

        # use Shanks' BSGS algorithm to find an element of order 2
        try:
            ambiguous = BSGS(sample, -disc)

        except Exception as ex:
            if ex.args[0] == "Element order not divisible by 2":
                continue
            raise

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

factor(214)