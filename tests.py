import random
import numpy as np

DEBUG = False

def repeatTest(test, repetitions = 1000):
    print("- " + str(test.__name__), end=": ")
    for i in range(repetitions):
        try:
            test()
        except NotImplementedError:
            print("\033[33munimplemented\033[0m")
            return
        except Exception as e:
            print("\033[31mfailed: \033[31m" + str(e) + "\033[0m")
            if DEBUG:
                raise e
            return
    print("\033[32msuccess\033[0m")

print("\nRunning all tests\n")

# TO DO: put checks in place for functions to ensure good data is passed


#========================================================================================================================#
#                                                       Algorithms                                                       #
#========================================================================================================================#

from algorithms import *

MAX_TERM = 10000

def testExtendedEuclid():
    a = random.randint(-MAX_TERM, MAX_TERM)
    a = 1 if a == 0 else a
    b = random.randint(-MAX_TERM, MAX_TERM)
    b = 1 if b == 0 else b
    u, v, gcd = extendedEuclid(a, b)
    if a * u + b * v != gcd:
        raise Exception(f"Extended Euclidean algorithm fails for {a}, {b}")

print("Testing Algorithms:")
repeatTest(testExtendedEuclid)
print()

#========================================================================================================================#
#                                                     Quadratic Forms                                                    #
#========================================================================================================================#

from QuadForm import QuadForm as Quad

MAX_DISC = 1000
MAX_TERM = 10000

def testRandomDiscriminants():
    disc = Quad.randomDisc(-MAX_DISC, MAX_DISC)
    if disc % 4 > 1 or disc == 0:
        raise Exception(f"Generated impossible discriminant {disc}")

def testReduction():
    disc = Quad.randomDisc(-MAX_DISC) # only test definite forms
    form = Quad.random(disc, ensurePositive=True)

    if form.disc() != disc:
        raise Exception(f"Discriminant mismatch, {form} generated to have discriminant {disc}, computed to be {form.disc()}")

    reduced = form
    reduced.reduce()
    if reduced.disc() != disc:
        raise Exception(f"Discriminant ({disc}) is not preserved by reduction of form {form}")

    doubleReduced = reduced
    doubleReduced.reduce()
    if reduced != doubleReduced:
        raise Exception(f"Repeated reductions yields {form} -> {reduced} -> {doubleReduced}")
    
    NUM_STEPS = 10
    transMatrix = np.matrix([[1, 0], [0, 1]])
    for _ in range(10):
        transMatrix *= np.matrix(random.choice([[[1, 0], [0, 1]], [[0, -1], [1, 0]], [[-1, 0], [0, -1]], [[0, 1], [-1, 0]]])) # powers of S = 0, -1 ; 1, 0
        transMatrix *= np.matrix([[1, random.randint(1, MAX_TERM // NUM_STEPS) * random.choice([1, -1])], [0, 1]]) # powers of T = 1, 1 ; 0, 1
    formMatrix = np.matrix([[form.a, form.b / 2], [form.b / 2, form.c]])
    relatedMatrix = transMatrix.T * formMatrix * transMatrix
    related = Quad(round(relatedMatrix[0, 0]), round(relatedMatrix[0, 1]), round(relatedMatrix[1, 1]))
    reducedRelated = related
    reducedRelated .reduce()
    if reducedRelated != related:
        raise Exception(f"Reductions of equivalent forms yeild {form} -> {reduced}, {related} -> {reducedRelated}")
    
def testReductionEdgeCases():
    a = random.randint(1, MAX_TERM - 1)

    # test case |b| = a <= c
    f = Quad(a, random.choice([a, -a]), random.randint(a + 1, MAX_TERM))
    r = f
    r.reduce()
    if f.a != r.a or abs(f.b) != r.b or f.c != r.c:
        raise Exception(f"Edge case {f} incorrectly reduces to {r}")

    # test case |b| < a = c
    f = Quad(a, random.randint(0, a - 1) * random.choice([1, -1]), a)
    r = f
    r.reduce()
    if f.a != r.a or abs(f.b) != r.b or f.c != r.c:
        raise Exception(f"Edge case {f} incorrectly reduces to {r}")

def testIdentityComposition():
    disc = Quad.randomDisc(-MAX_DISC)
    f = Quad.randomComposable(disc)
    id = Quad.identity(disc)
    if (f * id != f) or (id * f != f):
        raise Exception(f"Identity composition fails on {f}")

def testInverseComposition():
    disc = Quad.randomDisc(-MAX_DISC)
    f = Quad.randomComposable(disc)
    id = Quad.identity(disc)
    inv = Quad.inverse(f)
    if f * inv != id or inv * f != id:
        raise Exception(f"Inverse composition fails on {f} * {inv}")

def testCompositionAssociativity():
    disc = Quad.randomDisc(-MAX_DISC)
    f = Quad.randomComposable(disc)
    g = Quad.randomComposable(disc)
    if f * g != g * f:
        raise Exception(f"Associativity fails on {f} * {g}")

def testCompositionCommutativity():
    disc = Quad.randomDisc(-MAX_DISC)
    f = Quad.randomComposable(disc)
    g = Quad.randomComposable(disc)
    h = Quad.randomComposable(disc)
    if f * (g * h) != (f * g) * h:
        raise Exception(f"Commutativity fails on {f} * {g} * {h}")

def testCubeComposition():
    disc = Quad.randomDisc(-MAX_DISC)
    f = Quad.randomComposable(disc)
    if f * (f * f) != (f * f) * f:
        raise Exception(f"Commutativity fails when cubing {f}")

def testKnownCompositions():
    raise NotImplementedError

def testOrderFunction():
    disc = Quad.randomDisc(-MAX_DISC)
    f = Quad.randomComposable(disc)
    order = f.order()
    if f.exp(order) != Quad.identity(disc):
        raise Exception(f"Order computation on {f} failed")

print("Testing quadratic forms:")
repeatTest(testRandomDiscriminants)
repeatTest(testReduction)
repeatTest(testReductionEdgeCases)
repeatTest(testIdentityComposition)
repeatTest(testInverseComposition)
repeatTest(testCompositionAssociativity)
repeatTest(testCompositionCommutativity)
repeatTest(testCubeComposition)
repeatTest(testKnownCompositions, repetitions=1)
repeatTest(testOrderFunction)
print()