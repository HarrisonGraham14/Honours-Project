import random
import math
from algorithms import *

class QuadForm:
    def __init__(self, a, b, c):
        self.a = int(a)
        self.b = int(b)
        self.c = int(c)
    
    def reduce(self): # assumes forms are positive definite

        # following Cohen Algorithm 5.4.2 (pg 238)
        while True:
            if self.b <= -self.a or self.b > self.a:
                r = self.b % (2 * self.a)
                r = r if r <= self.a else r - 2 * self.a
                q = (self.b - r) // (2 * self.a)
                
                self.c = self.c - (self.b + r) * q // 2
                self.b = r
            
            if self.a > self.c:
                self.a, self.b, self.c = self.c, -self.b, self.a
                continue
            break

        if self.a == self.c and self.b < 0:
            self.b = -self.b

    def __mul__(self, other):
        return QuadForm.compose(self, other)

    # TO DO: 
    @staticmethod
    def compose(X, Y): # assumes forms are primitive positive definite

        # following Cohen Algorithm 5.4.7 (pg 242)
        # step 1
        if X.a > Y.a:
            X, Y = Y, X
        s = (X.b + Y.b) // 2
        n = Y.b - s

        # step 2
        if Y.a % X.a == 0:
            y1 = 0
            d = X.a
        else:
            u, v, d = extendedEuclid(Y.a, X.a)
            y1 = u
        
        # step 3
        if s % d == 0:
            y2 = -1
            x2 = 0
            d1 = d
        else:
            u, v, d1 = extendedEuclid(s, d)
            x2 = u
            y2 = -v

        # step 4
        v1 = X.a // d1
        v2 = Y.a // d1
        r = (y1 * y2 * n - x2 * Y.c) % v1
        Z = QuadForm(v1 * v2, Y.b + 2 * v2 * r, (Y.c * d1 + r * (Y.b + v2 * r)) // v1)
        Z.reduce()
        return Z

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b and self.c == other.c

    @staticmethod
    def identity(disc):
        return QuadForm(1, 0, -disc // 4) if disc % 4 == 0 else QuadForm(1, 1, (1 - disc) // 4)

    @staticmethod
    def inverse(form):
        return QuadForm(form.a, -form.b, form.c)

    # randomly generates a valid (0 or 1 mod 4), nonzero discriminant in the given range
    @staticmethod
    def randomDisc(min, max = 0):
        if min > max:
            min, max = max, min
        positiveOptions = 0 if max <= 0 else ((max + 3) // 4) + ((max - 4) // 4) - ((min + 2) // 4) - ((min - 5) // 4)
        negativeOptions = 0 if min >= 0 else ((1 - min) // 4) + ((0 - min) // 4) - ((0 - max) // 4) - ((-1 - max) // 4)
        if positiveOptions + negativeOptions == 0:
            raise Exception(f"No valid discriminants on the range [{min}, {max}]")
        selection = random.randint(0, positiveOptions + negativeOptions - 1)
        if selection < negativeOptions:
            disc = -4 * ((selection + 2) // 2) + (1 if selection % 2 == 0 else 0)
        else:
            selection -= negativeOptions
            disc = 4 * ((selection + 1) // 2) + (1 if selection % 2 == 0 else 0)
        return disc

    # for testing with small discriminants
    @staticmethod
    def random(disc, ensurePositive = False):
        b = random.randint(0, abs(disc) // 2) * 2 + (disc % 2)
        if not ensurePositive:
            b = random.choice([b, -b])
        N = (b**2 - disc) // 4
        
        divisors = []
        for i in range(1, math.isqrt(abs(N)) + 1):
            if N % i == 0:
                divisors.append(i)
                if i * i != abs(N):
                    divisors.append(abs(N) // i)
        
        if N == 0:
            a, c = 0, random.randint(0, abs(disc))
            if random.choice([True, False]):
                a, c = c, a
        else:
            a = random.choice(divisors)
            c = N // a
        
        if (disc > 0 or not ensurePositive) and random.choice([True, False]):
            a, c = -a, -c
        return QuadForm(a, b, c)
    
    # gives a random reduced, positive, definite, primitive form
    @staticmethod
    def randomComposable(disc):
        if disc >= 0:
            raise Exception("expected negative discriminant")
        gcd = 0
        while gcd != 1: # TO DO: improve method for finding primative forms
            form = QuadForm.random(disc, ensurePositive=True)
            form.reduce()
            _, _, gcd = extendedEuclid(form.a, form.b)
            _, _, gcd = extendedEuclid(gcd, form.c)
        return form

    def __str__(self):
        return f"({self.a}, {self.b}, {self.c})"
    
    def form(self):
        return self.a, self.b, self.c
    
    def disc(self):
        return pow(self.b, 2) - 4 * self.a * self.c
    
    def __lt__(self, other):
        return self.a < other.a or (self.a == other.a and self.b < other.b)
    
    # TO DO: optimise with powers of 2
    def exp(self, exponent):
        power = self
        for i in range(1, exponent):
            power *= self
        return power
    
    # TO DO: improve using BSGS
    def order(self):
        exponent = 1
        power = self
        identity = QuadForm.identity(self.disc())
        while power != identity:
            power *= self
            exponent += 1
        return exponent