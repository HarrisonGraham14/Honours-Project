import random
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

    @staticmethod
    def compose(X, Y): # assumes forms are primitive positive definite

        # following Cohen Algorithm 5.4.7 (pg 242)
        # step 1
        if X.a > Y.a:
            X, Y = Y, X
        s = (X.b + Y.b) // 2
        n = Y.b - s

        # step 2
        if Y.b % X.b == 0:
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
            u, v, d1 = extendedEuclid(Y.a, X.a)
            x2 = u
            y2 = -v

        # step 4
        v1 = X.a // d1
        v2 = Y.a // d1
        r = (y1 * y2 * n - x2 * Y.c) % v1
        Z = QuadForm(v1 * v2, Y.b + 2 * v2 * r, (Y.c * d1 + r * (Y.b + v2 + r)) // v1)
        Z.reduce()
        return Z

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b and self.c == other.c

    @staticmethod
    def identity(disc):
        return QuadForm(1, 0, disc // 4) if disc % 4 == 0 else QuadForm(1, 1, (1 - disc) // 4)

    @staticmethod
    def inverse(form):
        return QuadForm(form.a, -form.b, form.c)

    @staticmethod
    def random(disc):
        p = random.randint(1, 10000)
        q = random.randint(1, 10000)
        f = QuadForm(p, disc - 2 * q, (pow(disc, 2) - 3 * disc * q + 4 * pow(q, 2)) // (4 * p))
        f.reduce()
        return f

    def __str__(self):
        return f"({self.a}, {self.b}, {self.c})"
    
    def form(self):
        return self.a, self.b, self.c
    
    def disc(self):
        return pow(self.b, 2) - 4 * self.a * self.c
    
    def __lt__(self, other):
        return self.a < other.a or (self.a == other.a and self.b < other.b)