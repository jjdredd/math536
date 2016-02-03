#! /usr/bin/env python3

import numpy

r = r1 = r2 = 0

while r <= r1 + r2 + 0.5:
    A = numpy.random.random_integers(1, 10, (2, 2))
    B = numpy.random.random_integers(1, 10, (2, 2))

    r1 = numpy.amax([abs(x) for x in numpy.linalg.eigvals(A)])
    r2 = numpy.amax([abs(x) for x in numpy.linalg.eigvals(B)])
    r = numpy.amax([abs(x) for x in numpy.linalg.eigvals(A + B)])


print(r, r1 + r2, r1, r2)
print (A)
print (B)
