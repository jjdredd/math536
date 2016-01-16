#! /usr/bin/env python3

import numpy

A = numpy.array([[1, 0.35, 0],
                 [0.35, 1, 0.35],
                 [0, 0.35, 1]])

print(A)

print("=========")

print(numpy.linalg.cholesky(A))
