#! /usr/bin/env python3

import numpy

# A = numpy.array([[1, 0.35, 0],
#                  [0.35, 1, 0.35],
#                  [0, 0.35, 1]])

# print(A)

# print("=========")

# print(numpy.linalg.cholesky(A))

A = numpy.ndarray((5, 5))

for i in range(0, 5):
    for j in range(0, 5):
        A[i, j] = 1.0/( (i + 1) + (j + 1) - 1)

print(A)

print("=========")

C = numpy.linalg.cholesky(A)

print(C)

print("=========")


b = numpy.array([5, 3.55, 2.81428571428571,
		 2.34642857142857, 2.01746031746032])


print(numpy.linalg.solve(A, b))

A[4, 0] = A[0, 4] = 0.20001

print(numpy.linalg.solve(A, b))
