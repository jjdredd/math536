#! /usr/bin/env python3

import numpy


N = 200
A = numpy.ndarray((N, N))

for i in range(0, N - 1):
    A[i, i] = 1;

for i in range(0, N - 1):
    A[i, i + 1] = -0.1

for i in range(0, N - 1):
    A[i + 1, i] = -0.1

A[N - 1, N - 1] = 1e+6

b = numpy.array([0.0 for i in range(0, N)])

b[0] = 0.9
b[N - 1] = 999999.9
for i in range(1, N - 1):
    b[i] = 0.8

print(b)

print (numpy.linalg.solve(A, b))

