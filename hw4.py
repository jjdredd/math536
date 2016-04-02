#! /usr/bin/env python3

import math
import numpy

Pe = 10
h = 0.01
k = 0.2

N = int(1/h)

def FillMat(k):
    A = numpy.zeros((N - 1, N - 1))
    for i in range(0, N - 2):
        A[i, i] = 1 + (2 * k) / (Pe * h * h)
        A[i, i + 1] = (1/2 - 1 / (Pe * h)) * k/h
        if i != 0:
            A[i, i - 1] = -(1/2 + 1 / (Pe * h)) * k/h

        A[N - 2, N - 2] = 1 + (2 * k) / (Pe * h * h)
        A[N - 2, N - 3] = -(2 * k) / (Pe * h * h)

    return A

for k in [0.05, 0.1, 0.2]:
    A = FillMat(k)
    b = numpy.zeros(N - 1)
    n = 0
    while n * k < 0.8:
        b[0] += (1/2 + 1 / (Pe * h)) * k/h
        b = numpy.linalg.solve(A, b)
        n += 1

    fout = open("hw4_sln_2_" + str(k) + ".txt", mode = "w")
    for i, x in enumerate(b):
        fout.write ("{}\t{}\n".format(i*h, x))
