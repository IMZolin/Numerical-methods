from math import *
from matplotlib import pyplot as plt
import numpy as np

eps = 10**-6
b = exp(1)
y0 = 0.01105170918

def main_part(N,y0,  end=1.1):
    Error = []
    x0 = 0.1
    z0 = 0.23208589279
    h = (end - x0) / N
    array_x = []
    array_y = []

    while x0 < 1.1:
        k1 = h * function(x0, y0, z0)
        q1 = h * gfuncti0n(x0, y0, z0)

        k2 = h * function(x0 + h / 2, y0 + q1 / 2, z0 + k1 / 2)
        q2 = h * gfuncti0n(x0 + h / 2, y0 + q1 / 2, z0 + k1 / 2)

        k3 = h * function(x0 + h / 2, y0 - q1 + 2 * q2, z0 - k1 + 2 * k2)
        q3 = h * gfuncti0n(x0 + h / 2, y0 - q1 + 2 * q2, z0 - k1 + 2 * k2)

        z1 = z0 + (k1 + 4 * k2 + k3) / 6
        y1 = y0 + (q1 + 4 * q2 + q3) / 6

        z0 = z1
        y0 = y1
        array_y.append(y0)
        array_x.append(x0)

        x0 += h
    return array_y


def gfuncti0n(x, y, z):
    return z

def function(x, y, z):
    return (2+1/x)*z - (1+2/x)*y+x*exp(x)

def main():

    return
