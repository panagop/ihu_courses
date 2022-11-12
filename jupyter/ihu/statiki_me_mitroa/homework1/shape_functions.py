import numpy as np
import sympy as sym
from sympy import Matrix


def f1(w: float) -> float:
    return 1 - w

def f2(w: float) -> float:
    return 1 - 3 * w**2 + 2 * w**3

def f3(w: float, x: float) -> float:
    return x * (1 - w)**2

def f4(w: float) -> float:
    return w

def f5(w: float) -> float:
    return 3*w**2 + 2*w**3

def f6(w: float, x: float) -> float:
    return x * (-w + w**2)


def F(L: float, x: float) -> Matrix:
    w = x/L
    return Matrix([[f1(w), 0, 0, f4(w), 0, 0],
                [0, f2(w), f3(w, x), 0, f5(w), f6(w, x)]])

def Y(F: Matrix, U:Matrix) -> Matrix:
    return F * U