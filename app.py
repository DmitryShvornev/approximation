import numpy as np
from math import pi, cos
import matplotlib.pyplot as plt

K = np.array([[0.1951, -0.5556, 0.8315, -0.9808],
              [0.8315, -0.1951, -0.9808, -0.5556],
              [0.9808, 0.8315, 0.5556, 0.1951],
              [0.5556, -0.9808, 0.1951, 0.8315]])
b = np.array([[0.0549], [0.7385], [0.969], [0.3713]])
print(b)
a = np.dot(np.linalg.inv(K),b)
A = [a[0][0],a[1][0],a[2][0],a[3][0]]

def f(x):
    res = 0
    for n in range(1, 5):
        res += A[n - 1] * cos(pi*x*(2*n - 1) / 2)
    return res

F = [f(x) for x in np.linspace(-1, 1, 1001)]
F1 = [(x**2 -1)**2 for x in np.linspace(-1, 1, 1001)]
plt.plot(np.linspace(-1, 1, 1001), F)
plt.plot(np.linspace(-1, 1, 1001), F1)

# графики разностей не забыть
from sympy import *
x = symbols('x')
X = [-1, -5/8, -1/8, 3/8, 7/8, 1]
K = np.zeros((4,4))
b = np.zeros((4))
for n in range(1, 5):
    for m in range(0, 4):
        K[m][n-1] = integrate(cos(pi*x*(2*n-1)/2), (x, X[m], X[m+1]))
        b[m] = integrate((x**2 - 1)**2, (x, X[m], X[m+1]))
print(K)
print(b)