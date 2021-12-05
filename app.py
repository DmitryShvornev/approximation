import numpy as np
from math import pi, cos
import matplotlib.pyplot as plt
from sympy import *

# задание отрезка интегрирования, точек сетки, приближаемой функции, базисных функций
x = symbols('x')
N= lambda x,n: cos(pi * x * (2 * n - 1) / 2)
psi = lambda x: (x**2 -1)**2
psi_0 = 0
bound_d = -1
bound_u = 1
X = [-1, -5/8, -1/8, 3/8, 7/8, 1]

# задание графиков
gridsize = (1, 3)
fig = plt.figure(figsize=(18, 16))
ax1 = plt.subplot2grid(gridsize, (0, 0))
ax2 = plt.subplot2grid(gridsize, (0, 1))
ax3 = plt.subplot2grid(gridsize, (0, 2))
ax1.set_title("Метод поточечной коллокации")
ax2.set_title("Метод коллокаций по областям")
ax3.set_title("Метод Галеркина")


fig1 = plt.figure(figsize=(18, 16))
ax1_1 = plt.subplot2grid(gridsize, (0, 0))
ax2_1 = plt.subplot2grid(gridsize, (0, 1))
ax3_1 = plt.subplot2grid(gridsize, (0, 2))
ax1_1.set_title("Метод поточечной коллокации")
ax2_1.set_title("Метод коллокаций по областям")
ax3_1.set_title("Метод Галеркина")

# Метод поточечной коллокации
print("Метод поточечной коллокации")
K = np.array([[0.1951, -0.5556, 0.8315, -0.9808],
              [0.8315, -0.1951, -0.9808, -0.5556],
              [0.9808, 0.8315, 0.5556, 0.1951],
              [0.5556, -0.9808, 0.1951, 0.8315]])
b = np.array([[0.0549], [0.7385], [0.969], [0.3713]])
print(K)
print(b)
a = np.dot(np.linalg.inv(K),b)
A = [a[0][0],a[1][0],a[2][0],a[3][0]]

# функция, выводящая  на график результат аппроксимации
# а также исходную функцию
def form_result(ax_i, ax_j, A):
    def f(x):
        res = 0
        for n in range(1, 5):
            res += A[n - 1] * N(x,n)
        return res
    F = np.array([f(x) for x in np.linspace(bound_d, bound_u, 1001)])
    F1 = np.array([psi(x) for x in np.linspace(bound_d, bound_u, 1001)])
    ax_i.plot(np.linspace(bound_d, bound_u, 1001), F)
    ax_i.plot(np.linspace(bound_d, bound_u, 1001), F1)
    ax_j.plot(np.linspace(bound_d, bound_u, 1001), F1-F)

form_result(ax1,ax1_1,A)

# Метод коллокаций по областям
print("Метод коллокаций по областям")
K = np.zeros((4,4))
b = np.zeros((4))
for n in range(1, 5):
    for m in range(0, 4):
        K[m][n-1] = integrate(N(x,n), (x, X[m], X[m+1]))
        b[m] = integrate(psi(x)-psi_0 , (x, X[m], X[m+1]))
print(K)
print(b)
a = np.dot(np.linalg.inv(K),b)
print(a)
A = [a[0],a[1],a[2],a[3]]
form_result(ax2,ax2_1, A)

# Метод Галеркина
print("Метод Галеркина")
K = np.zeros((4,4))
b = np.zeros((4))
for n in range(1, 5):
    for m in range(1, 5):
        K[m-1][n-1] = integrate(N(x,n)*N(x,m), (x, bound_d,bound_u))
        b[m-1] = integrate((psi(x)-psi_0)*N(x,m), (x,bound_d, bound_u))
print(K)
print(b)
a = np.dot(np.linalg.inv(K),b)
print(a)
A = [a[0],a[1],a[2],a[3]]
form_result(ax3,ax3_1, A)

plt.show()

