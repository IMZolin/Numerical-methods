import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
import math
import pandas as pd


n = 10
a = 0.2
b = 1

A = np.zeros(4*n)
B = np.zeros(2*n)
C = np.zeros(n)

x1 = np.zeros(4*n)
x2 = np.zeros(2*n)
x3 = np.zeros(n)

B_Mistake = np.zeros(2*n-1)
C_Mistake = np.zeros(n)

STEP = np.zeros(14)
STEP_COORDINATE = np.zeros(14)





data = []
with open("CW_MKR.txt") as f:
    for line in f:
        data.append([float(x) for x in line.split()])
data = pd.DataFrame(data)
print(data)

for i in range(0, 4*n, 1):
    A[i] = data[i][0]
    x1[i] = a + (b-a)/(4*n - 1)*i

for i in range(0, 2*n-1, 1):
    B[i] = data[i][1]
    x2[i] = a + (b-a)/(2*n - 2)*i

for i in range(0, n, 1):
    C[i] = data[i][2]
    x3[i] = a + (b-a)/(n - 1)*i

for i in range(0, 2*n-1, 1):
    B_Mistake[i] = data[i][3]


for i in range(0, n, 1):
    C_Mistake[i] = data[i][4]

for i in range(0, 14, 1):
    STEP[i] = data[i][5]
    STEP_COORDINATE[i] = data[i][6]




print(B_Mistake)

''''------------------SECOND PART 5th------------------'''

n_1 = 9
a_1 = 0.2
b_1 = 1

A_1 = np.zeros(4*n_1+1)
B_1 = np.zeros(n_1+1)
C_1 = np.zeros(2*n_1+1)

B_mistake_1 = np.zeros(n_1+1)
C_mistake_1 = np.zeros(2*n_1+1)


x1_1 = np.zeros(4*n_1+1)
x2_1 = np.zeros(n_1+1)
x3_1 = np.zeros(2*n_1+1)

STEP_1 = np.zeros(14)
STEP_COORDINATE_1 = np.zeros(14)
xSTEP_1 = np.zeros(14)

print()
data_1 = []
with open("CW_RK.txt") as f:
    for line in f:
        data_1.append([float(x) for x in line.split()])
data_1 = pd.DataFrame(data_1)
print(data_1)

for i in range(0, n_1 + 1, 1):
    B_1[i] = data_1[i][1]
    B_mistake_1[i] = data_1[i][3]
    x2_1[i] = a_1 + (b_1-a_1)/n_1*i

for i in range(0, 4*n_1+1, 1):
    A_1[i] = data_1[i][0]
    x1_1[i] = a_1 + (b_1-a_1) / (4*n_1) * i

for i in range(0, 2*n_1+1, 1):
    C_1[i] = data_1[i][2]
    C_mistake_1[i] = data_1[i][4]
    x3_1[i] = a_1 + (b_1-a_1) / (2*n_1) * i

for i in range(0, 14, 1):
    STEP_1[i] = data_1[i][5]
    STEP_COORDINATE_1[i] = data_1[i][6]
    xSTEP_1[i] = (b-a)/(2**(i+1))


print(STEP_COORDINATE)


''''-------------------GRAPHICS-------------------'''

plt.subplot(1, 2, 1)
plt.title("Графики функций")
plt.plot(x1, A, label='Exact fun', color='r')
plt.plot(x2[:-1], B[:-1], label=f'h = {(b-a)/(2*n)} Краевая задача', color='b', marker='o', markersize=4)
plt.plot(x3, C, label=f'h = {(b-a)/(n)} Краевая задача', color='g', marker='o', markersize=4)
plt.plot(x2_1, B_1, label=f'h = {(b-a)/(n_1 + 1)} Задача Коши', color='c', marker='o', markersize=4)
plt.plot(x3_1, C_1, color='y', label=f'h = {(b-a)/(2*n_1 + 2)} Задача Коши', marker='o', markersize=4)
#plt.plot(px, py)
plt.legend()
plt.xlabel("x")
plt.ylabel("y")



plt.subplot(1, 2, 2)
plt.title("Фактические ошибки")
plt.plot(x2[:-1], B_Mistake, label=f'h = {(b-a)/(2*n)} Краевая задача', color='b', marker='o', markersize=4)
plt.plot(x3, C_Mistake, label=f'h = {(b-a)/(n)} Краевая задача', color='g', marker='o', markersize=4)
plt.plot(x2_1, B_mistake_1, color='c', marker='o', label=f'h = {(b-a)/(n_1 + 1)} Задача Коши', markersize=4)
plt.plot(x3_1, C_mistake_1, color='y', marker='o', label=f'h = {(b-a)/(2*n_1 + 2)} Задача Коши', markersize=4)
plt.xlabel("x")
plt.ylabel("mistake")
plt.legend()
plt.show()


plt.subplot(1, 2, 1)
plt.title("Погрешность - Шаг")
plt.loglog(xSTEP_1[1:], STEP[1:], label='Краевая задача', color='b')
plt.loglog(xSTEP_1, STEP_1, label='Задача Коши', color='c')
plt.loglog(xSTEP_1, xSTEP_1*xSTEP_1, label='h^2', color='r')
plt.loglog(xSTEP_1, xSTEP_1*xSTEP_1*xSTEP_1*xSTEP_1, label='h^4', color='y')
plt.legend()
plt.xlabel("Шаг")
plt.ylabel("Погрешность")

plt.subplot(1, 2, 2)
plt.title("Координата max погрешности (отрезок [0.2;1])")
plt.semilogx(xSTEP_1[1:], STEP_COORDINATE[1:], label='Краевая задача', color='b')
plt.semilogx(xSTEP_1[1:-1], STEP_COORDINATE_1[1:-1],label='Задача Коши', color='c')
plt.xlabel("Шаг")
plt.ylabel("Координата")
plt.legend()
plt.show()



''''-------------------GRAPHICS-2-------------------'''

plt.subplot(1, 2, 1)
plt.title("Задача Коши")
plt.plot(x1_1, A_1, label='Exact fun', color='r')
plt.plot(x2_1, B_1, label=f'h = {(b-a)/(n_1 + 1)}', color='b', marker='o', markersize=4)
plt.plot(x3_1, C_1, color='g', label=f'h = {(b-a)/(2*n_1 + 2)}', marker='o', markersize=4)
#plt.plot(px, py)
plt.legend()
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(1, 2, 2)
plt.title("Задача Коши")
plt.plot(x2_1, B_mistake_1, color='b', marker='o', markersize=4)
plt.plot(x3_1, C_mistake_1, color='g', marker='o', markersize=4)
plt.show()

plt.subplot(1, 2, 1)
plt.title("Задача Коши")
plt.loglog(xSTEP_1, STEP_1, label='Погрешность - Шаг', color='b')
plt.xlabel("шаг")
plt.ylabel("погрешность")
plt.legend()

plt.subplot(1, 2, 2)
plt.title("Задача Коши")
plt.semilogx(xSTEP_1[:-1], STEP_COORDINATE_1[:-1],label='Координата max погрешности (отрезок [0.2;1])', color='g')
plt.xlabel("шаг")
plt.ylabel("координата")
plt.legend()
plt.show()