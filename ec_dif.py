import numpy as np

def calcular_thetadpt(ang, vel):
    return - k/I * ang - (A*B)**2 / R * (np.cos(ang)) ** 2 * vel

def f(t):
    # Valores iniciales
    theta = 1.
    thetapt = 0.
    dt = 0.01
    iteraciones = 0
    for tiempo in np.arange(0,t,dt):
        iteraciones = iteraciones + 1
        thetadpt = calcular_thetadpt(theta, thetapt)
        thetapt = thetapt + thetadpt * dt
        theta = theta + thetapt * dt
        if iteraciones % 10 == 0:
            THETA.append(theta)

# Parámetros constantes:
k = 1.
A = 1.
B = 1.
R = 1.
I = 1.

THETA = []

tiempo = int(input("Tiempo: "))
f(tiempo)
print(THETA)
