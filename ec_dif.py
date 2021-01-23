import numpy as np

# Parámetros constantes:
k = 1.
A = 1.
B = 1.
R = 1.
I = 1.

def calcular_thetadpt(ang, vel):
    ec_theta = (- k * ang )/I - ((A*B)**2)* (np.cos(ang)) ** 2 * vel / R 
    return ec_theta

THETA = []

def f(t):
    # Valores iniciales
    theta = 1.
    thetapt = 0.
    dt = 0.01
    iteraciones = 0
    for tiempo in np.arange(0,t,dt):
        iteraciones += 1
        thetadpt = calcular_thetadpt(theta, thetapt)
        thetapt = thetapt + thetadpt * dt
        theta = theta + thetapt * dt
        if iteraciones % 10 == 0:
            THETA.append(theta)

tiempo = int(input("Tiempo: "))
f(tiempo)
print(THETA)
