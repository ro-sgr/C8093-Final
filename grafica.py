import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
## Desde aquí hasta el siguiente marcador es el código de ec_dif.py
k = 100.
A = 1.
B = 1.
R = 1.
I = 1.

def calcular_thetadpt(ang, vel):
    ac_ang = (- k * ang )/I - (((A*B)**2)* ((np.cos(ang)) ** 2) * vel) / R # Cálculo de la aceleración angular según la ecuación diferencial
    return ac_ang

TIEMPO = []
THETA = []

def f(tiempo):
    # Valores iniciales
    theta = 1.
    thetapt = 0.
    dt = 0.01 # Paso
    for i in np.arange(0,tiempo,dt):
        TIEMPO.append(i) # Se adjunta el segundo a una lista
        thetadpt = calcular_thetadpt(theta, thetapt) # Se calcula la aceleración angular según las variables
        thetapt += thetadpt * dt # Como thetapt = d/dt (thetapt), theta_f = thetapt_i + thetadpt * dt, aprox
        theta += thetapt * dt # Como thetapt = d/dt (theta), theta_f = thetapt_i + thetadpt * dt, aprox
        THETA.append(theta)
    return theta # Le agregué esto
    
## Aquí termina el código de ec_dif.py

fig, ax = plt.subplots()

ax.grid(linestyle = "--")

f(10)
cuadros = 1. / 20

def animacion(i):
    j = int(i / (1./20))
    if j < 101:
        t = TIEMPO[0:j]
        th = THETA[0:j]
        a = TIEMPO[0]
        b = TIEMPO[200]
    else: 
        t = TIEMPO[j-100:j]
        th = THETA[j-100:j]
        a = TIEMPO[j-100]
        b = TIEMPO[j+100]
    ax.set_xlim(a,b) # Límites de la gráfica
    ax.set_ylim(-5,5)
    linea, = ax.plot(t, th, "r") 
    return linea

animation = FuncAnimation(fig, func=animacion, frames=np.arange(0, 100, (1./20)),interval = 10) # Altérese el denominador del paso para cambiar los fps

plt.show()
