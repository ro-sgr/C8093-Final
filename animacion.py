import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
## Desde aquí hasta el siguiente marcador es el código de ec_dif.py
k = 10. ## Rodrigo es ingeniero industrial
A = 1.
B = 1.
R = 1.
I = 1.

def calcular_thetadpt(ang, vel):
    ac_ang = (- k * ang )/I - (((A*B)**2)* ((np.cos(ang)) ** 2) * vel) / (R*I) # Cálculo de la aceleración angular según la ecuación diferencial
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

# Se definen dos listas que serán las entradas en X y en Y de los extremos del oscilador
x = [0,0]
y = [0,0]

fig, ax = plt.subplots()

# Límites de la gráfica
ax.set_xlim(-2, 2)
ax.set_ylim(-2,2)

# Recta que hará de oscilador
linea, = ax.plot(0, 0, "r")

# Recta horizontal
plt.axhline(0,-2,2,color="gray",linestyle="--")

def animacion(i):
    # Se definen las coordenadas de los extremos
    x1 = np.cos(f(i))
    y1 = np.sin(f(i))
    x2 = -np.cos(f(i))
    y2 = -np.sin(f(i))
    
    # Se alteran las entradas de las listas que definen los extremos de acuerdo a los valores recién calculados

    x[0] = x1
    y[0] = y1

    x[1] = x2
    y[1] = y2

    linea.set_xdata(x) # Se dan los extremos en x del oscilador
    linea.set_ydata(y) # Se dan los extremos en y del oscilador
    
    return linea

animation = FuncAnimation(fig, func=animacion, frames=np.arange(0, 100, (1./60)),interval = 10) # Altérese el denominador del paso para cambiar los fps
plt.show()
