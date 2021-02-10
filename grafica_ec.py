from os import system, name
from time import sleep

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def texto_pausado(texto,tiempo):
    for i in range(len(texto)):
        print(texto[i], sep='', end='', flush=True); sleep(tiempo)

def clear(): 
    # windows 
    if name == 'nt': 
        _ = system('cls') 
    # mac and linux
    else: 
        _ = system('clear')

clear()

def entrada_real(numero):
    try:
        entrada = float(numero)
    except:
        print("Por favor introduzca una entrada válida.")
        exit()
    return entrada


## Desde aquí hasta el siguiente marcador es el código de ec_dif.py
k = 100.
A = 1.
B = 1.
R = 2.
I = 1.

_run = True
while _run:

    TIEMPO = []
    THETA = []

    def calcular_thetadpt(k, A, B, R, I, t, ang, vel):
        ac_ang = (- k * ang )/I - (((A*B)**2)* ((np.cos(ang)) ** 2) * vel) / R # Cálculo de la aceleración angular según la ecuación diferencial
        return ac_ang

    def f(k, A, B, R, I, t):
        # Valores iniciales
        theta = 1.
        thetapt = 0.
        dt = 0.01 # Paso

        for i in np.arange(0,t,dt):
            TIEMPO.append(i) # Se adjunta el segundo a una lista

            thetadpt = calcular_thetadpt(k, A, B, R, I, t, theta, thetapt) # Se calcula la aceleración angular según las variables
            thetapt += thetadpt * dt # Como thetapt = d/dt (thetapt), theta_f = thetapt_i + thetadpt * dt, aprox
            theta += thetapt * dt # Como thetapt = d/dt (theta), theta_f = thetapt_i + thetadpt * dt, aprox
            THETA.append(theta)

        return theta # Le agregué esto
    
    texto_pausado("----- GRAFICADORA -----\n",0.025)

    texto_pausado("(1) Graficar\n",0.025)
    texto_pausado("(2) Salir\n",0.025)

    opcion = input("> ")

    clear()

    if opcion == '1':
        print("----- GRAFICADORA -----\n")
        k = entrada_real(input("Constante del resorte de torsión: "))
        A = entrada_real(input("Área de la espira: "))
        B = entrada_real(input("Intensidad de campo magnético: "))
        R = entrada_real(input("Resistencia de la espira: "))
        I = entrada_real(input("Momento de incercia de la espira: "))
        t = entrada_real(input("Tiempo: "))
        ## Aquí termina el código de ec_dif.py

        fig, ax = plt.subplots()

        ax.grid(linestyle = "--")

        f(k, A, B, R, I, t)
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

        clear()

    elif opcion == '2':
        _run = False

    else:
        texto_pausado("Por favor, introduzca una opción válida\n", 0.05)
        texto_pausado(" ",0.8)
        clear()

clear()
exit()
