from os import system, name, path, kill, getppid
from time import sleep

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

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

        f(k, A, B, R, I, t)

        # Data for plotting
        x = TIEMPO
        y = THETA

        fig, ax = plt.subplots()
        ax.plot(x, y)

        ax.set(xlabel='Tiempo', ylabel='Ángulo',
            title='Gráfica Ec. Diferencial')
        ax.grid()

        plt.show()

        clear()
    
    elif opcion == '2':
        _run = False

    else:
        texto_pausado("Por favor, introduzca una opción válida\n", 0.05)

clear()
exit()
