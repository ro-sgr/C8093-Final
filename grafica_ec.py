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

_run = True
while _run:

    print("----- GRAFICADORA -----\n")

    texto_pausado("(1) Graficar tiempo real\n",0.025)
    texto_pausado("(2) Mover parámetros\n",0.025)
    texto_pausado("(3) Salir\n",0.025)

    opcion = input("> ")

    clear()

    if opcion == '1':

        TIEMPO = []
        THETA = []

        print("----- Graficar en tiempo real -----\n")

        k = entrada_real(input("Constante del resorte de torsión: "))
        A = entrada_real(input("Área de la espira: "))
        B = entrada_real(input("Intensidad de campo magnético: "))
        R = entrada_real(input("Resistencia de la espira: "))
        I = entrada_real(input("Momento de incercia de la espira: "))
        t = entrada_real(input("Tiempo: "))

        def calcular_thetadpt(ang, vel):
            ac_ang = (- k * ang )/I - (((A*B)**2)* ((np.cos(ang)) ** 2) * vel) / (R*I) # Cálculo de la aceleración angular según la ecuación diferencial
            return ac_ang

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

        # FIGURA 1

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.grid(linestyle = "--")

        f(t)
        cuadros = 1. / 20

        def animacion1(i):
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
            ax1.set_xlim(a,b) # Límites de la gráfica
            ax1.set_ylim(-5,5)
            linea, = ax1.plot(t, th, "r") 
            return linea

        animation1 = FuncAnimation(fig1, func=animacion1, frames=np.arange(0, 100, (1./20)),interval = 10) # Altérese el denominador del paso para cambiar los fps

        # FIGURA 2

        # Se definen dos listas que serán las entradas en X y en Y de los extremos del oscilador
        x = [0,0]
        y = [0,0]

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)

        # Límites de la gráfica
        ax2.set_xlim(-2,2)
        ax2.set_ylim(-2,2)

        # Recta que hará de oscilador
        linea, = ax2.plot(0, 0, "r")

        # Recta horizontal
        plt.axhline(0,-2,2,color="gray",linestyle="--")

        def animacion2(i):
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

        animation2 = FuncAnimation(fig2, func=animacion2, frames=np.arange(0, 100, (1./60)),interval = t)

        plt.show()

        clear()

    elif opcion == '2':

        from matplotlib.widgets import Slider, Button, RadioButtons

        print("----- Mover parámetros -----\n")

        temp = entrada_real(input("Tiempo: "))

        def f(k, A, B, R, I):
            TIEMPO =[]
            THETA = []
            # Valores iniciales
            theta = 1.
            thetapt = 0.
            dt = 0.01 # Paso

            for i in np.arange(0.0, temp, dt):
                TIEMPO.append(i)

                # Se calcula la aceleración angular según las variables
                thetadpt = (- k * theta )/I - (((A*B)**2)* ((np.cos(theta)) ** 2) * thetapt) / (R*I)

                thetapt += thetadpt * dt # Como thetapt = d/dt (thetapt), theta_f = thetapt_i + thetadpt * dt, aprox
                theta += thetapt * dt # Como thetapt = d/dt (theta), theta_f = thetapt_i + thetadpt * dt, aprox
                THETA.append(theta)

            return np.asarray(THETA)

        fig = plt.figure()
        ax = fig.add_subplot(111)

        fig.subplots_adjust(bottom=0.4)

        t = np.arange(0.0, temp, 0.01)
        k_0 = 1.
        A_0 = 1.
        B_0 = 1.
        R_0 = 1.
        I_0 = 1.

        # Gráfica inicial
        y = f(k_0, A_0, B_0, R_0, I_0)
        [linea] = ax.plot(t, y, linewidth=2, color='red')
        ax.set_xlim([0, temp])
        ax.set_ylim([-2, 2])

        # Add two sliders for tweaking the parameters

        # Define an axes area and draw a slider in it
        k_deslizador_ax  = fig.add_axes([0.1, 0.25, 0.65, 0.03])
        k_deslizador = Slider(k_deslizador_ax, 'k', 0.1, 10.0, valinit=k_0)

        # A
        A_deslizador_ax = fig.add_axes([0.1, 0.20, 0.65, 0.03])
        A_deslizador = Slider(A_deslizador_ax, 'A', 0.1, 10.0, valinit=A_0)

        # B
        B_deslizador_ax = fig.add_axes([0.1, 0.15, 0.65, 0.03])
        B_deslizador = Slider(B_deslizador_ax, 'B', 0.1, 10.0, valinit=B_0)

        # R
        R_deslizador_ax = fig.add_axes([0.1, 0.10, 0.65, 0.03])
        R_deslizador = Slider(R_deslizador_ax, 'R', 0.1, 10.0, valinit=R_0)

        # I
        I_deslizador_ax = fig.add_axes([0.1, 0.05, 0.65, 0.03])
        I_deslizador = Slider(I_deslizador_ax, 'I', 0.1, 10.0, valinit=I_0)

        # Modifica la línea cuando un valor de cualquier deslizador cambia
        def cambio_deslizadores(val):
            linea.set_ydata(f(k_deslizador.val, A_deslizador.val, B_deslizador.val, R_deslizador.val, I_deslizador.val))
            fig.canvas.draw_idle()

        k_deslizador.on_changed(cambio_deslizadores)
        A_deslizador.on_changed(cambio_deslizadores)
        B_deslizador.on_changed(cambio_deslizadores)
        R_deslizador.on_changed(cambio_deslizadores)
        I_deslizador.on_changed(cambio_deslizadores)

        # Botón para reiniciar los parámetros
        boton_reinicio_ax = fig.add_axes([0.85, 0.075, 0.1, 0.04])
        boton_reinicio = Button(boton_reinicio_ax, 'Reiniciar', hovercolor='0.975')
        def reinicio_boton(mouse_event):
            k_deslizador.reset()
            A_deslizador.reset()
            B_deslizador.reset()
            R_deslizador.reset()
            I_deslizador.reset()

        boton_reinicio.on_clicked(reinicio_boton)

        plt.show()

    elif opcion == '3':
        _run = False

    else:
        texto_pausado("Por favor, introduzca una opción válida\n", 0.05)
        texto_pausado(" ",0.8)
        clear()

clear()
exit()
