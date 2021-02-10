from os import system, name
from time import sleep

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider

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

    print("----- GRAFICADORA -----\n")

    texto_pausado("(1) Graficar tiempo real\n",0.025)
    texto_pausado("(2) Mover parámetros\n",0.025)
    texto_pausado("(3) Salir\n",0.025)

    opcion = input("> ")

    clear()

    if opcion == '1':

        TIEMPO = []
        THETA = []

        print("----- GRAFICADORA -----\n")

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
        THETA = []

        def calcular_thetadpt(k, A, B, R, I, ang, vel):
            ac_ang = (- k * ang )/I - (((A*B)**2)* ((np.cos(ang)) ** 2) * vel) / R # Cálculo de la aceleración angular según la ecuación diferencial
            return ac_ang

        def f(k, A, B, R, I, t):
            # Valores iniciales
            theta = 1.
            thetapt = 0.
            dt = 0.01 # Paso

            for i in np.arange(0,t,dt):
                thetadpt = calcular_thetadpt(k, A, B, R, I, theta, thetapt) # Se calcula la aceleración angular según las variables
                thetapt += thetadpt * dt # Como thetapt = d/dt (thetapt), theta_f = thetapt_i + thetadpt * dt, aprox
                theta += thetapt * dt # Como thetapt = d/dt (theta), theta_f = thetapt_i + thetadpt * dt, aprox
                THETA.append(theta)

            return theta # Le agregué esto
        
        print("----- Mover parámetros -----\n")

        k_0 = 1.
        A_inicial = 1.
        B_inicial = 1.
        R = 2.
        I = 1.

        vel_inicial = 0.
        theta_inicial = 1.

        temp = entrada_real(input("Tiempo: "))

        theta = 1.
        thetapt = 0.
        dt = 0.01 # Paso

        fig, ax = plt.subplots()

        f(k_0, A_inicial, B_inicial, R, I, temp)
        angulo = np.asarray(THETA)

        x = np.arange(0, temp, dt) # valores en x
        y = (- k_0 * angulo )/I - (((A_inicial*B_inicial)**2)* ((np.cos(angulo)) ** 2) * vel_inicial) / R # valores en y     

        l, = plt.plot(x, y, lw=1)

        # Deslizador k
        axk = plt.axes([0.25, .03, 0.50, 0.02])
        dk = Slider(axk, 'K', 0, 1, valinit=k_0)

        # Deslizador A
        axA = plt.axes([0.25, .03, 0.50, 0.02])
        dA = Slider(axk, 'K', 0, 1, valinit=k_0)

        # Controles de animación
        is_manual = True
        intervalo = 100 # ms, tiempo entre cada cuadro de animación
        loop_len = 5.0 # segundo por iteración
        scale = intervalo / 1000 / loop_len

        ax.set(xlabel='tiempo (s)', ylabel='theta (θ)', title='Tiempo vs Theta')
        ax.grid(linestyle = "--")

        # k

        def update_k(val):
            # actualizar curva
            l.set_ydata((- val * angulo )/I - (((A_inicial*B_inicial)**2)* ((np.cos(angulo)) ** 2) * vel_inicial) / R)
            # volver a dibujar el lienzo mientras está inactivo
            fig.canvas.draw_idle()

        def update_slider_k(val):
            global is_manual
            is_manual=True
            update_k(val)

        # Global

        def update_plot(num):
            global is_manual
            if is_manual:
                return l, # no modificar

            val_k = (dk.val + scale) % dk.valmax
            dk.set_val(val_k)

            is_manual = False # la línea anterior llamada update_slider, por lo que necesitamos restablecer esto
            return l,

        def on_click(event):
            (xm_k,ym_k),(xM_k,yM_k) = dk.label.clipbox.get_points()

            # k
            if xm_k < event.x < xM_k and ym_k < event.y < yM_k:
                # Evento ocurre en el deslizador, ignorar ya que está siendo tratado por update_slider
                return

            else:
                # se hizo clic en otro lugar del lienzo, i.e. reanudar la pausa
                global is_manual
                is_manual=False

        dk.on_changed(update_slider_k)

        fig.canvas.mpl_connect('button_press_event', on_click)

        animacion = FuncAnimation(fig, update_plot, interval=intervalo)

        plt.show()

        clear()

    elif opcion == '3':
        _run = False

    else:
        texto_pausado("Por favor, introduzca una opción válida\n", 0.05)
        texto_pausado(" ",0.8)
        clear()

clear()
exit()
