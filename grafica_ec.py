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
    
        print("----- GRAFICADORA -----\n")

        k = entrada_real(input("Constante del resorte de torsión: "))
        A = entrada_real(input("Área de la espira: "))
        B = entrada_real(input("Intensidad de campo magnético: "))
        R = entrada_real(input("Resistencia de la espira: "))
        I = entrada_real(input("Momento de incercia de la espira: "))
        t = entrada_real(input("Tiempo: "))

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

        k_inicial = 1.
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

        f(k_inicial, A_inicial, B_inicial, R, I, temp)
        angulo = np.asarray(THETA)

        x = np.arange(0, temp, dt) # valores en x
        y = (- k_inicial * angulo )/I - (((A_inicial*B_inicial)**2)* ((np.cos(angulo)) ** 2) * vel_inicial) / R # valores en y     

        l, = plt.plot(x, y, lw=1)

        # Deslizador k
        axk = plt.axes([0.25, .03, 0.50, 0.02])
        dk = Slider(axk, 'K', 0, 1, valinit=k_inicial)

        # Deslizador A
        axA = plt.axes([0.25, .001, 0.50, 0.02])
        dA = Slider(axA, 'A', 0, 1, valinit=A_inicial)

        # Deslizador B
        axB = plt.axes([0.25, .001, 0.50, 0.02])
        dB = Slider(axB, 'A', 0, 1, valinit=B_inicial)

        # Controles de animación
        is_manual = True
        intervalo = 100 # ms, tiempo entre cada cuadro de animación
        loop_len = 5.0 # second por iteración
        scale = intervalo / 1000 / loop_len
        
        # ax.plot(TIEMPO, THETA)

        ax.set(xlabel='tiempo (s)', ylabel='theta (θ)', title='Tiempo vs Theta')
        ax.grid(linestyle = "--")

        # k

        def update_k(val):
            # actualizar curva
            l.set_ydata((- val * angulo )/I - (((A_inicial*B_inicial)**2)* ((np.cos(angulo)) ** 2) * vel_inicial) / R)
            # redraw canvas while idle
            fig.canvas.draw_idle()

        def update_slider_k(val):
            global is_manual
            is_manual=True
            update_k(val)

        # A

        def update_A(val):
            # actualizar curva
            l.set_ydata((- k_inicial * angulo )/I - (((val*B)**2)* ((np.cos(angulo)) ** 2) * vel_inicial) / R)
            # redraw canvas while idle
            fig.canvas.draw_idle()

        def update_slider_A(val):
            global is_manual
            is_manual=True
            update_A(val)

        # A

        def update_B(val):
            # actualizar curva
            l.set_ydata((- k_inicial * angulo )/I - (((A_inicial*val)**2)* ((np.cos(angulo)) ** 2) * vel_inicial) / R)
            # redraw canvas while idle
            fig.canvas.draw_idle()

        def update_slider_B(val):
            global is_manual
            is_manual=True
            update_B(val)

        # Global

        def update_plot(num):
            global is_manual
            if is_manual:
                return l, # don't change

            val_k = (dk.val + scale) % dk.valmax
            dk.set_val(val_k)

            val_A = (dA.val + scale) % dA.valmax
            dA.set_val(val_A)

            val_B = (dB.val + scale) % dB.valmax
            dB.set_val(val_B)

            is_manual = False # the above line called update_slider, so we need to reset this
            return l,

        def on_click(event):
            (xm_k,ym_k),(xM_k,yM_k) = dk.label.clipbox.get_points()
            (xm_A,ym_A),(xM_A,yM_A) = dA.label.clipbox.get_points()
            # k
            if xm_k < event.x < xM_k and ym_k < event.y < yM_k:
                # Evento ocurre en el deslizador, ignorar ya que está siendo tratado por update_slider
                return
                
            # A
            elif xm_A < event.x < xM_A and ym_A < event.y < yM_A:
                # Evento ocurre en el deslizador, ignorar ya que está siendo tratado por update_slider
                return

            else:
                # user clicked somewhere else on canvas = unpause
                global is_manual
                is_manual=False

        dk.on_changed(update_slider_k)
        dA.on_changed(update_slider_A)
        dB.on_changed(update_slider_B)

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
