from os import system, name
from time import sleep

# Texto se "escribe" en pantalla
def texto_pausado(texto,tiempo):
    for i in range(len(texto)):
        print(texto[i], sep='', end='', flush=True); sleep(tiempo)

# Borra el contenido de la consola
def clear():
    # windows
    if name == 'nt':
        _ = system('cls')
    # mac y linux
    else:
        _ = system('clear')

clear()

_run = True
while _run:

    import numpy as np
    import matplotlib.pyplot as plt

    # Verificador que las entradas sean números reales
    def entrada_real(numero):
        try:
            entrada = float(numero)
        except:
            print("Por favor introduzca una entrada válida.")
            exit()
        return entrada

    print("----- GRAFICADORA -----\n")

    texto_pausado("(1) Graficar tiempo real\n",0.025)
    texto_pausado("(2) Mover parámetros\n",0.025)
    texto_pausado("(3) Salir\n",0.025)

    opcion = input("> ")

    clear()

    if opcion == '1':

        from matplotlib.animation import FuncAnimation

        print("----- Graficar en tiempo real -----\n")

        k_0 = entrada_real(input("Constante del resorte de torsión: "))
        A_0 = entrada_real(input("Área de la espira: "))
        B_0 = entrada_real(input("Intensidad de campo magnético: "))
        R_0 = entrada_real(input("Resistencia de la espira: "))
        I_0 = entrada_real(input("Momento de incercia de la espira: "))
        max_t = entrada_real(input("Tiempo: "))

        # Solución a la ecuación diferencial en un intervalo de tiempo
        def f(k, A, B, R, I, t):
            TIEMPO = []
            THETA = []

            # Valores iniciales
            theta = 1.
            thetapt = 0.
            dt = 0.01 # Paso

            # Iteramos para anexar valores a nuestra lista en pasos diferentes de tiempo
            for i in np.arange(0,t,dt):
                TIEMPO.append(i) # Se adjunta el segundo a una lista

                # Cálculo de la aceleración angular según la ecuación diferencial
                thetadpt = (- k * theta )/I - (((A*B)**2)* ((np.cos(theta)) ** 2) * thetapt) / (R*I)

                thetapt += thetadpt * dt # Como thetapt = d/dt (thetapt), theta_f = thetapt_i + thetadpt * dt, aprox
                theta += thetapt * dt # Como thetapt = d/dt (theta), theta_f = thetapt_i + thetadpt * dt, aprox
                THETA.append(theta)

            return np.asarray(TIEMPO), np.asarray(THETA), theta

        ### FIGURA 1

        fig1 = plt.figure() # generación de la figura
        fig1.canvas.set_window_title('Oscilador') # título de la ventana
        fig1.suptitle('Oscilador armónico', fontsize=16)

        ax1 = fig1.add_subplot(111) # tamaño de la gráfica
        ax1.grid(linestyle = "--") # cuadrícula
        ax1.set_facecolor('#000000') # color del plano
        ax1.set_xlabel(r't ($s$)') # eje X
        ax1.set_ylabel(r'theta ($ \theta $)') # eje Y

        sol = f(k_0, A_0, B_0, R_0, I_0, max_t) # solución a la ecuación diferencial
        x_fig1 = sol[0] # valores de tiempo
        y_fig1 = sol[1] # valores en y

        # Generación de las posiciones en x y en y del oscilador
        def animacion1(i):
            j = int(i / (1./60))
            if j < 101:
                t = x_fig1[0:j]
                th = y_fig1[0:j]
                a = x_fig1[0]
                b = x_fig1[200]
            else:
                t = x_fig1[j-100:j]
                th = y_fig1[j-100:j]
                a = x_fig1[j-100]
                b = x_fig1[j+100]

            # Límites de la gráfica
            # Observe que los valores en x irán cambiando, por lo que la gráfica se desplazará
            ax1.set_xlim(a,b)
            ax1.set_ylim(-2,2)

            linea, = ax1.plot(t, th, "#05d2ed") # Generación de la curva
            
            return linea

        # Se genera la animación haciendo uso de la función "FuncAnimation" nativa de matplotlib
        # Altérese el denominador de 'cuadros' para cambiar los fps
        animation1 = FuncAnimation(fig1, func=animacion1, frames=np.arange(0, 100, (1./60)),interval = max_t)

        ### FIGURA 2

        # Se definen dos listas que serán las entradas en X y en Y de los extremos del oscilador
        x = [0,0]
        y = [0,0]

        fig2 = plt.figure() # generación de figura
        fig2.canvas.set_window_title('Espira') # título de la ventana emergente
        fig2.suptitle('Vista superior de la simulación', fontsize=16)

        ax2 = fig2.add_subplot(111) # tamaño de la gráfica
        ax2.set_facecolor('#000000') # color del plano
        ax2.set_xticks([]) # eje X
        ax2.set_yticks([]) # eje Y

        # Límites de la gráfica
        ax2.set_xlim(-2,2)
        ax2.set_ylim(-2,2)

        # Recta que generará el oscilador
        linea, = ax2.plot(0, 0, "#ede505")

        # Recta horizontal
        plt.axhline(0,-2,2,color="gray",linestyle="--")

        def animacion2(i):

            # Se definen las coordenadas de los extremos
            x1 = np.cos(f(k_0, A_0, B_0, R_0, I_0, i)[2])
            y1 = np.sin(f(k_0, A_0, B_0, R_0, I_0, i)[2])
            x2 = -x1
            y2 = -y1

            # Se alteran las entradas de las listas que definen los extremos de acuerdo a los valores recién calculados
            x[0] = x1
            y[0] = y1
            x[1] = x2
            y[1] = y2

            linea.set_xdata(x) # Se dan los extremos en x del oscilador
            linea.set_ydata(y) # Se dan los extremos en y del oscilador

            return linea

        animation2 = FuncAnimation(fig2, func=animacion2, frames=np.arange(0, 100, (1./60)), interval=max_t)

        plt.show()

        clear()

    elif opcion == '2':

        from matplotlib.widgets import Slider, Button, RadioButtons

        print("----- Mover parámetros -----\n")

        # Valores iniciales del oscilador
        # Note que la combinación de ciertos valores pueden provocar que el oscilador no sea graficado
        k_0 = 1.
        A_0 = 1.
        B_0 = 1.
        R_0 = 1.
        I_0 = 1.
        max_t = entrada_real(input("Tiempo: "))

        # Función que da solución a la ecuación diferencial
        def f(k, A, B, R, I):
            THETA = []

            # Valores iniciales
            # Modifique para apreciar diferentes comportamientos del oscilador
            theta = 1.
            thetapt = 0.
            dt = 0.001 # Paso

            for _ in np.arange(0.0, max_t, dt):

                # Se calcula la aceleración angular según las variables
                thetadpt = (- k * theta )/I - (((A*B)**2)* ((np.cos(theta)) ** 2) * thetapt) / (R*I)

                thetapt += thetadpt * dt # Como thetapt = d/dt (thetapt), theta_f = thetapt_i + thetadpt * dt, aprox
                theta += thetapt * dt # Como thetapt = d/dt (theta), theta_f = thetapt_i + thetadpt * dt, aprox
                THETA.append(theta)

            return np.asarray(THETA)

        fig = plt.figure() # generamos la figura
        fig.canvas.set_window_title('Oscilador') # título de la ventana
        fig.suptitle('Oscilador armónico', fontsize=16)

        ax = fig.add_subplot(111) # espaciado de la figura
        ax.set_facecolor('#000000') # color del plano
        ax.set_xlabel(r't ($s$)') # eje X
        ax.set_ylabel(r'theta ($ \theta $)') # eje Y

        fig.subplots_adjust(bottom=0.4) # espaciado de la gráfica

        t = np.arange(0.0, max_t, 0.001)
        y = f(k_0, A_0, B_0, R_0, I_0)

        # Lista de valores en y
        # Nos ayudará a actualizar los parámetros a futuro
        [linea] = ax.plot(t, y, linewidth=2, color='red')

        ax.set_xlim([0, max_t]) # rango eje x
        ax.set_ylim([-max(y)-1, max(y)+1]) # rango eje y

        ## DESLIZADORES
        # En esta sección, se encuentran los deslizadores que nos permitirán modificar la gráfica en tiempo real, la estructura es la siguiente:
        # parámetro_deslizador_ax = posicionamiento del objeto en la pantalla
        # parámetro_deslizador = Slider(posición, etiqueta, valor mínimo, valor máximo, valor inicial)

        # k
        k_deslizador_ax  = fig.add_axes([0.1, 0.25, 0.65, 0.03])
        k_deslizador = Slider(k_deslizador_ax, r'$k$', 0.1, 100, valinit=k_0)

        # A
        A_deslizador_ax = fig.add_axes([0.1, 0.20, 0.65, 0.03])
        A_deslizador = Slider(A_deslizador_ax, r'$A$', 0.1, 50, valinit=A_0)

        # B
        B_deslizador_ax = fig.add_axes([0.1, 0.15, 0.65, 0.03])
        B_deslizador = Slider(B_deslizador_ax, r'$B$', 0.1, 50, valinit=B_0)

        # R
        R_deslizador_ax = fig.add_axes([0.1, 0.10, 0.65, 0.03])
        R_deslizador = Slider(R_deslizador_ax, r'$R$', 0.1, 100, valinit=R_0)

        # I
        I_deslizador_ax = fig.add_axes([0.1, 0.05, 0.65, 0.03])
        I_deslizador = Slider(I_deslizador_ax, r'$I$', 0.1, 100, valinit=I_0)

        # Modifica la línea cuando un valor de cualquier deslizador cambia
        def cambio_deslizadores(val):
            linea.set_ydata(f(k_deslizador.val, A_deslizador.val, B_deslizador.val, R_deslizador.val, I_deslizador.val))
            fig.canvas.draw_idle()

        # La función se ejecuta si se persive algún cambio en el deslizador
        k_deslizador.on_changed(cambio_deslizadores)
        A_deslizador.on_changed(cambio_deslizadores)
        B_deslizador.on_changed(cambio_deslizadores)
        R_deslizador.on_changed(cambio_deslizadores)
        I_deslizador.on_changed(cambio_deslizadores)

        # Botón para reestablecer los valores a las condiciones iniciales por defecto establecidas
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
