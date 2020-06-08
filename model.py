from math import pi, atan2, sin, cos, sqrt, radians

import numpy as np
import cmath
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from scipy.special import hankel2

## Datos geometricos y de pesos

m = 0.367  # [kg]
g = 9.8  # [m/s]
rho = 1.225  # [kg/m^3]
S = 0.324  # [m^2]
b = 1.2  # [m]
c = S / b  # [m]
AR = b ** 2 / S  # [ND]
S_t = 0.09  # [m^2]
b_t = 0.46  # [m]
AR_t = b_t ** 2 / S_t  # [ND]
k_aero = 2.5
c_e = 5
A = radians(26.5)

# xcg = 149.137e-3 #[m]
xcg = 119.137e-3  # [m]
zcg = 5.814e-3  # [m]

xw = 0.09
lw = xcg - xw
zw = -0.05
hw = zcg - zw

xt = 0.570
lt = xcg - xt
zt = 15e-3
ht = zcg - zt

Iyy = 0.008  # [kg*m^2]

### Datos aerodinamicos

CL_alpha = 2 * pi * AR / (AR + 2)
CL_alpha_t = pi / 2 * AR_t
Li = 0.0051
CD_0 = 0.018
CD_0t = 0.021
k = 1 / pi / AR
k_t = 1 / pi / AR_t
# delta_t = -3*pi/180
CT = 0
depsilon = 0.2
alpha_w = 0
dalpha_w = 0
ddelta_t = 0  ## Variable que cambia con la cola, estudiar!!

h0 = 0.4 * sin(10 * pi / 180)  # [m]
a0 = 0  # [rad]
a = -0.5  # [ND]
phi = 0  # rad
dalpha = 0

### Parametros adimensionales

Uc = sqrt(2 * m * g / (rho * S))
Lc = c / 2
tc = Lc / Uc
Lambda = S_t / S
L = lt / lw
H = ht / hw
RHL = hw / lw
M = 2 * m / rho / S / c
X = rho / 8 * S * c ** 2 * lw / Iyy
h0ad = h0 / Lc

def ode_dydt(delta_t, fa):
    '''
    Given:
        :param delta_t: tail angle
        :param fa: frequency value

    Return the corresponding ODE model
    '''

    def model(Y, t, dalpha = 0):
        '''
            :param Y: initial state
            :param t: a set of time steps value
            :param dalpha: initial model variable
        '''

        u = Y[0]
        w = Y[1]
        q = Y[2]
        theta = Y[3]
        x = Y[4]
        h = Y[5]

        alpha = atan2(w, u)
        Ub = sqrt(u ** 2 + w ** 2)

        if not fa:
            if alpha < radians(10):
                CLs = CL_alpha * (alpha + alpha_w)
                CLns = CL_alpha * (1.5 * (dalpha + dalpha_w) - 2 * lw / c * q) / Ub
                CL = CLs + CLns
            else:
                CLs = CL_alpha * radians(10)
                CLns = CL_alpha * (0.5 * (dalpha + dalpha_w)) / Ub
                CL = CLs + CLns

            depsilon = 0.2
            CT = 0
            CDi = k * CLs ** 2

        else:
            CL, CT, _, _ = AleteoBAad(fa, Ub, h0ad, a0, phi, a, alpha, AR, 0, 0, t)
            depsilon = 0
            CDi = k * CL ** 2
            CLs = CL

        if alpha + delta_t - depsilon * CLs / CL_alpha < radians(25):
            CLts = CL_alpha_t * (alpha + delta_t - depsilon * CLs / CL_alpha)
            if fa:
                CLtns = CL_alpha_t * (1.5 * (dalpha + ddelta_t) - lt / Lc * q) / Ub
            else:
                CLtns = CL_alpha_t * (1.5 * (dalpha + ddelta_t) - 2 * lt / c * q) / Ub
            CLt = CLts + CLtns
        else:
            CLts = CL_alpha_t * radians(25)
            CLtns = CL_alpha_t * (0.5 * (dalpha + ddelta_t)) / Ub
            CLt = CLts + CLtns

        CDit = k_t * CLts ** 2

        Fxw = sin(alpha) * CL - cos(alpha) * (CDi + CD_0 - CT)
        Fxt = Lambda * (sin(alpha) * CLt - cos(alpha) * (CDit + CD_0t))
        Fxb = -cos(alpha) * Li
        Fzw = -cos(alpha) * CL - sin(alpha) * (CDi + CD_0 - CT)
        Fzt = Lambda * (-cos(alpha) * CLt - sin(alpha) * (CDit + CD_0t))
        Fzb = -sin(alpha) * Li

        du = -q * w + 1 / 2 / M * ((u ** 2 + w ** 2) * (Fxw + Fxt + Fxb) - sin(theta))
        dw = q * u + 1 / 2 / M * ((u ** 2 + w ** 2) * (Fzw + Fzt + Fzb) + cos(theta))
        dq = X * (u ** 2 + w ** 2) * (-Fzw - L * Fzt - RHL * (Fxw + H * Fxt))
        dtheta = q
        dx = u * cos(theta) + w * sin(theta)
        dh = -w * cos(theta) + u * sin(theta)

        dalpha = (u * dw - w * du) / Ub ** 2

        return [du, dw, dq, dtheta, dx, dh]

    return model

def AleteoBAad(f, Ub, h0, a0, phi, a, alpha_m, AR, t0, phase0, t):
    '''
    :param f: Frecuencia
    :param h0: Amplitud de heaving
    :param Ub: Velocidad
    :param a0: Amplitud de pitching
    :param phi: Desfase entre heaving y pitching
    :param a: Punto sobre el que hace pitching (perfil entre -1 y 1)
    :param alpha_m: angulo de ataque medio
    :param AR: Densidad
    :param t: Instante de tiempo
    :param phase0:
    :param t0: tiempo inicial
    '''
    omega = 2 * pi * f
    k = omega / Ub

    zphi = complex(0, phi)
    alpha0 = a0 * cmath.exp(zphi)
    phase = omega * (t - t0) + phase0

    zomega = complex(0, omega)
    zphase = complex(0, phase)
    G0 = 2 * pi * (Ub * alpha0 - zomega * h0 - zomega * (a - 1 / 2) * alpha0)

    zbessel = complex(hankel2(1, k), hankel2(0, k))
    C = hankel2(1, k) / zbessel
    F = C.real
    G = C.imag

    CL = 2 * pi * ((alpha_m + k * h0 * (G * cos(phase) + F * sin(phase)) + a0 * (
            cos(phase + phi) * (F + G * k * (a - 0.5)) + sin(phase + phi) *
            (-G + F * k * (a - 0.5)))) * AR / (AR + 2) + (
                           k * h0 * k / 2. * cos(phase) + a0 * cos(phase + phi) * a * k ** 2. / 2 - a0 * sin(
                       phase + phi) * k / 2) * AR / (AR + 1))
    CLM = 2 * pi * alpha_m * AR / (AR + 2)

    zk = complex(0, k)
    zbessel = complex(hankel2(1, k), hankel2(0, k))
    C1 = 1 / k * cmath.exp(-zk) / zbessel
    F1 = C1.real
    G1 = C1.imag

    cT1 = pi * a0 * k * sin(phase + phi) * (k * h0 * sin(phase) + a * k * a0 * sin(phase + phi) + a0 * cos(phase + phi))
    cT2 = 4 * (a0 * (cos(phase + phi) + k * (a - 0.5) * sin(phase + phi)) + k * h0 * sin(phase)) * (
            F1 * (k * h0 * cos(phase) +
                  a0 * (-2. * sin(phase + phi) + k * (a - 1) * cos(phase + phi))) + G1 * (-k * h0 * sin(phase) +
                                                                                          a0 * (-2. * cos(
                phase + phi) - k * (a - 1) * sin(phase + phi))) - pi / 2 * a0 * (
                    F * cos(phase + phi) - G * sin(phase + phi)))

    zphase = complex(0, phase)
    zalpha = complex(alpha_m, alpha0 * cmath.exp(zphase))
    CT = -zalpha.real * CL + cT1 * AR / (AR + 1) + cT2 * AR / (AR + 2)
    CTM = -alpha_m * CLM + AR / (AR + 2) * (-2 * (k * h0) ** 2 * G1 - a0 ** 2 * (
            2 * pi * F + 2 * a * k * F1 + 2 * (2 + k ** 2 * (a - 1) * (a - 0.5)) * G1) +
                                            a0 * k * h0 * (-(2 * F1 + k * (4 * a - 3) * G1) * cos(phi) + (
                    2 * pi * F + k * F1 + 6. * G1) * sin(phi)))

    return CL, CT, CTM, CLM


if __name__ == "__main__":

    ### Integration model test

    tail_angle = radians(-3)    ## 5 degrees
    frequency = 5*tc            ## 5 Hz frequency

    model = ode_dydt(tail_angle, frequency)

    total_seconds = 1000        ## total planification seconds
    precision = 0.03            ## Model integration precision
    tSpan = [i*precision for i in range(int(total_seconds/precision))]

    Y0 = [1,0,0,0,0,0]          ## initial state (u,w,q,theta,x,h)

    Y = odeint(model, Y0, tSpan, args=(0,))   ## Integration result

    path = "images/"

    X = [y_i[-2] * c / 2 for y_i in Y]
    plt.plot(tSpan, X, "r-", linewidth=2, label="X(t)")
    plt.xlabel("time")
    plt.ylabel("X")
    plt.savefig(path+"X")
    plt.show()
    del X

    Z = [y_i[-1] * c / 2 for y_i in Y]
    plt.plot(tSpan, Z, "g-", linewidth=2, label="Z(t)")
    plt.xlabel("time")
    plt.ylabel("Z")
    plt.savefig(path+"Z")
    plt.show()
    del Z

    U = [y_i[0] for y_i in Y]
    plt.plot(tSpan, U, "r-", linewidth=2, label="X(t)")
    plt.xlabel("time")
    plt.ylabel("U")
    plt.savefig(path+"U")
    plt.show()
    del U

    W = [y_i[1] for y_i in Y]
    plt.plot(tSpan, W, "g-", linewidth=2, label="Z(t)")
    plt.xlabel("time")
    plt.ylabel("W")
    plt.savefig(path+"W")
    plt.show()
    del W

    Q = [y_i[2] for y_i in Y]
    plt.plot(tSpan, Q, "g-", linewidth=2, label="Z(t)")
    plt.xlabel("time")
    plt.ylabel("Q")
    plt.savefig(path+"Q")
    plt.show()
    del Q

    T = [y_i[3] for y_i in Y]
    plt.plot(tSpan, T, "g-", linewidth=2, label="Z(t)")
    plt.xlabel("time")
    plt.ylabel("Theta")
    plt.savefig(path+"Theta")
    plt.show()
    del T

