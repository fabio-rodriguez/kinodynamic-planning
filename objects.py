from model import Uc, c
from math import sqrt


class FlightState:
    '''
        Class for store flight states information:
            tail_angle: tail angle for the ornithopter (radians)
            fa: flap frequency of wings (adimentional)
            u: horizontal velocity
            v: vertital velocity
            theta: pitch angle
            omega: angular velocity
            x: position in X axis (XZ flat)
            z: position in Z axis (XZ flat)
    '''

    def __init__(self, tail_angle, fa, u, v, theta, omega, x, z):
        self.x = x
        self.z = z
        self.u = u
        self.v = v
        self.theta = theta
        self.omega = omega
        self.tail_angle = tail_angle
        self.fa = fa
        self.cost = 0

    def dimensionalXvalue(self):
        return self.x * c/2

    def dimensionalZvalue(self):
        return self.z * c/2

    def Ub_value(self):
        return sqrt(self.v**2 + self.u**2)

    def real_velocity_value(self):
        return self.Ub_value() * Uc

    def increment_cost(self, cost0, P):

        self.cost = cost0 + P
        return self.cost

    def __str__(self):
        return "X: "+str(self.x) + " Z: "+str(self.z) + "\nV: "+str(self.v) + " U: "+str(self.u) + \
               "\nTHETA: " + str(self.theta) + " OMEGA: " + str(self.omega) +"\nTAIL_ANGLE: " + str(self.tail_angle) + \
               " COST: " + str(self.cost)

