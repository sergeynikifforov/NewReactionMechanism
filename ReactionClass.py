# This class contains the parameters of reaction rate constant
import numpy as np
import math
k_b = 1.38e-23
r = 0.038
R = 1.987


class ReactionClass:

    def __init__(self, temperature_beta, temperature_alpha, energy_alpha):
        self.temperature_beta = temperature_beta
        self.t_alpha = temperature_alpha
        self.energy_alpha = energy_alpha

    def constant_rate_evaluate(self, temperature):
        return self.temperature_beta*temperature**self.t_alpha*math.exp(-self.energy_alpha/(R * temperature))


class ReactionClassWall:
    def __init__(self, epsilon, molar_mass):
        self.epsilon = epsilon
        self.molar_mass = molar_mass

    def constant_rate_evaluate_wall_k(self, temperature):
        return np.sqrt((8.*k_b*temperature)/(np.pi * self.molar_mass))

    def constant_rate_evaluate_wall(self, temperature):
        return 1./4.*self.constant_rate_evaluate_wall_k(temperature)*self.epsilon*3./r


class LimitsClass:
    def __init__(self, k1, k17b, k9_inf, k9_0, k_hydrogen, k_hydroperoxyl_radical):
        self.k1 = k1
        self.k17b = k17b
        self.k9_inf = k9_inf
        self.k9_0 = k9_0
        self.k_hydrogen = k_hydrogen
        self.k_hydroperoxyl_radical = k_hydroperoxyl_radical
        self.ksi_h2 = 2./3.
        self.ksi_o2 = 1./3.

    # coeff of cube equation
    def a_t(self):
        return 2.*self.ksi_h2*self.ksi_o2*(1.+self.k1/self.k9_inf)

    def c_t(self):
        return 2.*self.ksi_o2*(self.k1*self.k_hydroperoxyl_radical)/(self.k9_0*self.k17b) - self.ksi_h2*self.k_hydrogen/self.k9_0 - (self.k_hydroperoxyl_radical*self.k_hydrogen)/(self.k9_inf*self.k17b)

    def b_t(self):
        return self.ksi_o2*self.k_hydroperoxyl_radical/self.k17b*(2.*self.k1/self.k9_inf - 1.) - self.ksi_h2*self.k_hydrogen/self.k9_inf + 2.*self.ksi_o2*self.ksi_h2*self.k1/self.k9_0

    def d_t(self):
        return -(self.k_hydrogen*self.k_hydroperoxyl_radical)/(self.k9_0*self.k17b)

    def cubic_equation_solve(self, temperature):
        a_t = self.a_t()/(4.18*R*temperature*1e6)**3
        b_t = self.b_t()/(4.18*R*temperature*1e6)**2
        c_t = self.c_t()/(4.18*R*temperature*1e6)
        d_t = self.d_t()
        # coeff of cube equation
        coeff = [a_t, b_t, c_t, d_t]
        # solve cube equation
        res = np.roots(coeff)
        mask = np.isreal(res)
        new_res = list(map(lambda x: x.real, res[mask]))
        return new_res
