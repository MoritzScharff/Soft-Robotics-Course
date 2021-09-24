"""This module contains all functions regarding Steigenberger [2012] Eval."""

# import section

import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve

def v_lin(L, d, E, F):
    """Use this function to calculate the elastica of a whisker under tip load considering linear deformation behavior.

       :param L: Length of the whisker [m]
       :type L: float
       :param d: diameter of the whisker [m]
       :type d: float
       :param E: Young's modulus [N/m**2]
       :type E: float
       :param F: Force/load [N]

       :retrun: elastica [m]
       :rtype: numpy arrays for x- and y-coordinates

    """

    Iz = np.pi/64 * d**4 # Second moment of area [m**4]

    x_vals = np.linspace(0,L,100)

    y_vals = x_vals*(-F*L*x_vals/2 + F*x_vals**2/6)/(E*Iz)

    elastica = np.empty((100,2))
    elastica[:,0] = x_vals
    elastica[:,1] = y_vals

    return elastica

def v_nonlin(L, d, E, F):
    """Use this function to calculate the elastica of a whisker under tip load considering nonlinear deformation behavior.

       :param L: Length of the whisker [m]
       :type L: float
       :param d: diameter of the whisker [m]
       :type d: float
       :param E: Young's modulus [N/m**2]
       :type E: float
       :param F: Force/load [N]

       :retrun: elastica [m]
       :rtype: numpy arrays for x- and y-coordinates

    """

    # Integration Errors
    abserr = 1.0e-12
    relerr = 1.0e-10

    # Optimization Error
    tol=1.0e-8

    Iz = np.pi/64 * d**4 # Second moment of area [m**4]

    p = L, d, E, F, Iz # build parameter vector
    s = np.linspace(0,L,100) # arclength vector
    BC = 0, 0, 0 # Boundary conditions x(0) = 0, y(0) = 0, phi(0) = 0

    sol = fsolve(ode_sol, 0.7*L, args=(p, BC, s, abserr, relerr), xtol=tol) # Returns the roots of the (non-linear) equations

    elastica_nonlin = odeint(ode_sys, BC, s, args=(p, sol[-1:]), atol=abserr, rtol=relerr)

    return elastica_nonlin


def ode_sol(opti_para,p,BC,s,abserr,relerr):
    """Apply the shooting method by using a root finding optimization algorithm.

        :param opti_param: optimization parameter
        :type opti_param: progam defined
        :param p: set of constant belonging to the specific problem
        :type p: tuple consiting of anything
        :param BC: boundary conditions belonging to the specific problem
        :type BC: tuple
        :param s: arc lentgh
        :type s: vector
        :param abserr: absolute error for integration
        :type abserr: float
        :param relerr: relative error for integration
        :type relerr: float

        :return: optimization value / evaluated optimization criterion
    """

    # unknown parameter
    xL = opti_para[0]

    # numeric integration
    elastica_nonlin = odeint(ode_sys, BC, s, args=(p, xL), atol=abserr, rtol=relerr)

    # optimization criterion
    eq = elastica_nonlin[-1:,0]-xL

    return eq[0]

def ode_sys(BC, s, p, xL):
    """Get the elastica of the whisker.

        :param BC: boundary conditions belonging to the specific problem
        :type BC: tuple
        :param s: arc lentgh
        :type s: vector
        :param p: set of constant belonging to the specific problem
        :type p: tuple consiting of anything
        :param xL: unkown position of the tip of the whisker
        :type xL: flot

        :return: Defined ODE system
        :rytpe: ODE set

    """

    L, d, E, F, Iz = p
    x, y, phi = BC

    v = [np.cos(phi),
         np.sin(phi),
         F/(E*Iz)*(x-xL)]

    return v
