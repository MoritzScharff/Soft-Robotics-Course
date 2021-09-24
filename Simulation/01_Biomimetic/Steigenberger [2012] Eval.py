"""Soft Robotics 2021-2.
   Application example of Biomimetics in soft robotics: Steigenberger (2012)

   Author: Moritz Scharff
   Contact: moritz.scharff@pucp.edu.pe

   :param L: Length of the whisker [m]
   :type L: float
   :param d: diameter of the whisker [m]
   :type d: float
   :param E: Young's modulus [N/m**2]
   :type E: float
   :param F: Force/load [N]

   :retrun: elasticae [m]
   :rtype: numpy arrays for x- and y-coordinates

"""

# import section
import numpy as np
import scipy as sci
import custom_lib as get
import matplotlib.pyplot as plt

# defining variables
L = 0.1
d = 1E-3
E = 2.1E11
F = 1.5

# get solution of linear / small deflection
elastica_lin = get.v_lin(L, d, E, F)

# get solution of nonlinear / large deflection by optimization/shooting method
elastica_nonlin = get.v_nonlin(L, d, E, F)


# plotting the results to compare
fig, axs = plt.subplots(1, 2)
axs[0].plot(elastica_lin[:,0], elastica_lin[:,1])
axs[1].plot(elastica_lin[:,0], elastica_lin[:,1])
axs[1].plot(elastica_nonlin[:,0], elastica_nonlin[:,1])
axs[1].axis('equal')
axs[0].grid()
axs[1].grid()
plt.show()
