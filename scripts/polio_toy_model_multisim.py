'''
A scratch file for adjusting parameters from covasim to represent polio
Author: Steve Kroiss
'''

# Imports
import covasim as cv    # cv.__file__ reports the location of the covasim package
# import synthpops as sp
import sciris as sc
import numpy as np
import os
import matplotlib.pyplot as plt




if __name__ == '__main__':
    pars_default = dict(
        beta = 0.2,
        pop_size = 10e3,
        )
    pars_polio = dict(
        beta = 0.5,
        pop_size = 10e3,
        )
    # Define the simulations
    sim_default = cv.Sim(pars=pars_default, label='Default')
    sim_polio = cv.Sim(pars=pars_polio, label='Polio')
    # Combine and run the simulations
    msim = cv.MultiSim([sim_default, sim_polio])
    msim.run()
    msim.plot()
