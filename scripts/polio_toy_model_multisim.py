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
    # pars_default = dict(
    #     beta = 0.2,
    #     pop_size = 10e3,
    #     )
    # pars_polio = dict(
    #     beta = 0.5,
    #     pop_size = 10e3,
    #     )
    # # Define the simulations
    # sim_default = cv.Sim(pars=pars_default, label='Default')
    # sim_polio = cv.Sim(pars=pars_polio, label='Polio')
    # # Combine and run the simulations
    # msim = cv.MultiSim([sim_default, sim_polio])
    # msim.run()
    # msim.plot()




    # Set the location & use a hybrid pop approach
    pars = {}
    pars['location']        = 'pakistan'
    pars['pop_type']        = 'hybrid'
    pars['pop_size']        = 50e3      # Small pop
    pars['beta']            = 0.232     # Double
    pars['dur']     = {'exp2inf': {'dist': 'uniform', 'par1': 1.0, 'par2': 1.0},            # Assume shedding starts day after infection
                       'inf2sym': {'dist': 'lognormal_int', 'par1': 9.0, 'par2': 4.0},      # Assume symptoms occur ~10 days after infection
                       'sym2sev': {'dist': 'lognormal_int', 'par1': 6.6, 'par2': 4.9},
                       'sev2crit': {'dist': 'lognormal_int', 'par1': 3.0, 'par2': 7.4},
                       'asym2rec': {'dist': 'lognormal_int', 'par1': 50.0, 'par2': 400.0},  # HUGE boost assuming infections are immunologically naive
                       'mild2rec': {'dist': 'lognormal_int', 'par1': 50.0, 'par2': 400.0},  # HUGE boost assuming infections are immunologically naive
                       'sev2rec': {'dist': 'lognormal_int', 'par1': 14.0, 'par2': 2.4},
                       'crit2rec': {'dist': 'lognormal_int', 'par1': 14.0, 'par2': 2.4},
                       'crit2die': {'dist': 'lognormal_int', 'par1': 6.2, 'par2': 1.7}}
    prognoses               = cv.get_prognoses()
    prognoses['sus_ORs']    = np.array([1.00,  1.00, 1.00, 1.00, 1.00, 1.00, 1.0, 1.0, 1.0, 1.0]) # Equally susceptible
    prognoses['symp_probs'] = np.array([0.005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])     # Equally susceptible
    prognoses['severe_probs'] = np.array([0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])    # Zero'd out
    prognoses['crit_probs'] = np.array([0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])      # Zero'd out
    prognoses['death_probs'] = np.array([0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])     # Zero'd out
    pars['prognoses']       = prognoses
    pars['start_day']       = '2020-01-01'
    pars['n_days']          = 365
    pars['quar_period']     = 140

    sim                     = cv.Sim(pars=pars, label='Default')  # Setup sim
    sim.initialize()
    people                  = cv.make_people(sim)  # Make the actual people for the simulation
    inds_o5                 = sc.findinds(people.age>=10)
    people.rel_sus[inds_o5] = 0
    inds_u5                 = sc.findinds(people.age<10)
    imm_frac                = 0.8
    inds_u5imm              = inds_u5[sc.findinds(np.random.uniform(low=0.0, high=1.0, size=len(inds_u5)) < imm_frac)]
    people.rel_sus[inds_u5imm] = 0
    sim.people              = sc.dcp(people)

    sim_TTQ                 = cv.Sim(pars=pars, label='TTQ')
    sim_TTQ.initialize()
    people_TTQ                  = cv.make_people(sim_TTQ)  # Make the actual people for the simulation
    inds_o5_TTQ                 = sc.findinds(people_TTQ.age>=10)
    people_TTQ.rel_sus[inds_o5_TTQ] = 0
    inds_u5_TTQ                 = sc.findinds(people_TTQ.age<10)
    imm_frac                = 0.8
    inds_u5imm_TTQ              = inds_u5_TTQ[sc.findinds(np.random.uniform(low=0.0, high=1.0, size=len(inds_u5_TTQ)) < imm_frac)]
    people_TTQ.rel_sus[inds_u5imm_TTQ] = 0
    sim_TTQ.people          = sc.dcp(people_TTQ)

    # sim['interventions'] = [cv.test_prob(symp_prob=1.00, asymp_prob=0.00)] # Test 100% of symptomatics and 0% of asymptomatics
    sim_TTQ['interventions'] = [cv.test_prob(symp_prob=1.00, asymp_prob=0.00)]  # , cv.contact_tracing()

    # Combine and run the simulations
    msim = cv.MultiSim([sim, sim_TTQ])
    msim.run()
    to_plot = ['cum_infections', 'cum_diagnoses', 'cum_symptomatic']
    # sim.plot(to_plot=to_plot)
    msim.plot(to_plot=to_plot)


# to_plot = ['cum_infections', 'cum_diagnoses', 'cum_symptomatic']
# sim.plot(to_plot=to_plot)
