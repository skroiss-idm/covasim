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



# Set the location & use a hybrid pop approach
pars = {}
pars['location']        = 'pakistan'
pars['pop_type']        = 'hybrid'
pars['pop_size']        = 10e3      # Small pop
pars['beta']            = 0.032     # Double
pars['dur']     = {'exp2inf': {'dist': 'uniform', 'par1': 1.0, 'par2': 1.0},
                   'inf2sym': {'dist': 'lognormal_int', 'par1': 9.0, 'par2': 4.0},
                   'sym2sev': {'dist': 'lognormal_int', 'par1': 6.6, 'par2': 4.9},
                   'sev2crit': {'dist': 'lognormal_int', 'par1': 3.0, 'par2': 7.4},
                   'asym2rec': {'dist': 'lognormal_int', 'par1': 8.0, 'par2': 2.0},
                   'mild2rec': {'dist': 'lognormal_int', 'par1': 8.0, 'par2': 2.0},
                   'sev2rec': {'dist': 'lognormal_int', 'par1': 14.0, 'par2': 2.4},
                   'crit2rec': {'dist': 'lognormal_int', 'par1': 14.0, 'par2': 2.4},
                   'crit2die': {'dist': 'lognormal_int', 'par1': 6.2, 'par2': 1.7}}
prognoses = cv.get_prognoses()
prognoses['sus_ORs']    = np.array([1.00,  1.00, 1.00, 1.00, 1.00, 1.00, 1.0, 1.0, 1.0, 1.0])
prognoses['symp_probs'] = np.array([0.005 , 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0])
pars['prognoses']       = prognoses
pars['start_day']       = '2020-01-01'
pars['n_days']          = 60
pars['quar_period']     = 14
sim                     = cv.Sim(pars=pars)  # Setup sim
people                  = cv.make_people(sim)  # Make the actual people for the simulation
over5_indx              = people['age'] >= 5  # Find over 5's
under5_indx             = people['age'] < 5  # Find under 5's
imm_frac                = 0.8
under5_imm_indx         = list(under5_indx) and list(np.random.uniform(low=0.0, high=1.0, size=len(under5_indx)) < imm_frac) # Randomly determine immunity for each individual
people['recovered'][over5_indx] = True  # Make anyone over 5 recovered
people['recovered'][under5_imm_indx] = True  # Make anyone over 5 recovered
unique, counts = np.unique(people['recovered'], return_counts=True)  # tabulate number of recovereds
print(np.asarray((unique, counts)).T)
# sim = cv.Sim(pars = pars, popfile=people, load_pop=True)
# sim.people = people
# unique, counts = np.unique(sim.people['recovered'], return_counts=True)
# print(np.asarray((unique, counts)).T)
sim.run()
sim.plot()
sim.results.keys()
sim.results.r_eff
print(np.mean(sim.results.r_eff.values[sim.results.r_eff.values > 0]))  # Reff ~0.39; ignore zero values if disease goes extinct























###########################################
#
#
# WARNING - What follows below is scary, poorly written scratch code. Don't judge me.
#
#
###########################################

















par1 = 9.0
par2 = 4.0
size = 1000
mean  = np.log(par1**2 / np.sqrt(par2 + par1**2)) # Computes the mean of the underlying normal distribution
sigma = np.sqrt(np.log(par2/par1**2 + 1)) # Computes sigma for the underlying normal distribution
samples = np.random.lognormal(mean=mean, sigma=sigma, size=size)
samples = np.round(samples)
plt.hist(samples, bins=30)  # `density=False` would make counts


np.random.uniform(low=1.0, high=1.0, size=size)



# Set the location & use a hybrid pop approach
pars_default            = {}
pars_polio              = {}
pars_polio['location']        = 'pakistan'
pars_polio['pop_type']        = 'hybrid'
pars_polio['pop_size']        = 10e4
# pars_polio['dur'] = {}
# pars_polio['dur']['exp2inf'] = {'dist': 'lognormal_int', 'par1': 1.6, 'par2': 4.8}  # Duration from exposed to infectious; see Lauer et al., https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7081172/, subtracting inf2sim duration
# pars['dur']['inf2sym'] = {'dist': 'lognormal_int', 'par1': 2.0, 'par2': 0.9}  # DOUBLE --- Duration from infectious to symptomatic; see Linton et al., https://doi.org/10.3390/jcm9020538



sim = cv.Sim(pars=pars_polio)
sim.pars
sim.run()
sim.plot()
sim.results.keys()
sim.results.r_eff
np.mean(sim.results.r_eff.values[sim.results.r_eff.values > 0])  # Reff ~0.39; ignore zero values if disease goes extinct



# prognoses['sus_ORs']     = np.array([1.00,  1.00, 1.00, 1.00, 0.50, 0.0, 0.0, 0.0, 0.0, 0.0])
# p1 = dict(
#     beta = 0.2,
#     pop_size = 10e3,
#     )
# p2 = dict(
#     beta = 0.2,
#     pop_size = 10e3,
#     prognoses = prognoses,
#     )
# Define the simulations
sim_default = cv.Sim(pars=pars_default, label='Default')
sim_polio = cv.Sim(pars=pars_polio, label='Polio')
# Combine and run the simulations
msim = cv.MultiSim([sim_default, sim_polio])
msim.run()
msim.plot()


sim = cv.Sim(pars=pars)
sim.pars
sim.run()
sim.plot()
sim.results.keys()
sim.results.r_eff
np.mean(sim.results.r_eff.values[sim.results.r_eff.values > 0])  # Reff ~0.39; ignore zero values if disease goes extinct



# Example of how to extract household ids from synthpops - Parts of code copied from generate_contact_network_with_microstructure.py in synthpops/examples
sp.validate()
datadir = sp.datadir  # point datadir where your data folder lives
location = 'seattle_metro'
state_location = 'Washington'
country_location = 'usa'
sheet_name = 'United States of America'
n = 10000
verbose = False
plot = False
# this will generate a population with microstructure and age demographics that approximate those of the location selected
# also saves to file in:
#    datadir/demographics/contact_matrices_152_countries/state_location/
popdict = sp.generate_synthetic_population(n, datadir, location=location, state_location=state_location,
                                 country_location=country_location, sheet_name=sheet_name, verbose=verbose, plot=plot,
                                 return_popdict = True)
# load that population into a dictionary of individuals who know who their contacts are
options_args = {'use_microstructure': True}
network_distr_args = {'Npop': n}
# Extract individuals and their contacts
contacts = sp.make_contacts(location=location, state_location=state_location, country_location=country_location,
                            options_args=options_args, network_distr_args=network_distr_args)
# show_layers(contacts, show_ages=True)
uids = popdict.keys()                       # Extract keys
uids = [uid for uid in uids]                # Convert to list
popdict[uids[0]].keys()                     # See keys for an individual
popdict[uids[0]]['hhid']                    # See that individual's household id (hhid)
popdict[uids[0]]['contacts'].keys()         # See keys for that individual's contacts
hhids = [popdict[x]['hhid'] for x in uids]  # Extract household ids
ages = [popdict[x]['age'] for x in uids]    # Extract ages
idx_u5 = [i for i,v in enumerate(np.array(ages) < 5) if v]  # Find hhids with individuals under 5
hhids_u5 = np.unique(np.array(hhids)[idx_u5])






