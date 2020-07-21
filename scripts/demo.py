'''
A sandbox file for adjusting parameters & exploring features of covasim
Author: Steve Kroiss
'''

# Imports
import covasim as cv    # cv.__file__ reports the location of the covasim package
import synthpops as sp
import sciris as sc
import numpy as np
import os



# Run with defaults parameters
sim = cv.Sim()
sim.pars
sim.run()
sim.plot()
sim.results.keys()
sim.results.r_eff
np.mean(sim.results.r_eff.values[sim.results.r_eff.values > 0])  # Reff ~0.204; ignore zero values if disease goes extinct



# Modify prognoses for one sim, then run multiple sims (1st with default settings) to compare results
prognoses = cv.get_prognoses()
prognoses['age_cutoffs'] = np.array([0,     1,    2,    3,    4,    5,   15,  30,  45,  90])
prognoses['sus_ORs']     = np.array([1.00,  1.00, 1.00, 1.00, 0.50, 0.0, 0.0, 0.0, 0.0, 0.0])
p1 = dict(
    beta = 0.2,
    pop_size = 10e3,
    )
p2 = dict(
    beta = 0.2,
    pop_size = 10e3,
    prognoses = prognoses,
    )
# Define the simulations
s1 = cv.Sim(pars=p1, label='Default age bins & prognoses')
s2 = cv.Sim(pars=p2, label='Modified age bins & prognoses')
# Combine and run the simulations
msim = cv.MultiSim([s1, s2])
msim.run()
msim.plot()



# # Double duration of inf2sym period
# pars = {}
# pars['dur'] = {}
# pars['dur']['exp2inf'] = {'dist': 'lognormal_int', 'par1': 4.6, 'par2': 4.8}  # Duration from exposed to infectious; see Lauer et al., https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7081172/, subtracting inf2sim duration
# pars['dur']['inf2sym'] = {'dist': 'lognormal_int', 'par1': 2.0, 'par2': 0.9}  # DOUBLE --- Duration from infectious to symptomatic; see Linton et al., https://doi.org/10.3390/jcm9020538
# pars['dur']['sym2sev'] = {'dist': 'lognormal_int', 'par1': 6.6, 'par2': 4.9}  # Duration from symptomatic to severe symptoms; see Linton et al., https://doi.org/10.3390/jcm9020538
# pars['dur']['sev2crit'] = {'dist': 'lognormal_int', 'par1': 3.0, 'par2': 7.4}  # Duration from severe symptoms to requiring ICU; see Wang et al., https://jamanetwork.com/journals/jama/fullarticle/2761044
# # Duration parameters: time for disease recovery
# pars['dur']['asym2rec'] = {'dist': 'lognormal_int', 'par1': 8.0, 'par2': 2.0}  # Duration for asymptomatic people to recover; see Wölfel et al., https://www.nature.com/articles/s41586-020-2196-x
# pars['dur']['mild2rec'] = {'dist': 'lognormal_int', 'par1': 8.0, 'par2': 2.0}  # Duration for people with mild symptoms to recover; see Wölfel et al., https://www.nature.com/articles/s41586-020-2196-x
# pars['dur']['sev2rec'] = {'dist': 'lognormal_int', 'par1': 14.0, 'par2': 2.4}  # Duration for people with severe symptoms to recover, 22.6 days total; see Verity et al., https://www.medrxiv.org/content/10.1101/2020.03.09.20033357v1.full.pdf
# pars['dur']['crit2rec'] = {'dist': 'lognormal_int', 'par1': 14.0, 'par2': 2.4}  # Duration for people with critical symptoms to recover, 22.6 days total; see Verity et al., https://www.medrxiv.org/content/10.1101/2020.03.09.20033357v1.full.pdf
# pars['dur']['crit2die'] = {'dist': 'lognormal_int', 'par1': 6.2, 'par2': 1.7}  # Duration from critical symptoms to death, 17.8 days total; see Verity et al., https://www.medrxiv.org/content/10.1101/2020.03.09.20033357v1.full.pdf
# sim = cv.Sim(pars=pars)
# sim.pars
# sim.run()
# sim.plot()
# sim.results.keys()
# np.mean(sim.results.r_eff.values[sim.results.r_eff.values > 0])      # Reff ~0.206; ignore zero values if disease goes extinct



# Set the location & use a hybrid pop approach
pars = {}
pars['location']        = 'pakistan'
pars['pop_type']        = 'hybrid'
pars['pop_size']        = 40e4      # Double pop size
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






