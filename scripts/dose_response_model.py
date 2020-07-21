'''
This file aims to adapt the dose response model from Famulare 2018 to covasim
It calculates susceptibility, shedding duration, and shedding amount based on immune titers
Author: Steve Kroiss
Source: Famulare et al. 2018 - Assessing the stability of polio eradication after the withdrawal of oral polio vaccine
'''

# Imports
import numpy as np
import scipy.special as special
import matplotlib.pyplot as plt



pars = {}

# OPV-equivalent antibody titer
pars['nab_nieve'] = 1.0  # Individual correlate of immunity. Nieve has N_Ab = 1
pars['nab_maximm'] = 2 ** 11  # Individual correlate of immunity. Nieve has N_Ab = 1

# Probability of shedding duration (equation S1)
pars['mu_s'] = 30.3         # Sabin median shedding duration (N_Ab = 1)
pars['sigma_s'] = 1.86      # Sabin scale parameter
pars['mu_wpv'] = 43.0       # WPV median shedding duration (NAb = 1)
pars['sigma_wpv'] = 1.69    # WPV scale parameter
pars['delta'] = 1.16        # median reduction per log2(NAb)

# peak shedding vs age (equation S2)
pars['sc_max'] = 6.7        # maximum stool concentration at age 7 months
pars['sc_min'] = 4.3        # maximum stool concentration at older ages
pars['tau'] = 12.0          # decay time constant of peak concentration with age

# peak shedding vs immunity (equation S3)
pars['kappa'] = 0.056       # shedding reduction with log2(NAb)

# shedding concentration vs time (equation S4)
pars['eta'] = 1.65          # n - location parameter
pars['nu'] = 0.17           # v - scale parameter
pars['epsilon'] = 0.32      # E - time-dependent scale

# dose response (equation S5)
pars['alpha'] = 0.44        # shape parameter
pars['gamma'] = 0.46        # immunity-dependent shape parameter exponent
pars['beta_s1'] = 14.0      # Sabin 1 scale parameter
pars['beta_s2'] = 8.0       # Sabin 2 scale parameter
pars['beta_s3'] = 18.0      # Sabin 3 scale parameter
pars['beta_wpv'] = 2.3      # WPV scale parameter

# waning immunity against infection (equation S6)
pars['lambda'] = 0.87       # immunity decay exponent




# Eq A - shedding duration given infection
# We assumed a log-normal survival distribution for the shedding duration given infection
def prob_shed_at_t(t, mu, delta, nab, sigma):
    return 0.5 * (1 - special.erf((np.log(t) - (np.log(mu) - np.log(delta) * np.log2(nab))) / (np.sqrt(2) * np.log(sigma))))
# Fig SA
t = np.arange(1,60,1)
plt.figure()
plt.subplot(111)
plt.plot(t, [prob_shed_at_t(x, mu = pars['mu_wpv'], delta = pars['delta'], nab = pars['nab_nieve'], sigma = pars['sigma_wpv']) for x in t], label='naive, WPV', color = 'blue')
plt.plot(t, [prob_shed_at_t(x, mu = pars['mu_s'], delta = pars['delta'], nab = pars['nab_nieve'], sigma = pars['sigma_s']) for x in t], label='naive, OPV', color = 'red')
plt.plot(t, [prob_shed_at_t(x, mu = pars['mu_s'], delta = pars['delta'], nab = pars['nab_maximm'], sigma = pars['sigma_s']) for x in t], label='max immunity, OPV', color = 'purple')
# p_naive_opv = plt.plot(prob_shed_at_t_naive_opv, label='naive, OPV', color = 'red')
# p_maximm_opv = plt.plot(prob_shed_at_t_maximm_opv, label='max immunity, OPV', color = 'purple')
plt.xlabel('shedding duration (days)')
plt.ylabel('fraction shedding')
plt.legend(loc="upper right")
plt.show()




# Eq B - age dependent peak shedding concentration
def peak_shed_conc(age):
    if (age >= 6):
        return 10**((pars['sc_max'] - pars['sc_min']) * np.exp((6-age)/pars['tau']) + pars['sc_min'])
    else:
        return 10**((pars['sc_max']))
# Fig 2B- not shown
age = range(1,65*12)  # age in months
plt.figure()
plt.subplot(111)
plt.plot(age, [peak_shed_conc(x) for x in age], label='naive, WPV')
# plt.yticks([1,10,100] + [max(y)])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('age (months)')
plt.ylabel('Peak shedding concentration (peak CID50/g)')
plt.show()




# Eq C - age dependent shedding concentration as a function of pre-challenge immunity
age = range(1,5*12)  # age in months
def shed_conc(age, nab):
    return 10**((1 - pars['kappa'] * np.log2(nab)) * np.log10(peak_shed_conc(age)))
# Fig - not shown
plt.figure()
plt.subplot(111)
plt.plot([shed_conc(age = x, nab = pars['nab_nieve']) for x in age], label='naive, WPV')
plt.plot([shed_conc(age = x, nab = pars['nab_maximm']) for x in age], label='max immunity, WPV')
plt.xlabel('age (months)')
plt.ylabel('Shedding concentration (peak CID50/g)')
plt.yscale('log')
plt.legend(loc="upper right")
plt.show()




# Eq D - Shedding concentration timeseries - age dependent for immunologically-naive individuals
t = range(1,30)  # timeseries
def shed_conc_over_time(t, age, nab):
    return max(10 ** 2.6, shed_conc(age, nab) * ( np.exp(pars['eta'] - (pars['nu'] ** 2)/2 - ((np.log(t)-pars['eta']) ** 2) / (2 * (pars['nu'] + pars['epsilon'] * np.log(t)) ** 2)) / t))
# Fig 2C
plt.figure()
plt.subplot(111)
plt.plot([shed_conc_over_time(x, age = 5, nab = pars['nab_nieve']) for x in t], label='5yo, naive, WPV')
plt.plot([shed_conc_over_time(x, age = 5, nab = 256) for x in t], label='5yo, max imm, WPV', color = 'blue')

# plt.plot([shed_conc(age = x, nab = pars['nab_maximm']) for x in age], label='max immunity, WPV')
plt.yscale('log')
plt.xlabel('time (days)')
plt.ylabel('Shedding concentration (CID50/g)')
plt.show()




# Eq E - oral susceptibility to infection from OPV challenge
def prob_inf_given_dose(dose, nab, beta):
    return 1 - (1 + (dose/beta)) ** (-pars['alpha'] * (nab ** (-pars['gamma'])))
# Fig 3C - Fraction shedding at vx doses (10**5-6 CID50) decreases with increasing OPV-equivalent antibody titer
nab_range = np.arange(1,4096, 1)
plt.figure()
plt.subplot(111)
plt.plot(nab_range,
         [prob_inf_given_dose(dose = 10**5.5, nab = x, beta = pars['beta_wpv']) for x in nab_range],
         label='WPV')
plt.xscale('log')
plt.xlabel('OPV-equivalent antibody titer')
plt.ylabel('Fraction shedding at vx dose (10**5-6 CID50)')
plt.legend(loc="lower right")
plt.show()
# Fig 3D - Dose-response model for immunologically naive (Nab=1), heterotypic bOPV and upper-bound IPV-only (NAb = 8),
# and typical tOPV or post-IPV-boosting
dose_range = 10 ** np.arange(-1,5.5,0.10)
plt.figure()
plt.subplot(111)
plt.plot(dose_range,
         [prob_inf_given_dose(dose = x, nab = 1, beta = pars['beta_s2']) for x in dose_range],
         label='naive (Nab=1)', color = 'red')
plt.plot(dose_range,
         [prob_inf_given_dose(dose = x, nab = 8, beta = pars['beta_s2']) for x in dose_range],
         label='bOPV+IPV (Nab=8)', color = 'green')
plt.plot(dose_range,
         [prob_inf_given_dose(dose = x, nab = 256, beta = pars['beta_s2']) for x in dose_range],
         label='tOPV (Nab=256)', color = 'blue')
plt.xscale('log')
plt.xlabel('Dose (CID50)')
plt.ylabel('Fraction shedding ')
plt.legend(loc="lower right")
plt.title('Dose-response model for SL2')
plt.show()
# WPV version of Fig 3D - Dose-response model for immunologically naive (Nab=1), heterotypic bOPV and upper-bound IPV-only (NAb = 8),
# and typical tOPV or post-IPV-boosting
dose_range = 10 ** np.arange(-1,5.5,0.10)
plt.figure()
plt.subplot(111)
plt.plot(dose_range,
         [prob_inf_given_dose(dose = x, nab = 1, beta = pars['beta_wpv']) for x in dose_range],
         label='Nab=1 (naive)', color = 'red')
plt.plot(dose_range,
         [prob_inf_given_dose(dose = x, nab = 8, beta = pars['beta_wpv']) for x in dose_range],
         label='Nab=8', color = 'green')
plt.plot(dose_range,
         [prob_inf_given_dose(dose = x, nab = 256, beta = pars['beta_wpv']) for x in dose_range],
         label='Nab=256', color = 'blue')
plt.xscale('log')
plt.xlabel('Dose (CID50)')
plt.ylabel('Fraction shedding ')
plt.legend(loc="lower right")
plt.title('Dose-response model for WPV')
plt.show()


# Eq F - Waning immunity against infection
def waned_imm(t, nab):
    return max(1, nab * t ** - pars['lambda'])
# Fig 4 - Waning immunity against infection
t = np.arange(1,600,1)
plt.figure()
plt.subplot(111)
plt.plot(t,
         [waned_imm(t = x, nab = 256) for x in t],
         label='Starting N_Ab = 256', color = 'red')
# plt.xscale('log')
plt.yscale('log')
plt.xlabel('Months between last immunizing event and OPV challenge')
plt.ylabel('OPV-equivalent antibody titer')
# plt.title('Starting N_Ab = 256')
plt.legend(loc="upper right")
plt.show()




# Eq 1 - Probability of infection on day t=1 due to either mOPV or WPV exposure on day t=0
def prop_inf_in_setting(px, dose, nab, beta):
    return px * prob_inf_given_dose(dose=dose, nab=nab, beta=beta)
prop_inf_in_setting(px = 1, dose = 10**5.5, nab = 1, beta = pars['beta_wpv'])  # UP & Bihar WPV
prop_inf_in_setting(px = 0.79, dose = 10**5.5, nab = 1, beta = pars['beta_s1'])  # Houston Sabin 1
prop_inf_in_setting(px = 1, dose = 10**5.5, nab = 1, beta = pars['beta_wpv'])  # Louisiana WPV




# Prevalence after exposure
def prob_index_shedding_at_t(px, dose, nab, beta, t, mu, delta, sigma):
    return prop_inf_in_setting(px, dose, nab, beta) * prob_shed_at_t(t, mu, delta, nab, sigma)

# Daily household incidence
















# Estimate household member infection probabilities
def beta_t(t, age_index, nab_index, nab_household, beta, tih, dih):
    dose = tih * shed_conc_over_time(t, age_index, nab_index)
    return 1 - (1 - prob_inf_given_dose(dose, nab_household, beta)) ** dih
def prob_house_trans_at_t(t, age_index, nab_index, nab_household, beta, tih, dih):
    if t == 1:
        return beta_t(t, age_index, nab_index, nab_household, beta, tih, dih)
    else:
        return beta_t(t, age_index, nab_index, nab_household, beta, tih, dih) * np.prod([1-beta_t(x, age_index, nab_index, nab_household, beta, tih, dih) for x in range(1,t)])
def prob_house_inf_at_t():
    prob_house_trans_at_t(t, age_index, nab_index, nab_household, beta, tih, dih) *\
    prob_index_shedding_at_t(px, dose, nab, beta, t, mu, delta, sigma)
# Figure
t = 60
plt.figure()
plt.subplot(111)
plt.plot(np.arange(1,t),
         [prob_house_trans_at_t(t=x, age_index=1, nab_index=1, nab_household=1, beta=2.3, tih=5, dih=1) for x in range(1,t)],
         label='Index & household naive', color = 'red')
plt.plot(np.arange(1,t),
         [prob_house_trans_at_t(t=x, age_index=1, nab_index=8, nab_household=8, beta=2.3, tih=5, dih=1) for x in range(1,t)],
         label='Nab_index = 8 & Nab_house = 8', color = 'green')
plt.plot(np.arange(1,t),
         [prob_house_trans_at_t(t=x, age_index=1, nab_index=256, nab_household=256, beta=2.3, tih=5, dih=1) for x in range(1,t)],
         label='Nab_index = 256 & Nab_house = 256', color = 'blue')
plt.xlabel('Months between last immunizing event and OPV challenge')
plt.ylabel('OPV-equivalent antibody titer')
plt.legend(loc="lower right")
plt.show()




# Estimate household member prevalence by convolving daily incidence (assuming no reinfection) with the shedding duration ditribution
# Figure
t = 60
tmp = np.cumsum([prob_house_trans_at_t(t=x, age_index=1, nab_index=256, nab_household=256, beta=pars['beta_s1'], tih=5, dih=1) for x in np.arange(1,t)]) * [prob_shed_at_t(x, mu = pars['mu_s'], delta = pars['delta'], nab = 256, sigma = pars['sigma_s']) for x in range(1,t)]
plt.figure()
plt.subplot(111)
plt.plot(np.arange(0,t),
         [0] +[prob_index_shedding_at_t(px=1, dose=10**5.5, nab=256, beta=pars['beta_s1'], t=x, mu=pars['mu_s'], delta=pars['delta'], sigma=pars['sigma_s']) for x in np.arange(1,t)],
         label='Index child', color = 'red')
plt.plot(np.arange(1,t),
         tmp,
         label='Family member', color = 'blue')
plt.xlabel('Days since exposure')
plt.ylabel('Fraction shedding')
plt.title('Sabin 1 - Nab for index & family = 256')
plt.legend(loc="upper right")
plt.show()


[0] +[prob_index_shedding_at_t(px=1, dose=10**5.5, nab=256, beta=pars['beta_s1'], t=x, mu=pars['mu_s'], delta=pars['delta'], sigma=pars['sigma_s']) for x in np.arange(1,t)]


def calculate_shed_duration(u=30.3, delta=1.16, sigma=1.86, strain_type = 'WPV', prechallenge_immunity = 1):
    """probability of shedding given Nab at time t (days post infection);
    assumes that individual was infected at t = 0; time is measured in days
    Equation S1 in Famulare 2018 PLOS Bio paper
    delta_t = time (days) since last infection -- survival curve follows lognormal distribution"""
    if strain_type == 'WPV':
        u, sigma = 43.0, 1.69
    mu = np.log(u) - np.log(delta) * np.log2(prechallenge_immunity)
    std = np.log(sigma)
    # scipy stats has weird implementation of parameters
    # the shape parameter (s) is the same as the stdev
    # the scale parameter is the same as the e^(mu)
    q = np.random.uniform(0, 1)
    shed_duration = stats.lognorm.isf(q, s=std,
                                                 scale=np.exp(mu))
    return shed_duration
calculate_shed_duration()

# Repeatedly call function
tmp = [calculate_shed_duration() for _ in range(100)]  # The _ is convention for a variable whose value you don't care about.


# Plot histogram of shedding duration
n, bins, patches = plt.hist(x=[calculate_shed_duration(prechallenge_immunity = 1) for _ in range(100)], bins='auto', color='#0504aa',
                            alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Shedding duration')
plt.ylabel('Frequency')
plt.title('Naive')


# Plot histogram of shedding duration
n, bins, patches = plt.hist(x=[calculate_shed_duration(prechallenge_immunity = 256) for _ in range(100)], bins='auto', color='#0504aa',
                            alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Shedding duration')
plt.ylabel('Frequency')
plt.title('Nab = 256')