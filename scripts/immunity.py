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
import scipy.stats as stats
import math

class ImmunoInfection:
    '''typical infection/immune dynamics for poliovirus, this is also where pathogen-side genetics will be put (either that or split it into its own class)
    This is keyed towards polio -- currently completely avirulent!
    dt is measured in days!'''
    sabin_scale_parameters ={'S1': 14, 'S2': 8, 'S3': 18, 'WPV': 2.3}
    strain_map = {'S1': 0, 'S2': 1, 'S3': 2, 'WPV': 3}

    def __init__(self, strain_type):
        self.naive = True
        self.IPV_naive = False  # received an IPV dose, but does not shed
        self.is_shed = False
        self.strain_type = strain_type
        self.t_infection = 0  # days, time since last infection
        self.current_immunity = 1  # Naive
        self.prechallenge_immunity = 1
        self.postchallenge_peak_immunity = None
        self.shed_duration = 0

    def p_infection(self, dose=10 ** 6, alpha=0.44, gamma=0.46):  # dose values are for a typical vaccine dose
        '''Dose response model. Pr(infection | dose, immunity).
        Immunity is an abstraction of Nab * theta
        Immunity based on probability of infection'''
        beta = ImmunoInfection.sabin_scale_parameters[self.strain_type]
        return (1 - (1 + dose / beta) ** (-alpha * (self.current_immunity) ** -gamma))  # * 0.554

    def attempt_infection(self, dose=10 ** 6, modifier=1):
        if (np.random.random() < modifier * self.p_infection(
                dose)):  # this actually causes the probability of infection to go beyond 1.0
            # however, np.random.random() only draws from 0-1 so it is the same as maxing it out
            self.is_shed = True
            self.naive, self.IPV_naive = False, False
            self.prechallenge_immunity = self.current_immunity
            self.t_infection = 0  # infection is reset
            self.calculate_peak_immunity()
            self.calculate_current_immunity()
            self.calculate_shed_duration()
            return 1
        else:
            return 0

    def calculate_peak_immunity(self):
        '''immunity immediately post infection'''
        self.postchallenge_peak_immunity = self.prechallenge_immunity * max(1,
                                                                            self.theta_Nab())  # prevent immunity from decreasing due to challenge

    def calculate_current_immunity(self, rate=0.87):
        '''immunity after t months have passed since exposure'''
        if math.floor(self.t_infection / 30) != 0:
            self.current_immunity = max(1, self.postchallenge_peak_immunity * (
                        self.t_infection / 30) ** -rate)  # t in this equation is measured in months, not days
        else:
            self.current_immunity = max(1, self.postchallenge_peak_immunity)

    def theta_Nab(self, a=4.82, b=-0.30, c=3.31, d=-0.32):
        # prepopulated with the mle estimates from the mle model using the Gelfand data set- only biological theta is modeled
        Nab = self.prechallenge_immunity
        mean = a + b * np.log2(Nab)
        stdev = np.sqrt(max(c + d * np.log2(Nab), 0))
        return np.exp(np.random.normal(loc=mean, scale=stdev))



    def calculate_shed_duration(self, u=30.3, delta=1.16, sigma=1.86):
        """probability of shedding given Nab at time t (days post infection);
        assumes that individual was infected at t = 0; time is measured in days
        Equation S1 in Famulare 2018 PLOS Bio paper
        delta_t = time (days) since last infection -- survival curve follows lognormal distribution"""
        if self.strain_type == 'WPV':
            u, sigma = 43.0, 1.69
        mu = np.log(u) - np.log(delta) * np.log2(self.prechallenge_immunity)
        std = np.log(sigma)
        # scipy stats has weird implementation of parameters
        # the shape parameter (s) is the same as the stdev
        # the scale parameter is the same as the e^(mu)
        q = np.random.uniform(0, 1)
        self.shed_duration = stats.lognorm.isf(q, s=std,
                                                     scale=np.exp(mu))  # inverse lognormal survival curve sampling

    def viral_shed(self, age, eta=1.65, v=0.17, epsilon=0.32):
        '''virus shed per gram, time is in months!'''
        if self.is_shed == True:
            predicted_concentration = 10 ** self.peak_cid50(age * 12) * np.exp(
                eta - 0.5 * v ** 2 - ((np.log10(self.t_infection) - eta) ** 2) / (
                            2 * (v + epsilon * np.log10(self.t_infection)) ** 2)) / self.t_infection
            return max(10 ** 2.6, predicted_concentration)
        else:
            return 0

    def peak_cid50(self, age, k=0.056, Smax=6.7, Smin=4.3, tau=12):
        '''returns the peak log10(cid50/g) given prior immunity, age is in months!'''
        if age >= 6:
            peak_cid50_naive = (Smax - Smin) * np.exp((7 - age) / tau) + Smin
        else:
            peak_cid50_naive = Smax

        return (1 - k * np.log2(self.prechallenge_immunity)) * peak_cid50_naive

    def update(self, dt_years):
        # dt_years because the household model operates in terms of years
        if self.naive != True:
            self.t_infection += dt_years * 365  # days!
            self.calculate_current_immunity()
            if self.t_infection >= self.shed_duration:
                self.is_shed = False


