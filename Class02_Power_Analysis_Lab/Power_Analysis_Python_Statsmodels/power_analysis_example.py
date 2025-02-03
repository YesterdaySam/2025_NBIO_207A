#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 18:40:25 2021

@author: Shantanu
"""

import numpy as np

from statsmodels.stats.power import TTestIndPower
from statsmodels.stats.power import FTestAnovaPower

import matplotlib.pyplot as plt

# variables for power analysis

effect_size = 0.7

alpha = 0.05

power = 0.9

p_analysis = TTestIndPower()

sample_size = p_analysis.solve_power(effect_size=effect_size, alpha=alpha, power=power)

print("Required Sample Size: " + str(sample_size))

# power vs. number of observations

fig = plt.figure()

fig = TTestIndPower().plot_power(dep_var='nobs',

                                nobs= np.arange(10, 200),

                                effect_size=np.array([0.2, 0.5, 0.8]),

                                alpha=0.01,

                                title='Power of t-Test at variable effect sizes\n' + r'$\alpha = 0.01$')

plt.show()

# Anova
fig = plt.figure()

fig = FTestAnovaPower().plot_power(dep_var='nobs',

                                nobs= np.arange(10, 200),

                                effect_size=np.array([0.2, 0.5, 0.8]),

                                alpha=0.01,

                                title='Power of ANOVA F-test at variable effect sizes\n' + r'$\alpha = 0.01$')

plt.show()