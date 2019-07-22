# encoding: utf-8
"""
Transitions are always hard. 

Help us decide when we change from one state to another. 
"""

from numpy import pi, sqrt
import math as m
import numpy.random as random


def rate_to_prob(rate, obs_duration):
    """Convert a rate to a probability, based on an observation duration
    length and the assumption that the rate is for a Poisson process.
    We are asking, what is the probability that at least one Poisson
    distributed value would occur during the timestep.

    Parameters
    ----------
    rate : float
        a per ms rate to convert to probability
    obs_duration : float
        how many seconds we are watching over   

    Returns
    -------
    probability : float
        probability the event occurs during the ms we watched
    """
    return 1 - m.exp(-rate * obs_duration)


def reverse_rate(rate_12, free_energy_1, free_energy_2):
    """What is the balanced reverse rate from state 2 to state 1?
    
    Parameters
    ----------
    rate_12 : float
        transition rate from state 1 to state 2
    free_energy_1 : float
        free energy in state 1, state we are querying reversing into
    free_energy_2 : float
        free energy in state 2, state we are currently in

    Returns
    -------
    rate : float
        per ms rate of transitioning from state 1 back to state 2
    """
    try:
        energy_term = m.exp(free_energy_1 - free_energy_2)
    except OverflowError:
        energy_term = float("inf")
    try:
        rate = rate_12 / energy_term
    except ZeroDivisionError:
        rate = float("inf")
    return rate
