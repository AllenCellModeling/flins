# encoding: utf-8
"""
_units: To protect constants and serve unit conversions
CDW 2019
"""


def poise(input_in_poise):
    """Convert the input viscosity, in poise, to g/nm*s"""
    return input_in_poise * 10 ** -7


def joules(input_in_joules):
    """Convert energy, in joules, to pN*nm"""
    return input_in_joules * 10 ** 21


def kcal(input_in_kcal):
    """Convert energy, in kcal, to pN*nm"""
    return input_in_kcal * 4.184 * 10 ** 24


def milliseconds(input_in_ms):
    """Convert time, in ms, to seconds"""
    return input_in_ms * 0.001


class constants:
    temperature = 288  # Standard temperature, K
    eta = poise(0.0114)  # Viscosity of water at 288K
    timestep = milliseconds(1)  # 1ms timestep
    boltzmann = joules(1.38 * 10**-23)  # From J/K to pN*nm/K
    kT = boltzmann * temperature
    avogadro = 6.022 * 10**23
