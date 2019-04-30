# encoding: utf-8
"""
To protect constants and serve unit conversions

Gather all the parameters that define the units in our universe, and the
particular section of it we are interested in, into one spot. 
"""

## Conversion
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
    """Universal constants we don't expect to change across runs"""

    temperature = 288  # Standard temperature, K
    eta = poise(0.0114)  # Viscosity of water at 288K
    boltzmann = joules(1.38 * 10 ** -23)  # From J/K to pN*nm/K
    kT = boltzmann * temperature
    avogadro = 6.022 * 10 ** 23


class world:
    """Things that are constant for our (simulated) world
    
    * timestep is the number of ms we advance each tick
    * D_cyto_corr comes from the difference in diffusion between water and
    eukaryotic cytoplasm (Swaminathan_1997_).

    .. _Swaminathan_1997: https://doi.org/10.1016/S0006-3495(97)78835-0
    """

    timestep = milliseconds(1)  # 1ms timestep
    D_cyto_corr = 1 / 3.2  # Cytoplasmic crowding sub-diffusion
