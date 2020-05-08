import eos, activity
from chemsep_operation import get_chemical

def fugacity_vapor(components, temp, pressure, fractions, method=None, kij_input = None, kij_tune=None):
    """ Equation of state solver for vapor phase.
    :param components: Array that contains chemicals.
    :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
    :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
    """
    methods = {"SRK":eos.SRK.phi_vapor,"Ideal":eos.Ideal.phi_vapor,"PR76":eos.PR76.phi_vapor,"PR78":eos.PR76.phi_vapor,"RK":eos.RK.phi_vapor}
    if method != None:
        phi = methods[method](components, temp, pressure, fractions, kij_input, kij_tune)
    else:
        phi = methods["PR78"](components, temp, pressure, fractions, kij_input, kij_tune)

    return phi

def fugacity_liquid(components, temperature, fractions, method=None, kij_input = None, kij_tune=None):
    
    methods = {"SRK":eos.SRK.phi_liquid,"Ideal":activity.Ideal.gamma, "NRTL":activity.NRTL.gamma,"Uniquac":activity.Uniquac.gamma,"Unifac":activity.Unifac.gamma,"Dortmund":activity.Dortmund.gamma}
    if method != None and method != "SRK":
        gamma = methods[method](components, temperature, fractions)
    else:
        gamma = methods["SRK"](components, temperature, 101325, fractions, kij_input, kij_tune)

    return gamma


"""
***EXAMPLE***
Library index for Water:1921
Library index for Ethanol:1102

chem1 = get_chemical("1921")  --Get chemical objects
chem2 = get_chemical("1102")

fugacity_liquid([chem1,chem2], 320, [0.7,0.3], method="Unifac", kij_input = None, kij_tune=None)

"""