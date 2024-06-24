"""Definitions for the Aspherical Kilonova Constraints` class."""
import astropy.constants as c
import numpy as np

from mosfit.constants import C_CGS, M_SUN_CGS, KM_CGS, G_CGS
from mosfit.modules.constraints.constraint import Constraint

R_NS_CGS = 20e5 # NS radius cm == 15 km

# G_CGS = c.G.cgs.value


class AsphKNConstraints(Constraint):
    """BNS constraints.

    1. Mej = ((5 * R_NS * c^2) / (3* G)) * (1 / sqrt(1 - (Vej / c)^2) / M_sun

    Free parameters are blue ejecta mass, blue ejecta 
    velocity,red ejecta mass and red ejecta velocity. 
    
    Constraint penalises ejecta mass and velocity combinations 
    (within each dynamical ejecta component) for which the 
    ejecta kinetic energy is less than the binding 
    energy of a NS, assuming a relatively stiff EoS.
    """

    def __init__(self, **kwargs):
        """Initialize module."""
        super(AsphKNConstraints, self).__init__(**kwargs)

    def process(self, **kwargs):
        """Process module. Add constraints below."""
        self._score_modifier = 0.0
        self._mblue = kwargs[self.key('mejecta_blue')]
        self._vblue = kwargs[self.key('vejecta_blue')]
        self._mred = kwargs[self.key('mejecta_red')]
        self._vred = kwargs[self.key('vejecta_red')]

        # for Kinetic energy to be greater than binding energy, this is the max Mej

        # 1
        blue_bound_mass = ((5.0*R_NS_CGS*C_CGS**2)/(3.0*G_CGS))*(1.0/np.sqrt(1.0 - (self._vblue * 1e5 / C_CGS)**2) - 1.0)/M_SUN_CGS 
        if self._mblue > blue_bound_mass:
            self._score_modifier -= (100. * (self._mblue-blue_bound_mass))**2

        # 2
        red_bound_mass = ((5.0*R_NS_CGS*C_CGS**2)/(3.0*G_CGS))*(1.0/np.sqrt(1.0 - (self._vred * 1e5 / C_CGS)**2) - 1.0)/M_SUN_CGS
        if self._mred > red_bound_mass:
            self._score_modifier -= (100. * (self._mred-red_bound_mass))**2

        return {self.key('score_modifier'): self._score_modifier}
