import math
from project.constants import boltzmann_eV

class Defect_at_Site:

    """
    The Defect_at_Site class contains the information about each defect at each site, its valence, mole fraction and the segregation energy for that defect at that site. 
    This functions in this class combine to give the correct statistics for site occupation in solid electrolytes. It is a 'Fermi-Dirac' like approach as described in 'Origin of Space Charge in Grain boundaries of Protoin-Conductibng BaZrO3' by Helgee, Lindmann and Wahnstrom in 2012 (published in Fuel Cells). 
    The overall equation takes traditional Boltzmann statistics (def Boltzmann one), and adds normalisation so that each site cannot have a probability of being occupied that is any higher than 1 (at 1 the site is fully occupied) and accounts for any competition of like charged defects.
    There is also an option to fix any defect at their sites which is appropriate for defects which are not mobile.   
    """
    def __init__( self, label, valence, mole_fraction, energy, site, fixed = False ):
        self.label = label
        self.valence = valence
        self.mole_fraction = mole_fraction
        self.energy = energy
        self.site = site
        self.fixed = fixed 
   

    def potential_energy( self, phi ):
        """
        Potential energy for the defect at this site.

        Args:
            phi (float): electrostatic potential at this site

        Returns:
            (float): The electrochemical potential
        """
        return ( phi * self.valence ) + self.energy


    def boltzmann_one( self, phi, temp ):
        """
        Boltzmann statistics calculation

        (In LaTeX notation) \exp^{ \frac{ \Phi z + \Delta E }{ k_BT } }

        Args:
            phi (float): Electrostatic potential.
            temp (float): Temperature of calculation.
    
        Returns:
            (float): Boltzmann statistics
   
        """
        return math.exp( - self.potential_energy( phi ) / ( boltzmann_eV * temp ) )

    def boltzmann_two( self, phi, temp ):
        """
        Boltzmann statistics calculation
 
        (In LaTeX notation) x ( \exp^{\frac{\Phi z + \Delta E}{K_BT} } )

        Args:
            phi (float): Electrostatic potential.
            temp (float): Temperature of calculation.
    
        Returns:
            (float): Boltzmann statistics * mole fraction
   
        """
        return self.mole_fraction * self.boltzmann_one( phi, temp )

    def boltzmann_three( self, phi, temp ):
        """
        Boltzmann statistics calculation

        ( In LaTeX notation ) x (\exp^{\frac{ \Phi z + \Delta E}{ K_BT} } - 1 )

        Args:
            phi (float): Electrostatic potential.
            temp (float): Temperature of calculation.
    
        Returns:
            (float): ( Boltzmann statistics - 1 ) * mole fraction.
   
        """
        return self.mole_fraction * ( self.boltzmann_one( phi, temp ) - 1.0 )
