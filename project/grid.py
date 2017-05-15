import numpy as np
import math
from project.constants import boltzmann_eV
from bisect import bisect_left
from scipy.interpolate import griddata
from project.constants import boltzmann_const

def phi_at_x( phi, grid, x ):
    """
    Assigns each site x coordinate a grid point and returns the electrostatic potential at the grid point clostest to the x coordinate.

    Args:
        phi (np.array): electrostatic potential on 1D grid. 
        grid (np.array): 1D grid of ordered numbers over a region to compare with x coordinates.
        x (float): Site x coordinate.

    Returns:
        phi[index] (float): The electrostatic potential at the x coordinate with position [index].
    """
    index = index_of_grid_at_x( grid, x )
    return phi[ index ]

def energy_at_x( energy, grid, x ):
    """
    Assigns each site x coordinate a grid point and returns the segregation energy at the grid point clostest to the x coordinate.

    Args:
        energy (np.array): Segregation energies on 1D grid. 
        grid (np.array): 1D grid of ordered numbers over a region.
        x (float): Site x coordinate.

    Returns:
        energy[index] (float): The segregation energy at the x coordinate with position [index].
    """

    index = index_of_grid_at_x( grid, x )
    return energy[ index ]

def index_of_grid_at_x( grid, x ):
    """ 
    Assigns each site x coordinate to a position on a regularly or irregularly spaced grid. 
    Returns the index of the grid point clostest to the value x 

    Args:
        grid (np.array): 1D grid of ordered numbers over a region.
        x (float): Site x coordinate

    Returns:
        closest_index (int): Index of grid position closest to the site x coordinate.
    """
    return closest_index( grid, x )

def closest_index(myList, myNumber):
    """
    Assumes myList is sorted. Returns index of closest value to myNumber.
    If two numbers are equally close, return the index of the smallest number.

    Args:
        myList  (list): List of numbers to compare against. This should be sorted.
        myNumber (float): The number to compare against myList.

    Returns:
        pos (int): Index of position of number in myList which is closest to myNumber.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return (len(myList) - 1)
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return pos - 1

def volumes_from_grid( b, c, grid ):
    """
    Calculates the volume of each grid point based on the cell parameters and distance between grid points ( 1/2 ( deltax1 + deltax2 ) ).

    Args:
        b, c (float): Cell parameters in y and z directions.
        grid (np.array): 1D grid of ordered numbers over a region. 

    Returns:
        volumes (np.array): Volumes of each grid point on a 1D grid. 
    """
    return delta_x_from_grid( grid ) * b * c

def delta_x_from_grid( grid ):
    delta_x = np.zeros_like( grid )
    delta_x = ( grid[2:] - grid[:-2] ) / 2.0
    delta_x = np.insert( delta_x, (0, len( delta_x ) ), ( grid[1] - grid[0], grid[-1] - grid[-2] ) )
    return delta_x

class Grid_Point:
    """ The Grid_Point class contains the information and calculations for each grid point individually """

    def __init__( self, x, volume ):
        self.x = x
        self.volume = volume
        self.sites = []

    def average_site_energy( self, method = 'mean' ):
        """ 
   
        Returns the average segregation energy for all sites based on a specified method 

        Args: 
            method (str): The method in which the average segregation energies will be calculated.
                          'mean' - Returns the sum of all values at that site divided by the number of values at that site.
                          'min' - Returns the minimum segregation energy value for that site (appropriate for low temperature calculations).

        Returns:
            average site energies (np.array): Average segregation energies on a 1D grid.
  
        """

        if self.sites:
            return avg( np.array( [ s.energies() for s in self.sites ] ), method )
        else:
            return [None]

def avg( energies, method='mean' ):
    """ 
    Returns the average segregation energy for a site based on a specified method 

    Args: 
        energies (np.array): Segregation energies on 1D grid.  
        method (str): The method in which the average segregation energies will be calculated.
                      'mean' - Returns the sum of all values at that site divided by the number of values at that site.
                      'min' - Returns the minimum segregation energy value for that site (appropriate for low temperature calculations).

    Returns:
        average site energies (np.array): Average segregation energies on a 1D grid.
  
    """
    if method == 'mean':
        return [ np.mean( row ) for row in energies.T ] 
    elif method == 'min':
        return [ np.min( row ) for row in energies.T ]
    else:
        raise ValueError( "method: {}".format( method ) )

class Grid:
    def __init__( self, x_coordinates, b, c, site_set ):
       # x_coordinates need to be sorted for the delta_x calculation in volumes_from_grid
        """ 
        Grid objects contain information and methods for the calculation grid.
  
        Args:
            x_coordinates (np.array): The x-coordinates for each grid point.
            b (float):                b dimension for every grid-point.
            c (float):                c dimension for every grid-point.
            site_set (project.Set_of_Sites): The set of defect sites for populating the grid.
       
        Returns:
            None
        """
        self.volumes = volumes_from_grid( b, c, x_coordinates )
        self.points = [ Grid_Point( x, v ) for x, v in zip( x_coordinates, self.volumes ) ]
        self.x = x_coordinates
        for site in site_set:
            i = index_of_grid_at_x( self.x, site.x )
            self.points[ i ].sites.append( site ) 
            site.grid_point = self.points[i]

    def __getitem__( self, key ):
        return self.points[ key ]

    def charge( self, phi, temp ):
        """
        Calculates the overall charge at each point on a grid.
        
        Args:
            phi (np.array): Electrostatic potential on a 1D grid. 
            temp (float): Temperature of calculation.

        Returns:
            charge (np.array): Overall charge at each point on a 1D grid. 
        """
        charge = np.zeros_like( self.x )
        for i, point in enumerate( self ):
            charge[ i ] = sum( [ site.charge( phi_at_x( phi, self.x, site.x ), temp ) for site in point.sites ] )
        return charge

    def rho( self, phi, temp ):
        """
        Calculates charge density at each point on a grid, by dividing the overall charge at that grid point by the grid point volume.

        Args:
            phi (np.array): Electrostatic potential on a 1D grid.
            temp (float): Temperature of calculation.

        Returns:
            rho (np.array): The charge density at each point on a grid.
        """
        rho = np.zeros_like( self.x )
        rho = self.charge( phi, temp) / self.volumes
        return rho

    def defect_valences( self ):
        """ Returns an array of valences for each defect from self.defects """
        for site in self.points.sites:
            return np.array( [ d.valence for d in site.defects ] )

    def resistivity_ratio( self, defect_density, bulk_density ):
        """
        Calculates the grain boundary resistivity. This is the resistivity in the space charge region / the resistivity in the bulk crystal.

        Args:
            defect_density (array) : The defect density at each site.
            bulk_density (float) : The defect density in the bulk region.

        Returns:
            local_resistivity (float) : The grain boundary resistivity.  
        """
        grain_boundary_resistance = sum( delta_x_from_grid( self.x ) / defect_density ) 
        local_resistivity = ( bulk_density / sum( delta_x_from_grid( self.x ) ) ) * grain_boundary_resistance 
        return local_resistivity


    def average_site_energies( self, method = 'mean' ):
        """ 
        Returns the average segregation energy for all grid points based on a specified method 

        Args: 
            method (str): The method in which the average segregation energies will be calculated.
                          'mean' - Returns the sum of all values at that site divided by the number of values at that site.
                          'min' - Returns the minimum segregation energy value for that site (appropriate for low temperature calculations).

        Returns:
            average site energies (np.array): Average segregation energies on a 1D grid.
        """
        return np.array( [ p.average_site_energy( method ) for p in self.points ] ).T

    def interpolated_energies( self ):
        """Returns the average site energies linearly interpolated onto a regularly spaced grid"""

        energies = []
        for row in self.average_site_energies():
            x_new, e_new = np.array( [ ( x, e ) for x, e in zip( self.x, row ) if e ] ).T
            interpolated_energies = griddata( x_new, e_new, self.x, fill_value = 0.0 )
            energies.append( interpolated_energies )
