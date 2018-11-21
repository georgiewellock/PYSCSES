import numpy as np
from project.grid import Grid
from scipy.interpolate import griddata
import math
from project.set_up_calculation import site_from_input_file, load_site_data, mirror_site_data
from project.grid import index_of_grid_at_x, phi_at_x, energy_at_x
from project.site import Site
from project.constants import boltzmann_eV
from bisect import bisect_left

class Set_of_Sites:
    """ The Set_of_Sites object groups together all of the Site objects into one object and contains functions for the calculations that provide properties of all of the sites together rather than individually. """
    def __init__( self, sites ):
        self.sites = sites
    
    def __add__( self, other ):
        """ Allows the concatenation of multiple Set_of_Sites objects"""
        if type( other ) is not Set_of_Sites:
            raise TypeError
        return Set_of_Sites( self.sites + other.sites )

    def __getitem__( self, index ):
        """ Returns the site corresponding to a given index """
        return self.sites[ index ]

    def __iter__( self ):
        """
        Iterator over self.sites
        """
        return iter( self.sites )

    def subset( self, label ):
        """ Returns a list of all the sites which contain a particular defect """
        return [ s for s in self.sites if s.label == label ] 

    def get_coords( self, label ):
        """ Returns a list of the x coordinates for all the sites wich contain a particular defect """
        return [ s.x for s in self.sites if s.label == label ]

    def calculate_energies_on_grid( self, grid, phi ):
        """ 
        Returns an array of energies at their points on a one-dimensional grid.
    
        Args: 
            grid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            phi (array): electrostatic potential on a one-dimensional grid. 

        Returns:
            energies_on_grid (array): energies at their grid points
 
        """
        energies_on_grid = np.zeros_like( grid )
        for site in self.sites:
            energies_on_grid[index_of_grid_at_x( grid.x, site.x )] =+ energy_at_x( site.defect_energies, grid.x, site.x )
        return energies_on_grid

    def calculate_probabilities( self, grid, phi, temp ):
        """ 
        Calculates the probability of a site being occupied by its corresponding defect.
    
        Args: 
            grid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            phi (array): electrostatic potential on a one-dimensional grid. 
            temp (float): Absolute temperature.

        Returns:
            probability (array): probabilities of defects occupying each site using their grid points
 
        """
        probability = np.zeros_like( grid.x )
        for site in self.sites:
            probability[index_of_grid_at_x( grid.x, site.x )] = np.asarray( site.probabilities( phi_at_x( phi, grid.x, site.x ), temp ) )
        return probability 

    def calculate_defect_density( self, grid, phi, temp ):
        """ 
        Calculates the defect density at each site.
    
        Args: 
            grid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            phi (array): electrostatic potential on a one-dimensional grid. 
            temp (float): Absolute temperature.

        Returns:
            defect_density (array): defect density for each site using their grid points
 
        """
        defect_density = np.zeros_like( grid.x )
        for site in self.sites:
            i = index_of_grid_at_x( grid.x, site.x )
            defect_density[ i ] += np.asarray( site.probabilities( phi_at_x( phi, grid.x, site.x ), temp ) ) / grid.volumes[ i ]
        return defect_density  

    def subgrid_calculate_defect_density( self, sub_grid, full_grid, phi, temp ):
        """ 
        Calculates the defect density at each site for a given subset of sites.
    
        Args: 
            subgrid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates. For a given subset of sites. 
            full_grid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates. For all sites.
            phi (array): electrostatic potential on a one-dimensional grid. 
            temp (float): Absolute temperature.

        Returns:
            defect_density (array): defect density for each site using their grid points
 
        """
        defect_density = np.zeros_like( sub_grid.x )
        for site in self.sites:
            i = index_of_grid_at_x( sub_grid.x, site.x )
            defect_density[ i ] += np.asarray( site.probabilities( phi_at_x( phi, full_grid.x, site.x ), temp ) ) / sub_grid.volumes[ i ]
        return defect_density  


    def form_continuum_sites( all_sites, x_min, x_max, n_points, b, c, defect_species, limits_for_laplacian, site_labels, defect_labels ):
        """
        NEEDS UPDATING TO WORK WITH NEW CODE FORMAT.
        INCLUDE DOCSTRING
        """
        grid_1 = np.linspace( x_min, x_max, n_points )
        limits = [ grid_1[1]-grid_1[0], grid_1[1]-grid_1[0] ]
        Gd_scaling = len( all_sites.subset( site_labels[1] ) ) / len( grid_1 )
        Vo_scaling = len( all_sites.subset( site_labels[0] ) ) / len( grid_1 )
    
        Vo_continuum_grid = Grid( grid_1, b, c, limits, limits_for_laplacian, all_sites.subset( site_labels[0] ) )
        Gd_continuum_grid = Grid( grid_1, b, c, limits, limits_for_laplacian, all_sites.subset( site_labels[1] ) )
    
        Vo_average_energies = np.array( [ site.average_local_energy( method = 'mean' )[0] for site in all_sites.subset( site_labels[0] ) ] )
        Gd_average_energies = np.array( [ site.average_local_energy( method = 'mean' )[0] for site in all_sites.subset( site_labels[1] ) ] )
    
    
        Vo_new_energies = griddata( ( [ site.x for site in all_sites.subset( site_labels[0] ) ] ), Vo_average_energies, grid_1, method = 'nearest' )
        Gd_new_energies = griddata( ( [ site.x for site in all_sites.subset( site_labels[1] ) ] ), Gd_average_energies, grid_1, method = 'nearest' )
    
    
        Vo_new_sites = Set_of_Sites( [ Site( site_labels[0], x, [ defect_species[defect_labels[0]] ], [e], scaling = np.array( Vo_scaling ) ) for x, e in zip( grid_1, Vo_new_energies ) ] )
        Gd_new_sites = Set_of_Sites( [ Site( site_labels[1], x, [ defect_species[defect_labels[1]] ], [e], scaling = np.array( Gd_scaling ) ) for x, e in zip( grid_1, Gd_new_energies ) ] )   
    
        all_sites_new = Vo_new_sites + Gd_new_sites
    
        return all_sites_new, limits

def form_continuum_sites_new( all_sites, x_min, x_max, n_points, b, c, defect_species, limits_for_laplacian, site_labels, defect_labels ):
    
    limits = [ x_min, x_max ]
    grid = np.linspace( x_min, x_max, n_points )
    sites = []
    for label, d_label in zip(site_labels, defect_labels):
        scaling = len( all_sites.subset( label ) ) / len( grid )
        continuum_grid = Grid( grid, b, c, limits, limits_for_laplacian, all_sites.subset( label ) )
        average_energies = np.array( [ site.average_local_energy( method = 'mean' )[0] for site in all_sites.subset( label ) ] )
        new_energies = griddata( ( [ site.x for site in all_sites.subset( label ) ] ), Vo.average_energies, grid, method = 'nearest' )
        for x, e in zip( grid, new_energies:
            sites.append( Site( label, x, [ defect_species[ d_label ] ], [e], scaling = np.array( scaling ) ) )
    return Set_of_Sites( sites )    


    @ classmethod
    def set_of_sites_from_input_data( cls, input_data, limits, defect_species, site_charge, core, temperature, offset=0.0 ):
        site_data = load_site_data( input_data, limits[0], limits[1], site_charge, offset )
        energies = [ line[4] for line in site_data ]
        min_energy = min(energies)
        if core == 'single':
            for line in site_data:
                if line[4] > min_energy:
                    line[4] = 0.0
        #print(boltzmann_eV * temperature, flush=True)
        if core == 'multi_site':
            for line in site_data:
                if ( -boltzmann_eV * temperature) <= line[4] <= ( boltzmann_eV * temperature ):
                    line[4] = 0.0
        return Set_of_Sites( [ site_from_input_file( line, defect_species, site_charge, core, temperature ) for line in site_data ] )

    @ classmethod
    def mirrored_set_of_sites_from_input_data( cls, input_data, limits, defect_species ):
        site_data = load_site_data( input_data, limits[0], limits[1] )
        mirrored_data = mirror_site_data( site_data )


    @ classmethod
    def core_width_analysis( cls, input_data, limits, defect_species, site_charge, core, temperature ):
        site_data = load_site_data( input_data, limits[0], limits[1], site_charge )
        energies = [ line[4] for line in site_data ]
        min_energy = min(energies)
        if core == 'single':
            for line in site_data:
                if line[4] > min_energy:
                    line[4] = 0.0
        #print(boltzmann_eV * temperature, flush=True)
        if core == 'multi_site':
            for line in site_data:
                if ( -boltzmann_eV * temperature) <= line[4] <= ( boltzmann_eV * temperature ):
                    line[4] = 0.0
        energies = [line[4] for line in site_data ]
        x = [line[2] for line in site_data ]
        x_seg = np.column_stack(( x, energies ))
        minval = np.min(x_seg[:,0][np.nonzero(x_seg[:,1])])
        maxval = np.max(x_seg[:,0][np.nonzero(x_seg[:,1])])
        core_width = maxval-minval
        return(core_width)
