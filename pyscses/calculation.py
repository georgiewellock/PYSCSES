import numpy as np 
import mpmath
from bisect import bisect_left, bisect_right
from pyscses.set_of_sites import SetOfSites
from pyscses.matrix_solver import MatrixSolver
from pyscses.set_up_calculation import calculate_grid_offsets
from pyscses.constants import *
from pyscses.grid import delta_x_from_grid, Grid, phi_at_x
from scipy.optimize import minimize

class Calculation:
    """The Calculation class contains methods for calculating the relevant space charge properties for any given system, such as electrostatic potential, charge density, defect mole fractions and parallel and perpendicular grain boundary resistivities. 

    Args:
        grid (:obj:`pyscses.Grid`): A pyscses.Grid object. This contains properties of the grid including the x coordinates and the volumes.
        bulk_x_min (float): The minimum x coordinate defining a region of bulk.
        bulk_x_max (float): The maximum x coordinate defining a region of bulk.
        alpha (float): A damping parameter for updating the potential at each iteration.
        convergence (float): The convergence limit. The difference between the updated phi and the phi from the previous iteration must be below this for the calculation to be sufficently converged. 
        dielectric (float): The dielectric constant for the studied material.
        temp (float): The temperature that the calculation is run. 
        boundary_conditions (str): Specified boundary conditions for the matrix solver. Allowed values are `dirichlet` and `periodic`. Default = `dirichlet`.
 
    """

    def __init__(self, grid, bulk_x_min, bulk_x_max, alpha, convergence, dielectric, temp, boundary_conditions ):
        self.grid = grid
        self.bulk_x_min = bulk_x_min
        self.bulk_x_max = bulk_x_max
        self.alpha = alpha
        self.convergence = convergence
        self.dielectric = dielectric 
        self.temp = temp 
        self.boundary_conditions = boundary_conditions


    def mole_fraction_error( self, input_mole_fractions, target_mole_fractions, approximation):
        """Calculates the sum of sqaures error between the output mole fractions calculated from the input mole fractions, compared to the target mole fractions.

        Args:
            input_mole_fractions (list): Mole fractions for each of the species used in the iterative Poisson-Boltzmann solver.
            target_mole_fractions  (list): The value that the mole fractions should be in the bulk. 

        Returns:
            float: Sum of squares error between output and target mole fractions. 

        """
        input_mole_fractions = np.array([input_mole_fractions])
        output_mole_fractions = self.mole_fraction_output(input_mole_fractions, approximation)
        squares = []
        for mf in input_mole_fractions:
            for i in range(len(mf)):
                output = output_mole_fractions[0,i]
                target = target_mole_fractions[0,i]
                squares.append(( output - target )**2)
        return sum( squares )

    def mole_fraction_output(self, input_mole_fractions, approximation):
        """Calculates the output mole fraction for a given input mole fraction when solving the Poisson-Boltzmann equation.

        Args:
            input_mole_fractions (list): Mole fractions for each of the species used in the iterative Poisson-Boltzmann solver.
            approximation (str): The defect mobility approximation. Either 'mott-schottky' to enforce only a single mobile defect, or 'gouy-chapman' to allow all defect species to redistribute. 
   
        Returns:
            list: Mole fractions that are calculated from the iterative Poisson-Boltzmann solver.
        """
        for mf in input_mole_fractions:
            for i in range(len(mf)):
                for site in self.grid.set_of_sites.subset( self.site_labels[i] ):
                    for defect in site.defect_species:
                        defect.mole_fraction = input_mole_fractions[0,i]
                    for defect in site.defects:
                        defect.mole_fraction = input_mole_fractions[0,i]
        
        self.solve(approximation)
        species = []
        for mf in input_mole_fractions:
            for i in range(len(mf)):
                species.append(self.site_labels[i])
        self.form_subgrids( species )
        self.mole_fractions()
        average_mole_fractions = []
        for mf in input_mole_fractions:
            for i in range(len(mf)):
                self.mf[self.site_labels[i]] = [mf for mf in self.mf[self.site_labels[i]] if mf != 0.0 ]
                average_mf = self.calculate_average(self.subgrids[self.site_labels[i]], self.bulk_x_min, self.bulk_x_max, self.mf[self.site_labels[i]])
                average_mole_fractions.append(average_mf) 
        output_mole_fractions = np.array( [ average_mole_fractions ] )
        return output_mole_fractions

        
    def mole_fraction_correction( self, target_mole_fractions, approximation, initial_guess ):
        """Starting from an initial guess for the appropriate input mole fractions, minimises the error between the target bulk mole fraction and the output mole fraction from the iterative Poisson-Boltzmann solver. The output is stored as a Calculation attribute. Calculation.initial_guess (list): The optimum values to be used as the input mole fractions for the iterative Poisson-Boltzmann solver so that the output bulk mole fractions are the target bulk mole fractions. 

        Args:
            target_mole_fractions (list): The value that the mole fractions should be in the bulk. 
            approximation (str): The defect mobility approximation. Either 'mott-schottky' to enforce only a single mobile defect, or 'gouy-chapman' to allow all defect species to redistribute. 
            initial_guess (list): Values for an initial guess for the defect mole fractions used in the error minimisation. 

        """
        self.initial_guess = initial_guess
        bounds = [(0.0001, 1)]*len(target_mole_fractions)
        target_mole_fractions = np.array([target_mole_fractions])
        opt_mole_fractions = minimize( self.mole_fraction_error, self.initial_guess, args=(  target_mole_fractions, approximation ), bounds = (bounds) )
        opt_mole_fractions.x = np.array([opt_mole_fractions.x])
        for mf in target_mole_fractions:
            for i in range(len(mf)):
                for site in self.grid.set_of_sites.subset( self.site_labels[i] ):
                    for defect in site.defect_species:
                        defect.mole_fraction = opt_mole_fractions.x[0,i]
                    for defect in site.defects:
                        defect.mole_fraction = opt_mole_fractions.x[0,i]
        self.initial_guess = opt_mole_fractions.x

    def find_index( self, grid, min_cut_off, max_cut_off ):
        """Calculates the indices of the grid positions closest to a minimum and maximum value.

        Args:
	    grid (:obj:`pyscses.Grid`): pyscses.Grid object. This contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
	    min_cut_off (float): Minimum x coordinate value defining the calculation region. 
	    max_cut_off (float): Maximum x coordinate value defining the calculation region.

	Returns:
	    int, int: Index for minimum cutoff; index for maximum cutoff

	"""
        min_index = bisect_left( grid.x, min_cut_off )
        max_index = bisect_left( grid.x, max_cut_off )
        return min_index, max_index

    def calculate_offset( self, grid, min_cut_off, max_cut_off ):
        """Calculate the offset between the midpoint of the last x coordinate in the calculation region and the x coordinate outside of the calulation region.

        Args:
            grid (:obj:`pyscses.Grid`): Contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            min_cut_off (float): Minimum x coordinate value defining the calculation region. 
	    max_cut_off (float): Maximum x coordinate value defining the calculation region.

        Returns:
            float: Offset for the minimum x coordinate.
	    float: Offset fot the maximum x coordinate.

        """
        min_index, max_index = self.find_index( grid, min_cut_off, max_cut_off )
        min_offset = ( ( grid.x[ min_index + 1 ] - grid.x[ min_index ] ) / 2 ) + ( ( grid.x[ min_index ] - grid.x[ min_index - 1 ] ) / 2 )
        max_offset = ( ( grid.x[ max_index ] - grid.x[ max_index - 1 ] ) / 2 ) + ( ( grid.x[ max_index - 1 ] - grid.x[ max_index -2 ] ) / 2 )
        return min_offset, max_offset

    def calculate_delta_x( self, grid, min_cut_off, max_cut_off ):
        """Calculates the distance between the midpoints of each consecutive site. Inserts the calculated distance to the next grid point outside of the calculation region to the first and last position as the delta_x value for the endmost sites. 

        Args:
            grid (:obj:`pyscses.Grid`): Contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            min_cut_off (float): Minimum x coordinate value defining the calculation region. 
	    max_cut_off (float): Maximum x coordinate value defining the calculation region.

	Returns:
	    list: Distance between consecutive sites. 

	"""
        min_index, max_index = self.find_index( grid, min_cut_off, max_cut_off )
        min_offset, max_offset = self.calculate_offset( grid, min_cut_off, max_cut_off )
        return delta_x_from_grid( grid.x[ min_index+1 : max_index ], [min_offset, max_offset] )

    def calculate_average( self, grid, min_cut_off, max_cut_off, sc_property ):
        """Calculate the average of a given space chage property over a given region.

        Args:
            grid (:obj:`pyscses.Grid`): Contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            min_cut_off (float): Minimum x coordinate value defining the calculation region. 
            max_cut_off (float): Maximum x coordinate value defining the calculation region.
            sc_property (list): Value of space charge property at all sites.

        Returns:
            float: The average value for the property over the given sites. 

	"""
        min_index, max_index = self.find_index( grid, min_cut_off, max_cut_off )
        delta_x = self.calculate_delta_x( grid, min_cut_off, max_cut_off )
        return np.sum( sc_property[ min_index + 1 : max_index ] * delta_x ) / np.sum( delta_x )

    def solve( self, approximation ):
        """
        Self-consistent solving of the Poisson-Boltzmann equation. Iterates until the convergence is less than the convergence limit. The outputs are stored as Calculation attributes.
        Calculation.phi (array): Electrostatic potential on a one-dimensional grid. 
        Calculation.rho (array): Charge density on a one-dimensional grid.
        Calculation.niter (int): Number of iterations performed to reach convergence.


        Args:
            approximation (str): Approximation used for the defect behaviour.
                                 'mott-schottky' - Some defects immobile / fixed to bulk mole fractions.
                                 'gouy-chapman' - All defects mobile / able to redistribute.

        """
        poisson_solver = MatrixSolver( self.grid, self.dielectric, self.temp, boundary_conditions=self.boundary_conditions )

        phi = np.zeros_like( self.grid.x )
        rho = np.zeros_like( self.grid.x )

        conv = 1
        niter = 0
        while conv > self.convergence:
            predicted_phi, rho = poisson_solver.solve( phi )
            if approximation == 'gouy-chapman':
                average_phi = self.calculate_average( self.grid, self.bulk_x_min, self.bulk_x_max, predicted_phi )
                predicted_phi -= average_phi
            if approximation == 'mott-schottky':
                vo_predicted_phi = [ phi_at_x( predicted_phi, self.grid.x, x ) for x in self.grid.subgrid(self.site_labels[0]).x ]
                average_vo_predicted_phi = self.calculate_average( self.grid.subgrid( self.site_labels[0] ), self.bulk_x_min, self.bulk_x_max, vo_predicted_phi )
                predicted_phi -= average_vo_predicted_phi
            phi =  self.alpha * predicted_phi + ( 1.0 - self.alpha ) * phi
            conv = sum( ( predicted_phi - phi )**2) / len( self.grid.x )
            prob = self.grid.set_of_sites.calculate_probabilities( self.grid, phi, self.temp )
            niter += 1
#            if niter % 500 == 0.0:
#            if niter == 1:
#                print(conv)
#                print(phi, rho)
#                stop
        self.phi = phi
        self.rho = self.grid.rho( phi, self.temp )
        self.niter = niter

    def form_subgrids( self, site_labels ):
        """Creates a `pysces.Grid` object for each species in the system. The output is a dictionary of separate Grid classes for the different site species and is stored as Calculation.subgrids. 

        Args:
            site_labels (list): List of strings for the different site species.

        """
        self.site_labels = site_labels
        subgrids = {} 
        for label in site_labels:
            name = '{}'.format( label ) 
            subgrids[name] = self.grid.subgrid( label )
            subgrids[name].delta_x[0] = subgrids[name].delta_x[1]
            subgrids[name].delta_x[-1] = subgrids[name].delta_x[1]
            subgrids[name].volumes[0] = subgrids[name].volumes[1]
            subgrids[name].volumes[-1] = subgrids[name].volumes[1]
        self.subgrids = subgrids

    def create_subregion_sites( self, grid, min_cut_off, max_cut_off ):
        """Creates a `pyscses.SetOfSites` object for a defined region of the grid.

        Args:
            grid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            min_cut_off (float): Minimum x coordinate value defining the calculation region. 
            max_cut_off (float): Maximum x coordinate value defining the calculation region.

        Returns:
            :obj:`pyscses.SetOfSites`: Set of Sites object for a given subregion of the grid.

        """
        sites = []
        for site in grid.set_of_sites:
            if site.x > min_cut_off and site.x < max_cut_off:
                sites.append( site )
        sites = SetOfSites( sites )
        return sites

    def create_space_charge_region( self, grid, pos_or_neg_scr, scr_limit ):
        """Calculate the space charge region. The space charge region is defined as the region when the electrostatic potential is greater than a predefined limit.

	Args:
	    grid (:obj:`pyscses.Grid`): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
	    pos_or_neg_scr (str): 'positive' - for a positive space charge potential.
				  'negative' - for a negative space charge potential.
	    scr_limit (float): The minimum electrostatic potential that the electrostatic potential must exceed to be included in the space charge region.

	Returns:
	    list: List of x coordinates for sites within the space charge region.

	"""	
        space_charge_region = []
        self.phi_on_mobile_defect_grid = [ phi_at_x( self.phi, self.grid.x, x ) for x in grid.x ] 
        x_and_phi = np.column_stack( ( grid.x, self.phi_on_mobile_defect_grid ) )
        for i in range( len( x_and_phi ) ):
            if pos_or_neg_scr == 'positive':
                if x_and_phi[i, 1]-x_and_phi[0,1] > scr_limit:
                    space_charge_region.append( x_and_phi[i,0] )
            if pos_or_neg_scr == 'negative': 
                if x_and_phi[i,1]-x_and_phi[0,1] < scr_limit:
                    space_charge_region.append( x_and_phi[i,0] )
        return space_charge_region

    def calculate_mobile_defect_conductivities( self, pos_or_neg_scr, scr_limit, species, mobility_scaling=False ):
        """Calculate the conductivity ratio between the space charge region and the bulk both perpendicular and parallel to the grain boundary.

        A `SetOfSites` object is created for the sites in the space charge region, and the defect distributions calculated. The width of the space charge region is calculated and a bulk region of the same width is defined. A SetOfSites object for the bulk region is created and the defect distributions calculated. Taking each site as a resistor in series or parallel respectively, the conductivity is calculated and the ratio between the space charge region and the bulk is taken. 

        Args:
            pos_or_neg_scr (str): 'positive' - for a positive space charge potential.
				  'negative' - for a negative space charge potential.
            scr_limit (float): The minimum electrostatic potential that the electrostatic potential must exceed to be included in the space charge region.
            species (str): The species for which the conductivity is being calculated.
            mobility_scaling (bool): For particles on a lattice which only interact through volume exclusion, the mobility exhibits a blocking term. True if the blocking term is to be included, False if the blocking term is not to be included. Default = False.

        Returns:
            float: The perpendicular conductivity ratio. The conductivity ratio between the bulk and the space charge region perpendicular to the grain boundary.
            float: The parallel conductivity ratio. The conductivity ratio between the bulk and the space charge region parallel to the grain boundary. 

	"""
        space_charge_region = self.create_space_charge_region( self.subgrids[species], pos_or_neg_scr, scr_limit )
        space_charge_region_limits = self.calculate_offset( self.subgrids[species], np.min(space_charge_region), np.max(space_charge_region) )
        space_charge_region_sites = self.create_subregion_sites( self.subgrids[species], np.min(space_charge_region), np.max(space_charge_region) )
        for site in space_charge_region_sites:
            charge = site.defects[0].valence
            mobilities = site.defects[0].mobility
        space_charge_region_grid = Grid.grid_from_set_of_sites( space_charge_region_sites, space_charge_region_limits, space_charge_region_limits, self.grid.b, self.grid.c )
        space_charge_region_width = space_charge_region_grid.x[-1] - space_charge_region_grid.x[0]
        mobile_defect_density = SetOfSites( self.subgrids[species].set_of_sites ).subgrid_calculate_defect_density( self.subgrids[species], self.grid, self.phi, self.temp )
        space_charge_region_mobile_defect_mf = space_charge_region_sites.calculate_probabilities( space_charge_region_grid, self.phi, self.temp )
        space_charge_region_mobile_defect_density = space_charge_region_sites.subgrid_calculate_defect_density( space_charge_region_grid, self.grid, self.phi, self.temp )
        if mobility_scaling:
             mobile_defect_conductivity = space_charge_region_mobile_defect_density * ( 1 - space_charge_region_mobile_defect_mf ) * charge * mobilities 
        else:
            mobile_defect_conductivity = space_charge_region_mobile_defect_density * charge * mobilities
        bulk_x_max = self.bulk_x_min + space_charge_region_width
        min_bulk_index, max_bulk_index = self.find_index( self.subgrids[species], self.bulk_x_min, bulk_x_max )
        self.bulk_limits = self.calculate_offset( self.subgrids[species], self.bulk_x_min, bulk_x_max ) 
        bulk_mobile_defect_sites = self.create_subregion_sites( self.subgrids[species], self.bulk_x_min, bulk_x_max )
        bulk_mobile_defect_grid = Grid.grid_from_set_of_sites( bulk_mobile_defect_sites, self.bulk_limits, self.bulk_limits, self.grid.b, self.grid.c )
        bulk_mobile_defect_density = SetOfSites( bulk_mobile_defect_grid.set_of_sites ).subgrid_calculate_defect_density( bulk_mobile_defect_grid, self.grid, self.phi, self.temp )  
        bulk_region_mobile_defect_mf = bulk_mobile_defect_sites.calculate_probabilities( bulk_mobile_defect_grid, self.phi, self.temp )
        if mobility_scaling:
            bulk_mobile_defect_conductivity = bulk_mobile_defect_density * charge * mobilities
        else:
            bulk_mobile_defect_conductivity = bulk_mobile_defect_density * charge * mobilities * (1-bulk_region_mobile_defect_mf)
        space_charge_array = np.column_stack( ( mobile_defect_conductivity, space_charge_region_grid.x ) )
        bulk_array = np.column_stack( ( bulk_mobile_defect_conductivity, bulk_mobile_defect_grid.x ) )
        if mobilities != 0.0:
            space_charge_perpendicular = sum( space_charge_region_grid.delta_x / mobile_defect_conductivity )
            self.average_bulk_mobile_defect_density = sum(bulk_mobile_defect_grid.delta_x * bulk_mobile_defect_density ) / sum(bulk_mobile_defect_grid.delta_x)
            bulk_perpendicular = sum( bulk_mobile_defect_conductivity / bulk_mobile_defect_grid.delta_x )
            space_charge_parallel = sum( mobile_defect_conductivity / space_charge_region_grid.delta_x )
            bulk_parallel = sum( bulk_mobile_defect_grid.delta_x / bulk_mobile_defect_conductivity )
            perpendicular_conductivity_ratio = 1 / ( space_charge_perpendicular * bulk_perpendicular )
            parallel_conductivity_ratio = space_charge_parallel * bulk_parallel
        else:
            perpendicular_conductivity_ratio = 0.0
            parallel_conductivity_ratio = 0.0 
#        self.depletion_factor = 1 - ( mobile_defect_density / average_bulk )
        return perpendicular_conductivity_ratio, parallel_conductivity_ratio

    def calculate_resistivity_ratio( self, pos_or_neg_scr, scr_limit, mobility_scaling=False ):
        """
        Each species present in the system is looped over and the perpendicular and parallel conductivity ratios are calculated and appended to separate lists. Once all species have been added to the list, the conductivities are summed and the reciprocal taken to give the resistivity ratios perpendicular and parallel to the grain boundary. The outputs are stored as calculation attributes, Calculation.perpendicular_resistivity_ratio(float), Calculation.parallel_resistivity_ratio (float): :math:`r_\mathrm{gb}^\perp, r_\mathrm{gb}^\parallel`: The perpendicular and parallel grain boundary resistivity ratios.
        Args:
            pos_or_neg_scr (str): 'positive' - for a positive space charge potential.
				  'negative' - for a negative space charge potential.
            scr_limit (float): The minimum electrostatic potential that the electrostatic potential must exceed to be included in the space charge region.
            mobility_scaling (bool): For particles on a lattice which only interact through volume exclusion, the mobility exhibits a blocking term. True if the blocking term is to be included, False if the blocking term is not to be included. Default = False.

	"""
        full_parallel_conductivity_data = []
        full_perpendicular_conductivity_data = []
        for label in self.site_labels:
            c_per, c_par = ( self.calculate_mobile_defect_conductivities( pos_or_neg_scr, scr_limit, label, mobility_scaling  ))
            full_parallel_conductivity_data.append(c_par)
            full_perpendicular_conductivity_data.append(c_per)
        self.perpendicular_resistivity_ratio = 1 / sum(full_perpendicular_conductivity_data)
        self.parallel_resistivity_ratio = 1 / sum(full_parallel_conductivity_data)
            
    def solve_MS_approx_for_phi( self, valence ):
        """Calculate the space-charge potential, :math:`\phi_0`, from the grain-boundary resistivity ratio, within the Mott-Schottky approximation.
        Within the Mott-Schottky approximation the grain boundary resistivity is related to the space-charge potential (the electrostatic potential at the grain boundary core, compared to the bulk value) according to
        
        .. math:: r_\mathrm{gb} = \frac{\rho_{i,\mathrm{gb}}}{\rho_{i,\infty}} = \frac{\exp(z_i\phi_0 / V_\mathrm{th})}{2z_i\phi_0/V_\mathrm{th}}
       
        where

        .. math:: V_\mathrm{th} = \frac{k_\mathrm{B}T}{q}.

        (See e.g. `S. Kim, Phys. Chem. Chem. Phys. 18, 19787 (2016).`_)

        .. _S. Kim, Phys. Chem. Chem. Phys. 18, 19787 (2016).: https://dx.doi.org/10.1039/c6cp02177h

        This allows a Mott-Schottky space-charge potential, :math:`\phi_{0,\mathrm{MS}}`, to be calculated using the `LambertW`_ function:

        .. math:: \phi_{0,\mathrm{MS}} = -\mathrm{LambertW}\left(\frac{1}{2 r_\mathrm{gb}}\right)\frac{V_\mathrm{th}}{z_i}.
        .. _LambertW: https://en.wikipedia.org/wiki/Lambert_W_function

        The output is stored as a Calculation attribute. Calculation.ms_phi (float): :math:`\phi_{0,\mathrm{MS}}`. The space charge potential calculated from Mott-Schottky model.


        Args:
            valence( float ): Charge of the mobile defect species. 
   
        Raises:
            ValueError: If the calculated resistivity ratio is less than 1.36, the LambertW function returns a complex, non-physical value.

        """
        if self.perpendicular_resistivity_ratio < 1.36:
            raise ValueError( "Resistivity ratio < 1.36. Solution not on a real branch." )
        self.ms_phi = (-mpmath.lambertw(-1/(2*self.perpendicular_resistivity_ratio),k=-1)) * ( ( boltzmann_eV * self.temp ) / valence )

    def calculate_debye_length( self ):
        """Calculate the `Debye length`_.

        .. _Debye length: https://en.wikipedia.org/wiki/Debye_length

        The output is stored as a Calculation attribute. Calculation.debye_length (float): The Debye length as derived from Poisson-Boltzmann equation.

        Args:
            None

        """
        self.debye_length = np.sqrt( ( self.dielectric * vacuum_permittivity * boltzmann_eV * self.temp ) / ( 2 * ( fundamental_charge ** 2 ) * self.average_bulk_mobile_defect_density ) )

    def calculate_space_charge_width( self, valence ):
        """Calculate the approximate space charge width from the Debye length.
           The output is stores as a Calculation attribute. Calculation.space_charge_width (float): The approximate space charge width.

        Args:
	    valence (float): The charge of the mobile defect species.

        """
        self.space_charge_width = 2 * ( self.debye_length * math.sqrt( max(self.phi) / ( ( boltzmann_eV * self.temp ) / valence ) ) )


    def mole_fractions( self ):
        """Calculate the mole fractions (probability of defects occupation) for each site on the subgrid for each species. The output is stored as a Calculation attribute. Calculation.mf (dict): A dictionary of the defect species mole fractions for each site on the subgrid for each site species. 
        Args:
            None
 
        """
        mole_fractions = {}
        for label in self.site_labels:
            name = '{}'.format(label)
            mole_fractions[name] = SetOfSites(self.subgrids[name].set_of_sites).calculate_probabilities( self.grid, self.phi, self.temp )
        self.mf = mole_fractions

def diff_central(x, y):
    """Calculate the numerical derivative of x,y data using a central difference approximation.
   
    Args:
        x (numpy.array): x values.
        y (numpy.array): y values.

    Returns:
        numpy.array: numerical derivative of x,y

    """
    x0 = x[:-2]
    x1 = x[1:-1]
    x2 = x[2:]
    y0 = y[:-2]
    y1 = y[1:-1]
    y2 = y[2:]
    f = (x2 - x1)/(x2 - x0)
    return (1-f)*(y2 - y1)/(x2 - x1) + f*(y1 - y0)/(x1 - x0) 

def calculate_activation_energies( ratios, temp ):
    """Solves the Arrhenius equation using the calculated resistivity ratio for a series of temperatures to calculate the activation energy for single defect segregation.

    Uses a central difference approach so endpoints filled in with NaN to create an array with the same length as the args.
   
    Args:
        ratios (list): Resistivity ratios calculated for a range of temperatures.
        temp (list): Temperature (values used to calculate resistivity ratio values).

    Returns:
        numpy.array: The activation energy calculated for different temperatures.  

    """
    temp = np.array( temp )
    ratios = np.array(ratios)
    x = 1 / temp
    y = np.log( 1 / ratios )
    slopes = diff_central( x, y )
    Ea = slopes*-boltzmann_eV
    Ea = np.append( Ea, 0 )
    Ea = np.insert( Ea, [0], 0 )
    Ea2 = np.where(Ea == 0, np.nan, Ea )
    return Ea2           
