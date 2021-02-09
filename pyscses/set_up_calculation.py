from copy import copy
from operator import itemgetter
from pyscses.site import Site
import numpy as np
from pyscses.constants import boltzmann_eV
from bisect import bisect_left, bisect_right
from sklearn.cluster import AgglomerativeClustering # type: ignore

def site_from_input_file( site, defect_species, site_charge, core, temperature ):
    """
    Takes the data from the input file and converts it into a site.
    The input data file is a .txt file where each line in the file corresponds to a site. The values in each line are formatted and separated into the corresponding properties before creating a Site object for each site. 

    Args:
        site (str): A line in the input file.
        defect_species (object): Class object containing information about the defect species present in the system.
        site_charge (bool): The site charge refers to the contribution to the overall charge of a site given by the original, non-defective species present at that site. True if the site charge contribution is to be included in the calculation, False if it is not to be included.

    Returns:
        :obj:`Site`

    """
    label = site[0]
    if site_charge == True:
        valence = float(site[1])
    if site_charge == False:
        valence = 0.0
    x = float(site[2])
    defect_labels = site[3::2]
    defect_energies = [ float(e) for e in site[4::2] ]
    min_energy = min(defect_energies)
    if core == 'single':
        for d_e in defect_energies:
            if d_e > min_energy:
                d_e = 0.0
    if core == 'multi-site':
        for d_e in defect_energies:
            if ( -boltzmann_eV * temperature) <= d_e <= ( boltzmann_eV * temperature ):
                d_e = 0.0
    #defect_energies = [ 0.0 for e in site[4::2] ]
    return Site( label, x, [ defect_species[l] for l in defect_labels ], defect_energies, valence=valence )

def format_line( line, site_charge, offset = 0.0 ):
    """
    Each line in the input file is formatted into separate site properties. 

    Args:
        line (str): A line in the input file.
        site_charge (bool): The site charge refers to the contribution to the overall charge of a site given by the original, non-defective species present at that site. True if the site charge contribution is to be included in the calculation, False if it is not to be included.

    Returns: 
        list: line from the input file, split into individual values. 

    """
    # x coordinate
    if site_charge == True:
        line[1] = float( line[1] )
    if site_charge == False:
        line[1] = 0.0

    line[2] = float( line[2] )
    # defect energies
    for i in range( 4, len(line) ):
        line[i] = float( line[i] )
        if line[i] < 0.0:
            line[i] += offset
        if line[i] > 0.0:
            line[i] += offset
#        line[i] = 0.0
    return line

def load_site_data( filename, x_min, x_max, site_charge, offset = 0.0 ):
    """
    Reads in the input data and formats the input data if the x coordinate values are within the calculation region.

    Args:
        filename(string): Filename for importing the input data.
        x_min(float): Minimum x coordinate value defining the calculation region.
        x_max(float): Maximum x coordinate value defining the calculation region. 
        site_charge (bool): True if site charges are to be included in the calculation, false if they are not to be included.
   
    Return:
        list: formatted data for calculation.

    """
    with open( filename, 'r' ) as f:
        input_data = [ line.split() for line in f ]
    input_data = [ format_line( line, site_charge, offset ) for line in input_data if ( float(line[2]) > x_min and float(line[2]) < x_max ) ]
    input_data = sorted( input_data, key=lambda x: x[2] )
    input_data = cluster_similar_sites( input_data )
#    for line in input_data:
#        print(line[0],line[2], flush=True)

    return input_data

def cluster_similar_sites( input_data ):
    """Clusters data points the input data to be projected onto a site. Clustering criterion is set as 0.01 nm.
        
        Args:
        input_data (str): The input file.
        
        Returns:
        str: The input file, formatted with x coordinate values with correct clustering coordinates.
        
        """
    coordinates = np.array([[ion[2]] for ion in input_data]) # format data into correct format
    # Use sklearn's clustering functionality to assign site clusters
    cluster = AgglomerativeClustering(n_clusters=None, affinity='euclidean', compute_full_tree = True,linkage='ward',distance_threshold = 1e-10)
    cluster_data = cluster.fit_predict(coordinates)
    # Match clusters to sites
    sites = [ [] for i in range(0,np.max(cluster_data)+1)]
    for cluster in range(len(cluster_data)):
        sites[cluster_data[cluster-1]].append(input_data[cluster-1])
    sites = sorted(sites,key = lambda x: x[0][2])
    # Prepare data in the correct format.
    sites_format = []
    for cluster in sites:
        site_loc = np.mean(np.array(cluster)[:,2].astype('float64'))
        for ion in cluster:
            ion[2] = site_loc
        sites_format += cluster
    return sites_format

def calculate_grid_offsets( filename, x_min, x_max, system ):
    """Reads in the input data calculates the distance to the next site outside of the defined calculation region. Allows calculation of the delta_x and volume values for the endmost grid points.

    Args:
        filename(string): Filename for importing the input data.
        x_min(float): Minimum x coordinate value defining the calculation region.
        x_max(float): Maximum x coordinate value defining the calculation region. 
        system(str): 'single' for single grain boundary systems, 'double' for systems where the input data has been mirrored to give two grain boundaries. 

    Returns:
        list, list: distance between the midpoint of the endmost sites and the midpoint of the next site outside of the calculation region for the first and last sites respectively. Distance between the endmost sites and the next site outside of the calculation region for the first and last sites respectively.

    """
    with open( filename, 'r' ) as f:
        input_data = [ line.split() for line in f ]
        input_data = [ format_line( line, 0.0, 0.0 ) for line in input_data ]
        input_data = sorted( input_data, key=lambda x: x[2] )
        input_data = cluster_similar_sites( input_data )
        x_coords = np.unique( np.array( [ float(l[2]) for l in input_data ] ) )
        min_index = bisect_left( x_coords, x_min )
        max_index = bisect_right( x_coords, x_max )
        min_offset = ((x_coords[min_index+1]-x_coords[min_index])/2)+((x_coords[min_index]-x_coords[min_index-1])/2)
        max_offset = ((x_coords[max_index+1]-x_coords[max_index])/2)+((x_coords[max_index]-x_coords[max_index-1])/2)
        if system == 'single':
            limits = [ min_offset, max_offset ]
            limits_for_laplacian = [ (x_coords[min_index]-x_coords[min_index-1]), (x_coords[max_index]-x_coords[max_index-1]) ]
        if system == 'double':
            limits = [min_offset, min_offset]
            limits_for_laplacian = [ (x_coords[min_index]-x_coords[min_index-1]), (x_coords[min_index]-x_coords[min_index-1]) ]
    return limits, limits_for_laplacian

