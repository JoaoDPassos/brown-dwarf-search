import lsdb
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from itertools import combinations

class Star:
    def __init__(self, ra, dec, mags):
        self.ra = ra
        self.dec = dec
        self.mags = mags
        self.deviation = np.nan
        self.proj_ra = np.nan
        self.proj_dec = np.nan
    
    def __lt__(self, other):
        # If self is NaN and other is not, self is greater (goes later)
        if np.isnan(self.deviation) and not np.isnan(other.deviation):
            return False
        # If other is NaN and self is not, self is less (goes earlier)
        if not np.isnan(self.deviation) and np.isnan(other.deviation):
            return True
        # If both are NaN, treat as equal (but keep original order)
        if np.isnan(self.deviation) and np.isnan(other.deviation):
            return False
        # Normal comparison
        return self.deviation < other.deviation

    def __eq__(self, other):
        return ((self.ra == other.ra) and 
                (self.dec == other.dec) and
                (self.mags == other.dec) and
                (self.deviation == other.dec) and
                (self.proj_ra == other.dec) and
                (self.proj_dec == other.dec))


    def __repr__(self):
        return f'RA: {self.ra}, DEC: {self.dec}, mags: {self.mags}, deviation: {self.deviation}'

    def __str__(self):
        return repr(self)


def n_neighbors_filter(df, min_neighbors, id_col):
    
    '''map_partitions() compatible function which filters catalog for objects groups of a minimum size. Objects in the same group
    are considered "neighbors". It is assumed that this catalog has already been self-crossmatched.
    
    Args:
        - df: passed by map_partitions(), dataframe which is filtered by function
        - min_neighbors: minimum groupsize to filter for, any groups smaller than this will be filtered out
        - id_col: Column with unique object identifier. Is generally different in every catalog and is what we use to distinguish 
          between objects.
    
    Returns: Dataframe filtered for group sizes greater than or equal to min_neighbors.
    
    Note: Because we are crossmatching without margin cache, we will be losing some potential matches at the border of healpix pixels.
    '''
    
    neighbors = df.groupby(id_col)['_dist_arcsec'].count()
    neighbors -= 1 #Double counting adjustment
    neighbors.name = 'neighbors'
    df = df.join(neighbors, on=id_col)
    xmatch_large_groups = df.query(f'neighbors >= {min_neighbors}')

    return xmatch_large_groups


def distance_to_line(line_vector, PQ):
    
    '''Evaluates the magnitude of the cross product of two vectors. When one vector is a vector from a line in 2D space and the other
    is a vector pointing from a point on that line to some arbitrary point in 2D space, Q, this will return the distance of point Q 
    from that line. This is relevant to finding the deviation of objects from projections.
    
    Args:
        - line_vector: line vector of a projection
        - PQ: vector starting at point P and pointing to point Q
    
    Returns: Magnitude of the cross product of two vectors, or distance of point Q from line with line_vector.
    '''
    
    #Add zero for proper cross product
    PQ_3D = np.append(PQ, 0)
    line_vector_3D = np.append(line_vector, 0)
    return np.linalg.norm(np.abs(np.cross(PQ_3D, line_vector_3D)) / np.linalg.norm(line_vector_3D))


# Time complexity can be improved to O(nlogn) using convex hull: https://www.geeksforgeeks.org/maximum-distance-between-two-points-in-coordinate-plane-using-rotating-calipers-method/
def find_max_obj_distance(neighbors):
    
    ''' Takes a list of neighbors in alignment and calculates the maximum distance between two neighbors in the list.

    Args:
        - neighbors: list of Star class objects. 

    Returns: maximum distance as np.float64(?)
    '''
    
    # Does not account for edgecases of dec and ra (?)

    max_distance = 0
    for star_1, star_2 in combinations(neighbors, 2):
        coord_1 = SkyCoord(ra=star_1.ra * u.deg, dec=star_1.dec * u.deg, frame='icrs')
        coord_2 = SkyCoord(ra=star_2.ra * u.deg, dec=star_2.dec * u.deg, frame='icrs')
        distance = coord_1.separation(coord_2)
        if distance > max_distance:
            max_distance = distance

    max_distance_arcsec = max_distance.arcsecond  
    return max_distance_arcsec
    
def find_max_mag_diff(neighbors, passbands):
    ''' Takes a list of neighbors in alignment and calculates the maximum difference between magnitudes across passbands.
    We compare magnitudes between stars in the same passband.

    Args:
        - neighbors: list of Star class objects. 
        - passbands: array of strings decribing the passbands (i.e. ['g','r'])

    Returns: maximum mag difference as np.float64(?)
    '''
    
    return np.nan

def get_aligned_neighbors(sorted_neighbors, k):
    
    ''' Takes an array of Star class objects which is already sorted, and returns the subset of 
    neighbors which are in alignment.

    Args:
        - sorted_neighbors: array of Star class objects sorted using the class' less than function.
        - k: number of stars we seek to be in alignment with the projection.

    Returns: An array of Star class objects, each of which is in alignment.
    
    Explanation: The two stars forming the projection are never assigned a deviation value; hence,
    their .deviation attribute is nan. By the class definition, nan is sorted to the end of the list,
    meaning that the final two objects in sorted_neighbors are always the two forming the projection.
    The remaining neighbors in alignment are from [0, k].
    '''

    proj_stars = sorted_neighbors[-2:]
    remaining_neighbors = sorted_neighbors[:k+1]
    return proj_stars + remaining_neighbors

def reset_deviations(neighbors):

    ''' Takes in a list of star class objects and sets every neighbor's deviation value to nan.

    Args:
        - neighbors: Array of star class objects.

    Returns: None.
    '''
    
    for neighbor in neighbors: 
        neighbor.deviation = np.nan
    return None

def get_proj_coords(origin, neighbors):
    
    ''' Creates an array of 2D coordinates of the stars in neighbors, using the coordinates from origin as the origin
    of the system. Also stores these coordinates under .proj_ra and .proj_dec for each star in neighbors. We assume that
    the origin is one of the stars in neighbors so that the array lengths are in agreement in kth_star_min_distance.

    Args:
        - origin: Star class object, who's coordinates serve as the origin of our coordinate system.
        - neighbors: List of star class objects which we are obtaining the coordinates for in our new coordinate system.

    Returns: n x 2 array of coordinates for stars in coordinates with origin as the origin.
    '''

    neighbor_ra, neighbor_dec = np.array([star.ra for star in neighbors]), np.array([star.dec for star in neighbors])
    x_vals = (neighbor_ra - origin.ra) * np.cos(np.radians(origin.dec)) * 3600
    y_vals = (neighbor_dec - origin.dec) * 3600

    res = np.vstack((x_vals, y_vals)).T 
    
    for i in range(len(res)):
        neighbors[i].proj_ra, neighbors[i].proj_dec = res[i][0], res[i][1]

    return res

def kth_star_min_distance(group, k, cols_to_keep, mag_cols, max_obj_deviation):
    
    mag_cols_1 = [f'{col}_1' for col in mag_cols]
    origin_star = Star(group['RA_1'].iloc[0], group['DEC_1'].iloc[0], group[mag_cols_1].iloc[0])
    mag_cols_2 = [f'{col}_2' for col in mag_cols]
    neighbors = [Star(row['RA_2'], row['DEC_2'], row[mag_cols_2]) for _,row in group.iterrows()]

    proj_coords = get_proj_coords(origin_star, neighbors)
    print(f"Showing proj_coords={proj_coords}")

    kth_distances = [None] * len(proj_coords)
    max_distances = [None] * len(proj_coords)
    max_mag_diffs = [None] * len(proj_coords)

    for i in range(len(neighbors)):
        proj_vector = proj_coords[i]
        distances = []

        if np.all(proj_vector == 0): continue #ensures vector points to a point that is not the origin

        for j in range(len(neighbors)):
            if np.all(proj_coords[j] == 0) or (j == i): continue
            print(f"proj_vector = {proj_vector}, proj_coords[{j}] = {proj_coords[j]}, distance_to_line(proj_vector, proj_coords[{j}]), {distance_to_line(proj_vector, proj_coords[j])}")
            neighbors[j].deviation = distance_to_line(proj_vector, proj_coords[j])

        sorted_neighbors = sorted(neighbors)
        
        # If the k value is larger than the number of neighbors or the deviation is too great, 
        # skip postfiltering and declare the projection invalid.
        if ((k > len(neighbors) - 2) or (sorted_neighbors[k].deviation > max_obj_deviation)): 
            kth_distances[i] = np.nan
            max_distances[i] = np.nan
            max_mag_diffs[i] = np.nan
            reset_deviations(neighbors)
            continue
        
        # We have k + 2 objects in alignment, let's save maximum distance between objs and max mag diff!
        for i in range(len(sorted_neighbors)):
            print(f"{i}th sorted_neighbors value = {sorted_neighbors[i].deviation}")

        print('done')
        kth_distances[i] = sorted_neighbors[k].deviation
        aligned_neighbors = get_aligned_neighbors(neighbors, k)
        max_distances[i] = find_max_obj_distance(aligned_neighbors)
        max_mag_diffs[i] = find_max_mag_diff(aligned_neighbors, mag_cols_1)

        reset_deviations(neighbors)

        # Reset deviation value for next projection
    
    group['kth_min_proj_error'] = kth_distances
    group['max_obj_distance'] = max_distances
    group['max_mag_diff'] = max_mag_diffs

    return group[cols_to_keep + ['kth_min_proj_error'] + ['max_obj_distance']]

def apply_kth_star(df, k, id_col):

    '''map_partitions() compatible function which runs the kth star algorithm on a catalog already crossmatched
    and filtered for groups.
    
    Args:
        - df: passed by map_partitions(), dataframe which is modified by function
        - k: number of stars which we seek to be in aligment
    
    Returns: Catalog with "kth_min_proj_error" column, which is the kth smallest deviation of a star with alignment
    with other stars.
    '''
    
    if df.empty:
        df['kth_min_proj_error'] = pd.Series(dtype=float)
    else:
        df.reset_index(inplace=True)
        df['kth_min_proj_error'] = (
            df.groupby(id_col)
              .apply(kth_star_min_distance, k=k-1)
              .reset_index(drop=True, level=0)
        )
        df.set_index('_healpix_29')

    cols_to_keep = [col for col in df.columns if col.endswith('_2')] + ['kth_min_proj_error']
    return df[cols_to_keep]

def execute_pipeline(catalog, query_string,
                     xmatch_max_neighbors, max_neighbor_dist, min_neighbors,
                     k, max_obj_deviation, id_col):

    '''Executes the high PM star data pipeline.

    Args:
        - catalog: Catalog in HATS format we are seeking high PM stars in
        - query_string: string for initial filter of dataframe, in our case this is used for filtering for stars
        - xmatch_max_neighbors: when performing crossmatching, this is the maximum size of a group which the crossmatch will return
        - max_neighbor_dist: maximum distance between neighbors
        - k: The number of stars which we seek to be in alignment
        - max_obj_deviation: The maxmimum we want an object to deviate from alignment.
        - id_col: Column with unique object identifier. Is generally different in every catalog and is what we use to distinguish 
          between objects.
    
    Returns: Catalog filtered for high PM stars using the kth star algorithm.
    
    Note: Postfiltering is not yet implemented.
    '''
    
    # Filter for stars and objects with reasonable measurement errors
    star_filtered = catalog.query(query_string)

    # Crossmatch catalog with itself and filter for group sizes >= min_ neighbors
    xmatch = star_filtered.crossmatch(
                star_filtered,
                n_neighbors=xmatch_max_neighbors,
                radius_arcsec=max_neighbor_dist,
                suffixes=['_1', '_2']
    )
    neighbors_filtered = xmatch.map_partitions(n_neighbors_filter, min_neighbors, id_col)
    
    # Add column for kth minimum distance
    kth_star = neighbors_filtered.map_partitions(apply_kth_star, k, id_col)
    kth_star_filtered = kth_star.query(f'kth_min_proj_error < {max_obj_deviation}')
    
    return kth_star_filtered