import pandas as pd
import numpy as np
import random
from astropy.coordinates import SkyCoord
from astropy import units as u
from itertools import combinations
from dataclasses import dataclass
from IPython.display import display
from collections import namedtuple
import networkx as nx


@dataclass (eq=True)
class Star:
    def __init__(self, ra, dec, mags, mag_errs):
        self.ra = ra
        self.dec = dec
        self.mags = mags
        self.mag_errs = mag_errs
        self.deviation = np.inf

    def __repr__(self):
        return f'RA: {self.ra}, DEC: {self.dec}, mags: {self.mags}, mag_errs: {self.mag_errs}, deviation: {self.deviation}'

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
    
    neighbors = df.groupby(id_col)['_dist_arcsec'].transform("count") - 1 # -1 for double counting adjustment
    return df.loc[neighbors >= min_neighbors]



def distance_to_line(line_vector, PQs):
    
    '''Evaluates the magnitude of the cross product of two vectors in parallel. When one vector is a vector from a line in 2D space and the other
    is a vector pointing from a point on that line to some arbitrary point in 2D space, Q, this will return the distance of point Q 
    from that line. This is relevant to finding the deviation of objects from projections.
    
    Args:
        - line_vector: line vector of a projection
        - PQ: nx2 array of vectors starting at point P and pointing to point Q
    
    Returns: 1xn array of magnitudes of the cross products of two vectors, or distances of multiple poins Q from line with line_vector.
    '''
    
    # Add zero z-component to make vectors 3D for cross product
    PQ_3D = np.hstack((PQs, np.zeros((PQs.shape[0], 1))))
    line_vector_3D = np.append(line_vector, 0)

    # Compute cross product (n x 3), then get z-component and normalize
    cross_products = np.cross(PQ_3D, line_vector_3D)
    magnitudes = np.linalg.norm(cross_products, axis=1)
    del PQ_3D
    
    # Divide by norm of the line vector
    return magnitudes / np.linalg.norm(line_vector_3D)


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
    
def find_max_mag_diff(neighbors, mag_cols):
    
    ''' Takes a list of neighbors in alignment and calculates the maximum difference between magnitudes across passbands.
    We compare magnitudes between stars in the same passband.

    Args:
        - neighbors: list of Star class objects. 
        - passbands: array of strings decribing the passbands (i.e. ['G','R'])

    Returns: maximum mag difference as np.float64(?)
    '''
    
    all_mag_lists = [neighbor.mags for neighbor in neighbors]
    max_diff = np.nan

    for i in range(len(mag_cols)):
        mags = [mag_list[i] for mag_list in all_mag_lists]

        for idx_1, idx_2 in combinations(range(len(mags)), 2):
            mag_1, mag_2 = mags[idx_1], mags[idx_2]
            if (mag_1 == -99.0) or (mag_2 == -99.0): continue 

            mag_1_err, mag_2_err = neighbors[idx_1].mag_errs[i], neighbors[idx_2].mag_errs[i]
            err_norm = np.hypot(mag_1_err, mag_2_err)
                
            curr_diff = abs(mag_1 - mag_2) / max(err_norm, 0.01)
            if ((curr_diff > max_diff) or (np.isnan(max_diff))):
                max_diff = curr_diff
    del all_mag_lists
    return max_diff

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

    ''' Takes in a list of star class objects and sets every neighbor's deviation value to infinity.

    Args:
        - neighbors: Array of star class objects.

    Returns: None.
    '''
    
    for neighbor in neighbors: 
        neighbor.deviation = np.inf
    return 

def get_proj_coords(origin, neighbors):
    
    ''' Creates an array of 2D coordinates of the stars in neighbors, using the coordinates from origin as the origin
    of the system. We assume that the origin is one of the stars in neighbors so that the array lengths are in agreement in
    kth_star_min_distance.

    Args:
        - origin: Star class object, who's coordinates serve as the origin of our coordinate system.
        - neighbors: List of star class objects which we are obtaining the coordinates for in our new coordinate system.

    Returns: n x 2 array of coordinates for stars in coordinates with origin as the origin.
    '''

    neighbor_ra, neighbor_dec = np.array([star.ra for star in neighbors]), np.array([star.dec for star in neighbors])
    x_vals = (neighbor_ra - origin.ra) * np.cos(np.radians(origin.dec)) * 3600
    y_vals = (neighbor_dec - origin.dec) * 3600
    del neighbor_ra
    del neighbor_dec

    return np.vstack((x_vals, y_vals)).T 


# def kth_star_min_distance(group, k, mag_cols, max_obj_deviation, id_col, debug_mode=False):
#     mag_cols_2 = [f"{col}_2" for col in mag_cols]
#     mag_col_errs = [f"{col.replace('WAVG_MAG_PSF', 'WAVG_MAGERR_PSF')}_2" for col in mag_cols]
#     neighbors = [
#         Star(row.RA_2, row.DEC_2, 
#              [getattr(row, col) for col in mag_cols_2], 
#              [getattr(row, col) for col in mag_col_errs])
#         for row in group.itertuples()
#     ]
    
#     # Group order is sorted by distances; hence, the origin star is always index 0 (col_1[0] = col_2[0])
#     origin_star = neighbors[0]

#     proj_coords = get_proj_coords(origin_star, neighbors)

#     kth_distances = [None] * len(proj_coords)
#     max_distances = [None] * len(proj_coords)
#     max_mag_diffs = [None] * len(proj_coords)
#     # debug_col = [None] * len(proj_coords)

#     for i in range(len(neighbors)):
#         proj_vector = proj_coords[i]

#         if np.all(proj_vector == 0): continue #ensures vector points to a point that is not the origin

#         # Compute all distances in parallel
#         distances = distance_to_line(proj_vector, proj_coords)
#         # print("Distances returned from dist_to_line: ", distances)

#         for j in range(len(distances)):
#             if np.all(proj_coords[j] == 0) or (j == i): continue
#             neighbors[j].deviation = distances[j]
#             # print(f"neighbors[{j}].deviation = {distances[j]}")

#         sorted_neighbors = sorted(neighbors, key=lambda item: item.deviation)
        
#         # If the k value is larger than the number of neighbors or the deviation is too great, 
#         # skip postfiltering and declare the projection invalid.
#         if ((k > len(neighbors) - 2) or (sorted_neighbors[k].deviation > max_obj_deviation)): 
#             kth_distances[i] = np.inf
#             max_distances[i] = np.inf
#             max_mag_diffs[i] = np.inf
#             reset_deviations(neighbors)
#             del sorted_neighbors
#             continue
        
#         # We have k + 2 objects in alignment, let's save maximum distance between objs and max mag diff!
#         kth_distances[i] = sorted_neighbors[k].deviation
#         aligned_neighbors = get_aligned_neighbors(sorted_neighbors, k)
#         # debug_col[i] = aligned_neighbors
#         max_distances[i] = find_max_obj_distance(aligned_neighbors)
#         max_mag_diffs[i] = find_max_mag_diff(aligned_neighbors, mag_cols_2)

#         # Reset deviation value for next projection
#         reset_deviations(neighbors)
    
#     group['kth_min_deviation'] = kth_distances
#     group['max_obj_distance'] = max_distances
#     group['max_mag_diff'] = max_mag_diffs
#     # group['aligned_neighbors'] = debug_col

    
#     return_cols = ['kth_min_deviation', 'max_obj_distance', 'max_mag_diff']
#     # if debug_mode: return_cols.append('aligned_neighbors')
         
#     return group[return_cols]

def kth_star_min_distance(group, k, mag_cols, max_obj_deviation, id_col, debug_mode=False):
    mag_cols_2 = [f"{col}_2" for col in mag_cols]
    mag_col_errs = [f"{col.replace('WAVG_MAG_PSF', 'WAVG_MAGERR_PSF')}_2" for col in mag_cols]

    neighbors = [
        Star(row.RA_2, row.DEC_2, 
             [getattr(row, col) for col in mag_cols_2], 
             [getattr(row, col) for col in mag_col_errs])
        for row in group.itertuples()
    ]
    # Group order is sorted by distances; hence, the origin star is always index 0 (col_1[0] = col_2[0])
    origin_star = neighbors[0]
    proj_coords = get_proj_coords(origin_star, neighbors)

    # Track best
    best_dev = np.inf
    best_max_dist = None
    best_max_mag_diff = None
    best_index = None  # keep track of which neighbor gave the best

    for i, proj_vector in enumerate(proj_coords):
        if np.all(proj_vector == 0): #ensures vector points to a point that is not the origin
            continue
        # Compute all distances in parallel
        distances = distance_to_line(proj_vector, proj_coords)
        for j in range(len(distances)):
            if np.all(proj_coords[j] == 0) or (j == i):
                continue
            neighbors[j].deviation = distances[j]

        sorted_neighbors = sorted(neighbors, key=lambda item: item.deviation)
        if k > len(neighbors) - 2:
            reset_deviations(neighbors)
            continue

        kth_dev = sorted_neighbors[k].deviation
        if kth_dev <= max_obj_deviation and kth_dev < best_dev:
            aligned_neighbors = get_aligned_neighbors(sorted_neighbors, k)
            best_dev = kth_dev
            best_max_dist = find_max_obj_distance(aligned_neighbors)
            best_max_mag_diff = find_max_mag_diff(aligned_neighbors, mag_cols_2)
            best_index = i

        reset_deviations(neighbors)

    if best_index is None:
        # no valid candidate: return empty DataFrame
        return group.iloc[0:0]

    # pick one row to represent the group
    best_row = group.iloc[[best_index]].copy()
    best_row['kth_min_deviation'] = best_dev
    best_row['max_obj_distance'] = best_max_dist
    best_row['max_mag_diff'] = best_max_mag_diff

    return best_row


def apply_kth_star(df, k, id_col, mag_cols, max_obj_deviation, debug_mode=False):

    '''map_partitions() compatible function which runs the kth star algorithm on a catalog already crossmatched
    and filtered for groups.
    
    Args:
        - df: passed by map_partitions(), dataframe which is modified by function
        - k: number of stars which we seek to be in aligment
    
    Returns: Catalog with "kth_min_deviation" column, which is the kth smallest deviation of a star with alignment
    with other stars.
    '''
    original_index = df.index
    original_index_name = df.index.name
    df.reset_index(inplace=True)
    result = df
    if df.empty:
        df['kth_min_deviation'] = pd.Series(dtype=float)
        df['max_obj_distance'] = pd.Series(dtype=float)
        df['max_mag_diff'] = pd.Series(dtype=float)
        df['aligned_neighbors'] = pd.Series(dtype=float)
    else:
        result = (
            df.groupby(id_col)
              .apply(kth_star_min_distance, k=k-1, mag_cols=mag_cols, 
                     max_obj_deviation=max_obj_deviation, id_col=id_col, debug_mode=debug_mode)
              .reset_index(drop=True, level=0)
        )
        # print("df after join K*:")
        # display(df)

    df = merge_groups_and_reduce(result, id_col, max_obj_deviation)
    df.index = df[original_index_name]
    df = df.drop(original_index_name,axis=1)

    hpms_col_names = ['kth_min_deviation', 'max_obj_distance', 'max_mag_diff']
    if debug_mode: hpms_col_names.append('aligned_neighbors') 
    
    cols_to_keep = [id_col] + [col for col in df.columns if col.endswith('_2')] + hpms_col_names
    return df[cols_to_keep]

def num_missed_detections(df, passband_cols, missed_val):

    '''map_patitions() compatible function which counts the number of missed detections across all passbands, storing this in
    a new column.

    Args:
        - df: Passed by map_partitions method, and is the df from which we obtain the missed detection statistics.
        - passband_cols: List of passband columns names (strings) of which we will check for missed detections.
        - missed_val: Some catalogs, like DES-DR2, have unique values recorded when there is a missed detection. 
          Thus, we pass this value into this function to ensure that we can find the missed detections in a general catalog.
    Returns: df with n_missed_detections column describing the total number of missed detections.
    '''

    n_missed_detections=[]

    for idx, row in df.iterrows():
        passbands = row[passband_cols]
        n_missed_detections.append((passbands == missed_val).values.sum())
        
    df['n_missed_detections']=n_missed_detections
    return df['n_missed_detections']
    

def execute_pipeline(catalog, query_string,
                     xmatch_max_neighbors, max_neighbor_dist, min_neighbors,
                     k, max_obj_deviation, id_col, mag_cols, debug_mode=False):

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

    # Crossmatch catalog with itself and filter for group sizes >= min_neighbors
    xmatch = star_filtered.crossmatch(
                star_filtered,
                n_neighbors=xmatch_max_neighbors,
                radius_arcsec=max_neighbor_dist,
                suffixes=['_1', '_2']
    )
    # if debug_mode: 
    #     print(f"Length of self crossmatch: {len(xmatch.compute())}")
        
    neighbors_filtered = xmatch.map_partitions(n_neighbors_filter, min_neighbors, id_col)
    del xmatch

    if debug_mode: print(f"Length after neighbors filter: {len(neighbors_filtered.compute())}")
    
    # Add column for kth minimum distance and postfiltering parameters
    kth_star = neighbors_filtered.map_partitions(apply_kth_star, k, id_col, mag_cols, max_obj_deviation, debug_mode)
    del neighbors_filtered
    if debug_mode: print(f"Length after kth star filter: {len(kth_star.compute())}")
    kth_star_filtered = kth_star.query(f'kth_min_deviation < {max_obj_deviation}')
    del kth_star
    
    return kth_star_filtered


def sample_mag_diffs(catalog, band_col_name, arr_size):
    n_partitions=catalog.npartitions
    diffs = []

    for _ in range(arr_size):
        id_1 = random.randint(0, n_partitions - 1)
        obj_1 = catalog.sample(id_1, n=1)
            
        id_2 = random.randint(0, n_partitions - 1)
        obj_2 = catalog.sample(id_2, n=1)
        
        if (obj_1[band_col_name].item() == -99.0) or (obj_2[band_col_name].item() == -99.0):
            diffs.append(-1)
        else:
            diffs.append(abs(obj_1[band_col_name].item() - obj_2[band_col_name].item()))
    
    return np.array(diffs)

def merge_groups_and_reduce(df, id_1, max_obj_deviation):
    # Step 1: Build a graph where stars are nodes, connected if they are in the same group
    G = nx.Graph()
    id_2 = id_1[:-1] + '2'

    for group_id, group_df in df.groupby(id_1):
        #TODO: query here for kth_min_deviation which we want before adding them to the group
        filtered_group = group_df[
            (group_df['kth_min_deviation'].notna()) &
            (group_df['kth_min_deviation'] < max_obj_deviation)
        ]
        stars_2 = filtered_group[id_2].unique()
        star_1 = filtered_group[id_1].unique()
        # Connect all pairs of stars in this group
        for i in range(len(stars_2)):
            G.add_edge(star_1[0], stars_2[i])
            for j in range(i + 1, len(stars_2)):
                G.add_edge(stars_2[i], stars_2[j])
                
    # Step 2: Find connected components (merged groups)
    connected_components = list(nx.connected_components(G))

    # Step 3: Map each star_id to a merged_group_id (e.g., by index)
    star_to_component = {}
    for i, component in enumerate(connected_components):
        for star_id in component:
            star_to_component[star_id] = i

    # Step 4: Assign merged group id to original df
    df['merged_group_id'] = df[id_2].map(star_to_component)

    # Step 5: Reduce each merged group to the star with the smallest kth_min_deviation
    reduced_df = df.loc[df.groupby('merged_group_id')['kth_min_deviation'].idxmin()]

    return reduced_df
