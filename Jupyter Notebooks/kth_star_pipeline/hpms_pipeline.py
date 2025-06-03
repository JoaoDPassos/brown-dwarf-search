import lsdb
import pandas as pd
import numpy as np

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
def n_neighbors_filter(df, min_neighbors, id_col):

    neighbors = df.groupby(id_col)['_dist_arcsec'].count()
    neighbors -= 1 #Double counting adjustment
    neighbors.name = 'neighbors'
    df = df.join(neighbors, on=id_col)
    xmatch_large_groups = df.query(f'neighbors >= {min_neighbors}')

    return xmatch_large_groups

'''Evaluates the magnitude of the cross product of two vectors. When one vector is a vector from a line in 2D space and the other
is a vector pointing from a point on that line to some arbitrary point in 2D space, Q, this will return the distance of point Q 
from that line. This is relevant to finding the deviation of objects from projections.

Args:
    - line_vector: line vector of a projection
    - PQ: vector starting at point P and pointing to point Q

Returns: Magnitude of the cross product of two vectors, or distance of point Q from line with line_vector.
'''
def distance_to_line(line_vector, PQ):
    #Add zero for proper cross product
    PQ_3D = np.append(PQ, 0)
    line_vector_3D = np.append(line_vector, 0)
    return np.linalg.norm(np.abs(np.cross(PQ_3D, line_vector_3D)) / np.linalg.norm(line_vector_3D))

'''The kth star algorithm. Given a group of objects, this algorithm will project lines between two objects and calculate the 
kth smallest deviation of other stars within the group from this line. For example, if k = 4 and the kth minimum distance is zero,
then we know that 4 stars are in perfect alignment, with may be one high PM star.

Args:
    - group: groupby object which holds all the data of one group of objects
    - k: number of objects we seek to filter for alignment

Returns: group passed in as an argument with an additional column containing the kth smallest deviation from a line projection.
'''
# Consider making the RA and DEC columns arguments into this func
def kth_star_min_distance(group, k):
    origin_ra, origin_dec = group['RA_1'], group['DEC_1']
    ra2, dec2 = group[["RA_2", "DEC_2"]].to_numpy().T
    
    x_vals = (ra2 - origin_ra) * np.cos(np.radians(origin_dec)) * 3600
    y_vals = (dec2 - origin_dec) * 3600
    delta_coords = np.vstack((x_vals, y_vals)).T 

    kth_distances = [None] * len(delta_coords)

    for i in range(len(delta_coords)):
        proj_vector = delta_coords[i]
        distances = []

        if np.all(proj_vector == 0): continue #ensures vector points to a point that is not the origin

        for j in range(len(delta_coords)):
            if np.all(delta_coords[j] == 0) or (j == i): continue
            
            distances.append(distance_to_line(proj_vector, delta_coords[j]))

        
        kth_distances[i] = sorted(distances)[k]

    group['kth_min_proj_error'] = kth_distances

    return group['kth_min_proj_error']

'''map_partitions() compatible function which runs the kth star algorithm on a catalog already crossmatched
and filtered for groups.

Args:
    - df: passed by map_partitions(), dataframe which is modified by function
    - k: number of stars which we seek to be in aligment

Returns: Catalog with "kth_min_proj_error" column, which is the kth smallest deviation of a star with alignment
with other stars.
'''
def apply_kth_star(df, k, id_col):
    
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
def execute_pipeline(catalog, query_string,
                     xmatch_max_neighbors, max_neighbor_dist, min_neighbors,
                     k, max_obj_deviation, id_col):
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