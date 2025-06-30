#!/usr/bin/env python3
# Calls query method to filter data frame given certain data requirements
import numpy as np
from dask import delayed, compute

def createQueryString(band, classStar, spreadModel, magError, flag, mag):
    queries = []

    if classStar != None:
        queries.append(f'CLASS_STAR_{band} > {classStar}')
    if spreadModel != None:
        queries.append(f'SPREAD_MODEL_{band} < {spreadModel}')
    if flag:
        queries.append(f'FLAGS_{band} < 4')
    if magError != None:
        queries.append(f'WAVG_MAGERR_PSF_{band} < {magError}')
    if mag != None:
        queries.append(f'(WAVG_MAG_PSF_{band} > {mag} or WAVG_MAG_PSF_{band} < 0)')
    
    query_string = ' and '.join(queries)
    return query_string

def createFilterTuples(band, classStar, spreadModel, magError, flag, mag):
    queries = []

    if classStar != None:
        queries.append((f'CLASS_STAR_{band}','>',classStar))
    if spreadModel != None:
        queries.append((f'SPREAD_MODEL_{band}','<', spreadModel))
    if flag:
        queries.append((f'FLAGS_{band}','<', 4))
    if magError != None:
        queries.append((f'WAVG_MAGERR_PSF_{band}', '<', magError))
    if mag != None:
        queries.extend([(f'WAVG_MAG_PSF_{band}','>', mag), (f'WAVG_MAG_PSF_{band}', '<', 90)])

    return queries


#Filters through multiple bands, makes sure objects satisfy filters for ALL bands (and)

def bandFilterStrict(bandList, classStar=None, spreadModel=None, magError=None, flag=False, invalidMags=False):
    query_parts = []

    for band in bandList:
        user_params = createQueryString(band, classStar, spreadModel, magError, flag, invalidMags)
        if user_params != '':
            query_parts.append(f'{user_params}')

    return ' and '.join(query_parts)


#Filters through multiple bands, makes sure objects satisfy filters for AT LEAST 1 band (or)

def bandFilterLenient(bandList, classStar=None, spreadModel=None, magError=None, flag=False, mag=None):
    query_parts = []

    for band in bandList:
        user_params = createQueryString(band, classStar, spreadModel, magError, flag, mag)
        if user_params != '':
            query_parts.append(f'{user_params}')

    return ' or '.join(query_parts)

#checks if the dataframe contains our known PM star, also can add more as an argument

def contains_PM(df, PM_set=[10370986892068913152, 10370986891217469440, 10370986798997307392, 10370986798804369408]):
    index_array = np.array(df.index, dtype=object)
    l = np.isin(index_array, PM_set)
    print(l)
    if np.any(np.isin(index_array, PM_set)):
        return True
    return False

#convert PM set to uint64 not object

# assert contains_PM([10370986892068913152, 10370986891217469440, 10370986798997307392, 10370986798804369408, 0, 1, 22, 86]) == True
# assert contains_PM([10370986892068913152, 10370986798997307392, 10370986798804369408, 80, 390902390482039, 0, -2]) == True
# assert contains_PM([0,1,2,3,5]) == False
# assert contains_PM([]) == False

def create_filter_list(bandList, classStar=None, spreadModel=None, magError=None, flag=False, invalidMags=False):
    final_filters = []

    for band in bandList:
        user_params = createFilterTuples(band, classStar, spreadModel, magError, flag, invalidMags)
        final_filters.append(user_params)

    return final_filters
