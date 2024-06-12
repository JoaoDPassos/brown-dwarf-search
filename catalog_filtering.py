#!/usr/bin/env python3
# Calls query method to filter data frame given certain data requirements

def bandFilter(df, band, classStarMin = None, spreadMax = None, flag = False, maxMagError = None):
    queryString = createQueryString(band, classStarMin, spreadMax, flag, maxMagError)
    if queryString != '':
        queryString += f' and WAVG_MAG_PSF_{band} > 0 and WAVG_MAG_PSF_{band} < 90'
    else:
        queryString += f'WAVG_MAG_PSF_{band} > 0 and WAVG_MAG_PSF_{band} < 90'
    return df.query(queryString)

#Generates the string for .query method
    
def createQueryString(band, classStarMin, spreadMax, flag, maxMagError):
    queryString = ''
    if classStarMin != None:
        if queryString != '': queryString += ' and '
        queryString += f'CLASS_STAR_{band} > {classStarMin}'
    if spreadMax != None:
        if queryString != '': queryString += ' and '
        queryString += f'SPREAD_MODEL_{band} < {spreadMax}'
    if flag:
        if queryString != '': queryString += ' and '
        queryString += f'FLAGS_{band} < 4'
    if maxMagError != None:
        if queryString != '': queryString += ' and '
        queryString += f'WAVG_MAGERR_PSF_{band} < {maxMagError}'
    return queryString

#Filters through multiple bands, makes sure objects satisfy filters for ALL bands (and)

def multiBandFilterStrict(df, bandList, classStarMin = None, spreadMax = None, flag = False, maxMagError = None):
    for band in bandList:
        df = bandFilter(df, band, classStarMin, spreadMax, flag, maxMagError)
    return df

#Filters through multiple bands, makes sure objects satisfy filters for AT LEAST 1 band (or)

def multiBandFilterLenient(df, bandList, classStarMin = None, spreadMax = None, flag = False, maxMagError = None):
    queryString = ''
    strictFilters = ''
    for band in bandList:
        queryString += f'({createQueryString(band, classStarMin, spreadMax, flag, maxMagError)}) or '
        strictFilters += f'WAVG_MAG_PSF_{band} > 0 and WAVG_MAG_PSF_{band} < 90 and ' #False values must always be filtered
    queryString = queryString[:-3]
    finalString = strictFilters + f'({queryString})' if queryString != '' else strictFilters
    return df.query(finalString)