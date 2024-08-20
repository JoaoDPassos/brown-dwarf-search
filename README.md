# Purpose

This repository is a collection of various Jupyter Notebooks and python files which were created with the intention of exploring new ways of identifying high proper
motion stars in stacked image catalogs. We've performed this research in hopes of using these new methods to identify previously unidentified brown dwarves in existing catalogs.
Brown dwarves are 'failed stars' which due to their low mass were unable to sustain a fusion reaction in their cores. Their inactive cores mean they emit very little light 
making them inherently challenging to identify, meaning much is unknown about their origins and population in our cosmic neighborhood.

## Python Tools/Libraries we will be using

The main python libraries we will be taking advantage of are Pandas, Dask, and LSDB. The LSDB python library (documentation found at https://lsdb.readthedocs.io/en/stable/) uses
the LSDB framework to analyze large astronomical catalogs, which is useful for the methods we will be implementing to search for high proper motion stars (HPMS). Pandas 
(https://pandas.pydata.org/docs/) is a data analysis library, it's methods of filtering large datasets will be useful for working with our catalogs. Lastly, Dask 
(https://docs.dask.org/en/stable/) is built on pandas and also works with data analysis, but it allows for parallel and distributed computing methods which Pandas does not allow
for. given the size of the catalogs we intend to work with, this is crucial to performing larger than memory computations or saving time if this code were to be ran on
multiple systems.

## Basic Overview of Folders/Files

The [Jupyter Notebooks](https://github.com/JoaoDPassos/astrophysics/tree/main/Jupyter%20Notebooks) folder contains all of the notebooks used in this project, in no specific
order. Here is a simple explanation behind each of the Notebooks, in the order I would reccomend you follow:
- [des-gaia.ipnyb](https://github.com/JoaoDPassos/astrophysics/blob/main/Jupyter%20Notebooks/des-gaia.ipynb): Explains how to obtain catalog files and convert them to HiPSCat,
as well as how to perform catalog crossmatching. (See more on how to obtain catalog files later in this ReadMe)
- [DR2Exploration.ipynb](https://github.com/JoaoDPassos/astrophysics/blob/main/Jupyter%20Notebooks/DR2Exploration.ipynb): Explores how to graph data from specific columns in
DES_DR2.
- [OneDegExploration.ipynb](https://github.com/JoaoDPassos/astrophysics/blob/main/Jupyter%20Notebooks/OneDegExploration.ipynb): Further explores graphing values from DES_DR2,
in a conesearch focused around a known HPMS. Includes color graphs, coordinate graphs, magnitude graphs, etc. Also employs catalog filtering functions from
[catalog_filtering.py](https://github.com/JoaoDPassos/astrophysics/blob/main/catalog_filtering.py)
- [CrossMatchFIltering.ipynb](https://github.com/JoaoDPassos/astrophysics/blob/main/Jupyter%20Notebooks/CrossMatchFIltering.ipynb): This is one of the more important notebooks
as in it I have outlined the entire pipeline and explained the algorithm we created to filter for systems where there may contain HPMS. At the end of the
notebook we see that the algorithm identified the known HPMS in the dataset.
- [xmatch_dask.ipynb](https://github.com/JoaoDPassos/astrophysics/blob/main/Jupyter%20Notebooks/xmatch_dask.ipynb): This notebook takes from the previous one to convert the
entire pipeline to be compatible with the Dask python library. Furthermore, we modified parts of the algorithm such that we can filter for systems where there are n number of
aligned stars in a group to some error value.
- [proj_dist.ipynb](https://github.com/JoaoDPassos/astrophysics/blob/main/Jupyter%20Notebooks/proj_dist.ipynb): This notebook takes from the previous to plot the distribution
of maximum distances from a line projection. At one point we felt this would be a useful filtering parameter, but this way of filtering has many faults and so we decided to
not explore it further.
- [des_gaia_xmatch.ipynb](https://github.com/JoaoDPassos/astrophysics/blob/main/Jupyter%20Notebooks/des_gaia_xmatch.ipynb): This notebook contains all of the results from
our work. We crossmatched DES and Gaia to find HPMS, then ran the alogrithm of conesearches surrounding these verified stars. The graphs show stars which the algorithm
successfully or unsuccessfully identified. After this we filtered out some of the verified HPMS and plotted again, as those stars were too bright for DES.

The [Plots](https://github.com/JoaoDPassos/astrophysics/tree/main/Plots) folder contains some of the plots created throughout this process, each has titles and labels to 
describe what they are depicting.

The [algo_gaia_verified.csv](https://github.com/JoaoDPassos/astrophysics/blob/main/algo_gaia_verified.csv) is the same data found from the DES-Gaia cross match in 
des_gaia_xmatch.ipynb, but I added additional columns for Gaia found magnitudes for these stars. I inserted these columns manually because I did not have access to these columns
in my download of Gaia.
