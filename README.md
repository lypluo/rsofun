# README

This describes the steps to execute a standard P-model calibration. Implemented in R.

## What is this repository for?

* At this stage, this repository provides the functionality to calibrate the globally uniform and temporally constant **apparent quantum yield parameter** `q0_apparent` in the P-model. 
* The calibration method is linear least squares.
* The P-model provides a process-based prediction of photosynthetic parameters and light use efficiency as a function of only climate (temperature, VPD, atmospheric pressure) and CO2. In general, parameters are determined independently from data (see Wang Han et al., 2017 *Nature Plants*). For simulating gross primary productivity (GPP), the P-model requires vegetation greenness, or the fraction of absorbed photosynthetically active radiation (fAPAR), to be prescribed from data (see also [here](https://stineb.github.io/pmodel.html)). 
* Calibrating `q0_apparent` is motivated by the fact that there remains some unclarity about how the photosynthetic quantum yield efficiency relates to the quantity quantified by fAPAR. Hence, we aim for a calibration of this parameter, *given fAPAR data*. The parameter value is to be treated temporally constant and globally uniform, applicable for C3 photosynthesis and has the same role as our $\phi_0$ parameter in the [P-model](https://stineb.github.io/pmodel.html).
* The calibration described here as the default uses:
  - GPP based on the night-time flux decomposition (`NT_VUT_REF`) from 167 sites in the FLUXNET 2015 Tier 1 dataset. Data is read in and cleaned using XXX.
  - fAPAR quantified by MODIS FPAR MCD15A3H (Collection 6) data (4 days, 500 m). Data downloaded and interpolated to daily using XXX.
  - `q0_apparent` minimises the sum of square model errors, daily data, filtered to exclude cold days (< 5 deg C) and days with dry soil (relative soil water content < 0.5 of maximum water holding capacity). See XXX.

## Setup steps

### Summary of steps:
  1. Get forcing data
  2. Run P-model with arbitrary `q0_apparent`
  3. Read forcing data and P-model output into common a dataframe.
  4. Perform calibration.

### Get forcing data

Model forcing data is available on Imperial HPC's CX1 (`work/bstocker/data`) or can be downloaded or read from files using R scripts in repository [getin](https://bitbucket.org/labprentice/getin). All described below.

Here, we're using site-specific model forcing data. Sites used are defined in a meta information file, specific for a *simulation suite*, referring to a collection of sites that belong to a certain dataset, measurement campaign, or similar. Here, we're using the FLUXNET 2015 dataset. Therefore, we're setting:
```r
simsuite <- "fluxnet2015"
```
Respective meta information (required are columns `mysitename`, `lon`, `lat`, `elv`) is defined in the files available from the repository [getin](https://bitbucket.org/labprentice/getin): `getin/metafiles_<simsuite>/siteinfo_<simsuite>_sofun.csv`

#### fAPAR

When executing the calibration on a local machine, download the data from CX1 or process new data.

##### Downloading from CX1 
Site-scale subsets from the MODIS FPAR MCD15A3H (Collection 6) data (4 days, 500 m) for all FLUXNET 2015 Tier 1 sites are available on Imperial HPC's CX1:
`work/bstocker/data/fapar_MODIS_FPAR_MCD15A3H_fluxnet2015_gee_subset/fapar_MODIS_FPAR_MCD15A3H_<sitename>_gee_subset.csv`

Download the entire directory `fapar_MODIS_FPAR_MCD15A3H_fluxnet2015_gee_subset` and place it in `your_chosen_home/data/`. Do not modify the directory and file name structure from `your_chosen_home/data/` downwards. Specify the path in the header of file XXX.

##### Processing new 
Downloading site-scale data from Google Earth Engine and interpolating to daily data is done using `gee_subset.py` from the repository [gee_subset](https://github.com/stineb/gee_subset) by Koen Hufkens and `get_sitesubset_gee.R` from the repository [getin](https://bitbucket.org/labprentice/getin). Gapfilling and interpolation to daily values is done by filtering based on the MODIS quality flags (see `gapfill_modis.R` from the repository [getin](https://bitbucket.org/labprentice/getin)) and applying a spline to daily values.

`your_chosen_home` is the path where you chose to place the repositories getin and gee_subset and the `data` directory.

**1. Clone gee_subset**

In your shell, do:
```bash
cd your_chosen_home
git clone https://github.com/stineb/gee_subset 
```
Switch to branch `gee_stineb` by:
```bash
cd your_chosen_home/gee_subset
git checkout gee_stineb
```
To execute `gee_subset.py`, you must have a Google Earth Enginge login and authenticate yourself. To set this up, follow steps described in `setup_steps.md`.

**2. Clone getin**

In your shell, do:
```bash
cd your_chosen_home
git clone https://bitbucket.org/labprentice/getin
```
You must have a login on bitbucket and belong to the group 'labprentice' to get getin.

**3. Execute**

Change to the directory to where your local clone of the [getin](https://bitbucket.org/labprentice/getin) repository is located. In the header of the R script `get_sitesubset_gee.R` (see `MANUAL SETTINGS`), set then manual settings by hand, for example:
```r
##--------------------------------------------------------------------
## MANUAL SETTINGS
##--------------------------------------------------------------------
myhome               = "/alphadata01/bstocker/"
simsuite             = "fluxnet2015"
bundle               = "fapar"
overwrite_raw        = FALSE
overwrite_nice       = TRUE
do_plot_interpolated = TRUE
start_date           = "2000-01-01"
end_date             = "2017-11-27"
years_out            = 1980:2017   # creating SOFUN input data for these years
##--------------------------------------------------------------------
```

 Then execute it in R:
```r
setwd("your_local_path/getin")
source("get_sitesubset_gee.R")
```


#### Climate

Climate model forcing for site-scale simulations is extracted from site-specific meteo data and/or global fields (WATCH-WFDEI and CRU TS implemented), and written to CSV files containing all days of all simulation years, as well as annual ascii text files for each year. All files are stored in a pre-defined directory structure:

* CSV files: 
  `your_chosen_home/sofun/input_<simsuite>_sofun/sitedata/climate/<sitename>/clim_daily_<sitename>.csv`
* Annual ascii text files: 
  `your_chosen_home/sofun/input_<simsuite>_sofun/sitedata/climate/<sitename>/<year>/d<var>_<sitename>_<year>.txt`

All climate data for FLUXNET 2015 simulations can be read from data distributed by the FLUXNET 2015 dataset. Download the full dataset from CX1 to your local machine (alternatively, run all the following scripts on CX1). The dataset with daily data is on CX1 under:
`/work/bstocker/data/FLUXNET-2015_Tier1`. Place this directory in `your_local_path/data/`.

Processing new climate input files is done using functionality from the repository [getin](https://bitbucket.org/labprentice/getin). Clone this repo as described above. In `get_climate2.R`, adjust the header manually. For example:
```r
##--------------------------------------------------------------------
## MANUAL SETTINGS
##--------------------------------------------------------------------
myhome           = "/alphadata01/bstocker/"
simsuite         = "fluxnet2015"  
overwrite        = TRUE
overwrite_byst   = TRUE
in_ppfd          = TRUE
in_netrad        = FALSE
startyr_override = 2010
add_cru          = FALSE
add_watch        = FALSE
add_meteo        = TRUE
##--------------------------------------------------------------------
```

Then execute it in R:
```r
setwd("your_local_path/getin")
source("get_climate2.R")
```


* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact