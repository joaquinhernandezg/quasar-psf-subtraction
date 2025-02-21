# quasar-psf-subtraction

For development please install with ``pip intall -e .``

# HSC config file
Credentials for accessing HSC data must be stored in a 2-line .txt file with username and password. This filename will be set in the configuration .yml file.

# Galfit multiband

Install [GALFITM](https://www.nottingham.ac.uk/astronomy/megamorph/).

Make sure you have the ``galfitm``script in your shell.

# UNAGI

You need the unagi package to query HSC data, but the original version does not work for querying PDR3 data. This forked version of the package does

https://github.com/joaquinhernandezg/unagi.git

Please install clone and install this version.



# Command Line Scripts

These are installed in your shell when installing the software

## ``make_config_file``

Run this script to create a default configuration file in the current directory.

## ``filter_catalog_hsc_footprint``

Run this script to filter a catalog with the entries inside a given HSC SSP footprint

### Arguments

- ``input_catalog``: filename of the catalog to filter
- ``hsc_release``: can be pdr3_wide, pdr3_dud, pdr2_wide, pdr2_dud
- ``ra_column``: name of the ra column
- ``dec_column``: name of the dec column
- ``output_catalog``: name of the output catalog to save, optional.
- ``catalog_format``: format of the input catalog. Uses astropy.table formats.
- ``radec_units``: units of the ra dec column. hourangle or deg.


## ``fit_catalog``
Runs the galfit fitting for all sources in a given catalog reading a configuration file.
Recomended before running the fitting routine for filtering sources that are already out of the footprint.

### Arguments

- ``config_file``: name of the yml file that stores the configuration for the fitting

# Configuration file entries


## General parameters

- `project_name`: "Substract quasar" # Just a project name
- `catalog_origin`: "SDSS DR16 QUASARS" # Where the catalog that will be used comes from
- `comments`: "Substract quasar from HSC image" # Just a comment

- `catalog_filename`: "data/catalogs/sdss_dr16_quasars.csv" # filename of the catalog to fit
- `catalog_format`: "csv" # format of the catalog. Uses astropy.table formats.
- `ra_column`: "ra" # name of the ra column in the catalog
- `dec_column`: "dec" # name of the dec column in the catalog
- `targname_column`: null # name of the target name column. When null creates a SDSS Jhhmmss.ss+/-ddmmss.ss name

- `ra_unit`: "deg" # units of the ra. Can be deg or hourangle
- `dec_unit`: "deg" # units of the dec. Can be deg or hourangle


- `output_dir`: "data/results/" # where the results of this run will be stored

- `n_processes`: 20 # number of paralel processes to run. The more the faster.
- `fit_first`: 10 # just fit the first entries in the catalog for testing. -1 to fit everything,

### Query parameters
- `hsc_query`: # this block defines the data query parameters
  - `output_hsc`: "hsc_data" # directory where the HSC cutouts will be stored.
  - `data_release`: "pdr2" # data release, pdr1, pdr2, pdr3
  - `rerun`: "pdr2_wide" # rerun to use e.g. pdr2_wide, pdr2_dud, pdr3_wide, pdr3_dud
  - `config_file_credentials`: "config_file.txt" # file with the HSC credentials
  - `filters`: "all" # filters to query, do not change
  -` cutout_size`: 5.0 # arcsec # size of the cutput to query in arcseconds
  - `pixel_scale`: 0.262 # arcsec/pix # plate scale of the cutouts, do not change.

### GALFIT settings
- `galfit`: # this block defines the galfit fitting parameters
  - `substract_sky`: True # Calculates and subsctract sky previous to the fitting
  - `fit_sky`: True # Fit the sky with galfit. 
  - `r_mask_central_source`: 10.0 # pix. Always unmask the central source for fitting
  - `only_mask_central_source`: False # if True Only fit the central pixels. If false the whole image.
  - `n_sigma_mask_sources`: 3.0 # when fitting the whole image, mask pixels above n sigmas
  - `run_single_band`: True # run the single band fitting? Do not turn off
  - `run_multi_band`: True # run the multiband fitting? Do not turn off
  - `plot`: True # make the plots? Recomened


# Ouputs

Each source output will be written in a directory ``output_dir/targname`` following this structure

- ``output_dir/``
    - ``config.yml``: copy of the configuration file used to run the code
    - ``fit_results.pdf``: summary pdf with the plots for all source fitted succesfully
    - ``log_fits.txt``: table with the fit status of each source. The status says then there was no data for this source, the fitting failed, or the fitting run sucesfully.
    - ``targname1/``
        - ``galfit_multiband/``: 
            - ``config.galfit/``: galfit configuration file for the multiband fitting
            - ``output.fits/``: galfit output block for the multiband fitting
        - ``galfit_singleband/``:
            - ``config_{band}.galfit/``: galfit configuration file for the band {band}
            - ``output_{band}.fits/``: galfit output block for the band {band}

        - ``cutout_{band}.fits``: cutout for the band {band}
        - ``mask_{band}.fits``: mask for the band {band}
        - ``psf_{band}.fits``: psf model for the band {band}


# Example

In the ``example/`` directory there is a catalog called ``double_quasars_makarov+2023.fits``. 

You can run 

``filter_catalog_hsc_footprint double_quasars_makarov+2023.fits pdr2_wide RAdeg DEdeg --radec_units deg --catalog_format fits ``

inside the ``example/`` dir to produce a catalog filtering sources inside the pdr2_wide footprint. This reduces the catalog size from 3140 to 309 sources (saving a lot of time).

Then you can run 

``quasar_subtract_from_config config_double_quasars_makarov_2023.yaml``

to run the pipeline after setting ``config_file_credentials`` inside the configuration file with you hsc credentials.