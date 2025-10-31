# BTO_Tetrad_Analysis
 
This respository contains the work flow to analyse the relaive abundance of pheasants within 2km^2 tetrads in Berkshire, Cornwall, Devon, and Hertfordshire.
All of the code needed can be found within the "Code" folder. The first two scripts extract the per-tetrad protected area coverage and habitat coverage, and
bind all of this in the dataframe for each tetrad. The third script runs an occupancy-abundance model using stan (the stan script can be found in 
'Code/Stan scripts/'), and outputs the results. 

The shapefiles of the SPAs, SACS, RAMSAR, and SSSIs could not be uploaded to this repository due to their size, but can be downloaded at:
https://jncc.gov.uk/our-work/uk-protected-area-datasets-for-download/

The UKCEH habitat rasters are also too large to upload to this repository, but can be downloaded from:
https://www.ceh.ac.uk/data/ukceh-land-cover-maps

Once downloaded, these shapefiles and rasters should be placed into their respective folders in the 'Data' folder. The scripts should then pull them from
these locations without issue.

All outputs should be set up to save to the 'Output' folder. 

This respository was written by Luke Ozsanlav-Harris and Joe A. Wilde, and was created and published by Joe A. Wilde..