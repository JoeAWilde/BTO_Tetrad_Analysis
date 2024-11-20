# BTO_Tetrad_Analysis
 
This respository contains the work flow to analyse the relaive abundance of pheasants within 2km^2 tetrads in Berkshire, Cornwall, Devon, and Hertfordshire.
All of the code needed can be found within the "Code" folder. The first two scripts extract the per-tetrad protected area coverage and habitat coverage, and
bind all of this in the dataframe for each tetrad. The third script runs an occupancy-abundance model using stan (the stan script can be found in 
'Code/Stan scripts/'), and outputs the results. 

All of the data needed to run these scripts can be found in the 'Data' folder. All outputs should be set up to save to the 'Output' folder. 

This respository was written by Luke Ozsanlav-Harris and Joe A. Wilde, and was created and published by Joe A. Wilde.