# Set-up
This project was run using an old version of the pipeline and USEARCH 9.0  
The main difference is the later pipeline removes primers before estimating error rate in the preprocessing step.  
The newer versions of the preprocessing scripts use the final two command line variables to set the length of the forward and reverse primers. Setting both of these to 0 should give the same results as the old verion of these scripts, but this is untested.

Removal of the primers can still be specified in the UPARSE step.
