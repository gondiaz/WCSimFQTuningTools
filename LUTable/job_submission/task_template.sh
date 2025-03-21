#!/bin/bash

source $HOME/.bashrc
setup_root

python $HOME/Software/WCSimFQTuningTools/LUTable/get_charge_distributions.py -v --rmax 200 --Nmax 500 --cprofile $T2K_LUSTRE/CProfiles/e-/tuning/cprofiles_fits.root --outfile outfilename --infiles input_files
