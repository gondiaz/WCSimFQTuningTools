## Cherenkov Profiles

Follow the steps to generate the Cherenkov profiles tuning file:

**Important comments**:
- Before producing the Cherenkov Profiles, check the 2D histogram binning is reasonable in the WCSimFQTuner (_WCSimFQTunerCherenkovProfile.cc_).
  This can be done by opening any of the produced WCSim files and look at the `trg` TH2D histogram; this can be done using the  
- The tools assume that the distance units are `mm` in the simulation files, and transformed to `cm` in the step 1).
- Make sure ROOT is setup (in order to be able to `import ROOT`) 

**Step 1**: Run `merge_cprofiles_parallel.py` script: `python merge_cprofiles_parallel.py /path/to/input/files/ tr|wt [-v]`, where the input files must be named as **out_{particle}\_{momentum}_{integer}.root**. It outputs **cprofiles_merged.root** which merges the cherenkov profiles for all the simulated momenta. The second argument stands for **tr** (true) or **wt** (weighted) profiles.

**Step 2** Run `integrate_cprofiles.py`: integrates the cherenkov profiles for each momentum to produce the $I_n = \int g(s, \cos\theta) s^n ~~ (n=0, 1, 2)$. Reads the **cprofiles_merged.root** (default) created at 2. Returns a **cprofiles_integrals.root** file.

Run it using `python integrate_cprofiles.py r0max nr0bins nth0bins [-i merged_file.root] [-v]` where:
- **r0max**: maximum $r_0$ distance in mm, ie distance between PMT and track vertex (typically 5000).
- **nr0bins**: number of $r_0$ bins (typically 50).
- **nth0bins**: number of $\cos \theta_0$ bins, with $\cos \theta_0 \in [-1, +1]$ (typically 50).

**Step 3** Run `fit_cprofile_integrals.py`. For each pair $(r_0, \cos \theta_0)$, it performs a polynomial fit of the integral $I_n$ vs the logarithm of the momentum $\log p[\text{MeV/c}]$.

Run it using `python fit_cprofile_integrals.py npars [-i integrals_file.root] [--nsectmax nsectmax] [--Rmin Rmin]` where `npars` is the number of parameters of the polynomial ($I_n = \sum^{n}_{j=1} J_i~p^i$); `nsectmax` the maximum number of sections used in each fit; and `Rmin` the minimum value of the coefficient of determination for a fit result. Reads the **cprofiles_integrals.root** (default) created at 3. The output is **cprofiles_fits.root** containing the fitted parameters $j_i$ for each pair of $(r_0, \cos \theta_0)$.

The output file format is the same as the required by fiTQun, but must be renamed.

**Caution**: the script `rename_files.jl` helps to rename the files from the simulation output to the naming convention required at step 2), but it needs hardcoding, thus using it is highly not recommended.

### **Analysis**

Contains jupyter-notebooks to check distributions for each production step and compare with previous implementation.

### Check 2D histograms with WCSimFQTuner/WCSimFQTunerCherenkovProfile

Todo: add instructions
