## Charge PDF

First of all remember to setup ROOT (in order to be able to `import ROOT`).

The procedure involves the following steps:

1) Create the 2D histograms of true/predicted charge ($\mu$) vs measured charge ($q$) and fit the unhit probability vs $\mu$.

    Run it through: `python create_2Dhistos_and_unhitP.py [-v] [--qbins qbins-filename] [--wcsimlib [WCSIMLIB]] indir`
    
    where `v` activates verbosity (recommended to see the unhit probability fit); `qbins` is a file containing the bin edges desired for the $q$ variable; `wcsimlib` is the directory of the WCSim library; and `indir` is the directory of the simulated data for all the $\mu$ values. The files in `indir` **must** be named **out_{$\mu$}_{index}.root**.

2) Perform the $q$ vs $f(q|\mu)$ fits.

    Run it through: `python fit_charge_pdfs.py [-v] [--infile [INFILE]] [--npars NPARS [NPARS ...]] [--qranges QRANGES [QRANGES ...]]`
    
    where `v` activates verbosity; `infile` point to the file produced at step 1 (default); `npars` is the number of parameters for each range of $q$; `qranges` defines the boundaries of the $q$ ranges.

    The meaning of the last two parameters is the following: if `qranges` is `0 10 20` and `npars` is `5 6`, then the fits for $q$ values between `0` and `10` are `5` parametric and for $q$ values between `10` and `20` are `6` parametric.
