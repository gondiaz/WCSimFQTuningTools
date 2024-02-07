## Angular response

The angular response for a PMT is the photon collection efficiency as a function of the angle of incidence. 
In fiTQun, the angular response factor is $\epsilon(\eta)$, and to compute it we follow the procedure:

1) Compute the most general incidence factor $J(R, \eta) = \Omega(R) T(R) \epsilon(\eta)$ where R is the distance between the emission point and $\eta$ is the incidence angle. This factor is computed as an 2D histogram using the "electron bomb" simulation files also used in the computation of the scatterring table. The output is a `.root` file containing three 2D histograms, which must be the same in the absence of bias: using all the PMTs, using only side PMTs and using the caps PMTs. The histograms are normalized such that $J(R, 0) = 1$ (normal incidence). Thus the computed 2D histogram is $J(R, \eta) / J(R, 0) = \epsilon(\eta) / \epsilon(\eta=0)$

    First of all remember to setup ROOT (in order to be able to `import ROOT`), then run it through:
    `python compute_angular_responses.py [-v] [--nbins Rnbins etanbins] [--vaxis 0|1|2 (default 2)] [--zedge activeLength] [--redge activeRad] [--wcsimlib [WCSIMLIB]] indir`

    where `nbins` are the number of bins for each dimension of the histogram; `vaxis` is the vector component of the vertical axis (ie 0: X, 1: Y, 2: Z); `zedge` and `redge` represent the vertical (half) and horizontal active area of the detector in cm (ie half-length and radius); `wcsimlib` is the directory of the WCSim library; and `indir` is the directory containing the input files.


2) Compute the angular response function $\epsilon (\eta)$. Taken into account that $\Omega (R=0) = 1,~ T(R=0) = 1$, the angular response is computed using the values of the first radial bin in the 2D histogram. The values are fitted using a polynomial function and saved to the output file. 

    Run it through: `python fit_angular_response.py [-v] [-i infile (default angular.root)] polydeg`

    where `infile` is the input file and `polydeg` is the degree of the polynomial used to fit the angular response. To check the fit output (plot), run the script using the `-v` flag.