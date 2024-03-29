# Cosmic Birefringence

This code reproduces the results of Eskilt & Komatsu (2022) where Planck and WMAP polarization data are jointly sampled
to measure the cosmic birefringence angle. The details behind the method of analysis can be found in 

* Y. Minami et al., Progress of Theoretical and Experimental Physics 2019, 083E02 (2019), arXiv:1904.12440 [astro-ph.CO]
* Y. Minami and E. Komatsu, Progress of Theoretical and Experimental Physics 2020, 103E02 (2020),
arXiv:2006.15982 [astro-ph.CO]
* Y. Minami and E. Komatsu, Phys. Rev. Lett. 125, 221301 (2020), arXiv:2011.11254 [astro-ph.CO]
* Diego-Palazuelos et al., Phys. Rev. Lett. 128, 091302 (2022), arXiv:2201.07682 [astro-ph.CO]
* J. R. Eskilt, A&A 662, A10 (2022), arXiv:2201.13347 [astro-ph.CO]

# How Do I Run It?

There is a lot you need before being able to reproduce the results. Here, we outline all the steps that need to be done.

## 1) Getting the data
First, you need the Planck and WMAP maps. The detector split maps of Planck Data Release 4 (NPIPE) can be found on NERSC.
Unfortunately, it does not seem like detector split maps are available on the Planck Legacy Archive.

The WMAP maps are available here: https://lambda.gsfc.nasa.gov/product/wmap/dr5/maps_da_r9_i_9yr_get.html

Second, you need beam transfer functions. For NPIPE, they are also avalable on NERSC.
And for WMAP, they can be found here https://lambda.gsfc.nasa.gov/product/wmap/dr5/beam_xfer_get.html

Third, you need the masks. They can be found here: https://drive.google.com/file/d/1uZzBdv4eICmgWu7Kw7BXF6TioFY4S6FH/view?usp=sharing
Put them in the `generate_observed_power_spectra/masks` folder. The sky fraction f_sky is already produced in the `f.npy` file. But feel free to double check.

## 2) Install PolSpice

Then you need to do the hardest part: Install PolSpice. http://www2.iap.fr/users/hivon/software/PolSpice/
I used version v03-07-03.

## 3) Get observed power spectra, psi_ell and LCDM spectra

Now, you need to run PolSpice. This is done through the file `make_cl_files_polspice.py` which runs
through all combinations for both masks.

I highly recommend that you first only analyse HFI. This runs faster and is easier to debug.
To do that, just switch `from hfi_lfi_wmap_eb import maps_param` to `from hfi_eb import maps_param` in `make_cl_files_polspice.py`.
This needs to be done in `beam_corrected_lcdm_spectra.py` and `correct_format_observed_cl.py` as well.
These parameter files can be found in the `parameter_files` folder. You need to give the parameter file
the right path to the maps.

Once you have the observed power spectra in .dat files, you can generate psi_ell through the file `calculate_psi_ell.py`

Then you can run `beam_corrected_lcdm_spectra.py` and `correct_format_observed_cl.py`. These files
generate the beam smoothed LCDM spectra and the observed spectra in file formats that `cb.py` can read.
These files are stored in the `pre_proc` folder.

## 4) Run it!

Now you are ready to run it!

The file `run_cosmic_birefringence_analysis.py` runs the analysis on the two masks.
As mentioned, please run HFI-only first. Doing Planck+WMAP takes time,
so make sure that HFI-only works first.

Please also read the comments in `run_cosmic_birefringence_analysis.py` to get a better
understanding of what the different parameters in the `params` variable do.

`cb.py` outputs a chain file in the `chains` folder that you are free to do all your analysis on.
It also creates a corner plot that can be found in `img/`. If you let the analysis run on both masks,
it outputs a file `results.py` that gives the mean and error bars of all sampled parameters.

## Is something not working?
Please let me know if you run into any problems or if anything is unclear!

You can either open a GitHub issue, or contact me on j.r.eskilt@astro.uio.no

## Citation

Feel free to use the code as you see fit, but if you use it for published results, please cite
* J. R. Eskilt and E. Komatsu, Phys. Rev. D 106, 063503 (2022), arXiv:2205.13962 [astro-ph.CO]
