# P28_codes

pro28createcsv_p3.py runs using Py3, remaining codes run using Py2. Tested that they run on my machine, if you can't get the python packages working give me a shout.




pro28createcsv_p3.py
- creates the csv overview file of the P28 spectra (based on your original script)


pro28statistics.py
- generates the statistics plots.
- made in 1st year of PhD, so apologies for the "quality" of the code
- code asks a series of y/n prompts to make typical plots
   - filter for good-quality spectra
   - work with rejected spectra (i.e bad, once filtered)
   - write out csv filtered by spectral type i.e. O-type stars
   - filter by spectral type (i.e. O-type) for subsequent plots (useless for barchart plot)
   - plot SNR line stats (U/B/V/R/I: count vs S/N)
   - No. observations per unique star
   - bar chart plot


p28_sky_plotting.py
- plots the galactic distribution of stars
- CP'd the quality filter code in, filters for good P28 spectra.
- Some extra CP'd stuff, unused
- color seems to be busted since I last used it. Remember messing around with many colorbars when doing P28 vs EDIBLES distrib plot.
- for the magnitude-icon size it looks like a do a quick renormalisation of Mv (for pythons sake), then clip the limits (for plottings sake), then feed the Mv values into the circle size which I believe is a radius


blend_depth_distributions.py

- scatter plot + histograms
- has a chunk of manual figure size stuff due to 2x scatter+hists layout, shout if its too confusing to copy
- do_plots() is the function you want!
- scatter plot is trivial
- histograms were manually defined in terms of bins and limits, sorry! Tweak until pretty!



Happy plotting!
