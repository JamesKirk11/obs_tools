## Tools for observing with ACAM, now working with Python 3

*quick_look.py* 
quickly extracts the spectral traces from a single file for either 1 or 2 windows and plots them. This is pretty rough and does not 
account for curvature of the traces but is useful when selecting a spectral line for tracking with xy_centroids.py. This produces 3/4 figures (in interactive pop-up windows which need to be closed when finished with) which are: 

1) plots of the trace and extraction regions (1 figure if a single window, 2 figures if 2 windows)
2) plot of the raw spectra
3) plot of the normalised, raw spectra. Useful if there is a large brightness difference between the stars
4) plot of the ratio of the 2 spectra

*xy_centroids.py* 
tracks the movement in x and y (given a spectral line) throughout all science FITS files within the current directory. The only required arguments are the name of the object *as it appears* in the science fits headers and the -yleft argument which defines the location (in the dispersion direction) of the spectral feature needed for Y tracking.

*ds9_finder_tool.py* 
loads up either 2MASS or DSSSAO images around a certain given object and provides suitbale comparison objects in addtion to overlaying the slit.
NOTE: the path to the local ds9 installation needs to be manually replaced on line 290.
NOTE: The precession of coordinates is not fully tested and should be checked independently.

*interactive_column_check.py* 
allows for testing of locations of spectral traces for bad columns, pixels and vignetting, given a flat field. 

*27as_flat_standard_window.fit* A flat field using the standard ACAM window and taken with the 27 arcsec slit. Useful when using interactive_column_check.py

*27as_through_slit_standard_window.fit* A through slit image using the standard ACAM window and taken with the 27 arcsec slit. Useful when using interactive_column_check.py

*40as_flat_standard_window.fit* A flat field using the standard ACAM window and taken with the 40 arcsec slit. Useful when using interactive_column_check.py

*40as_flat_standard_window.fit* A through slit image using the standard ACAM window and taken with the 40 arcsec slit. Useful when using interactive_column_check.py

### Dependencies

Numpy
Scipy
Astropy
Matplotlib

If using ds9_finder_tool.py:
DS9 v7.4+ for full compatibility
Astroquery
