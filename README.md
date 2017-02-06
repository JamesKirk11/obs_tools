## Tools for observing with ACAM

_quick_look.py_ \nquickly extracts the spectral traces from a single file for either 1 or 2 windows and plots them. This is pretty rough and does not 
account for curvature of the traces but is useful when selecting a spectral line for tracking with xy_centroids.py. This produces 3/4 figures (in interactive pop-up windows which need to be closed when finished with) which are: 

1) plots of the trace and extraction regions (1 figure if a single window, 2 figures if 2 windows)
2) plot of the raw spectra
3) plot of the normalised, raw spectra. Useful if there is a large brightness difference between the stars
4) plot of the ratio of the 2 spectra

_xy_centroids.py_ \ntracks the movement in x and y (given a spectral line) throughout all science FITS files within the current directory. The only required arguments are the name of the object *as it appears* in the science fits headers and the -yleft argument which defines the location (in the dispersion direction) of the spectral feature needed for Y tracking.

_ds9_finder_tool.py_ \nloads up either 2MASS or DSSSAO images around a certain given object and provides suitbale comparison objects in addtion to overlaying the slit.
NOTE: The manual entry of coordinates at the command line is not fully tested and may break down.

_interactive_column_check.py_ \nallows for testing of locations of spectral traces for bad columns, pixels and vignetting, given a flat field. 

_27as_flat_standard_window.fit_ A flat field using the standard ACAM window and taken with the 27 arcsec slit. Useful when using interactive_column_check.py\n
_27as_through_slit_standard_window.fit_ A through slit image using the standard ACAM window and taken with the 27 arcsec slit. Useful when using interactive_column_check.py\n
_40as_flat_standard_window.fit_ A flat field using the standard ACAM window and taken with the 40 arcsec slit. Useful when using interactive_column_check.py\n
_40as_flat_standard_window.fit_ A through slit image using the standard ACAM window and taken with the 40 arcsec slit. Useful when using interactive_column_check.py\n

### Dependencies

Numpy
Scipy
Astropy
Matplotlib

DS9 (if using finding chart tool) v7.4 for full compatibility
