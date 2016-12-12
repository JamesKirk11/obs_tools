## Tools for observing with ACAM

_quick_look.py_ quickly extracts the spectral traces from a single file for either 1 or 2 windows and plots them. This is pretty rough and does not 
account for curvature of the traces but is useful when selecting a spectral line for tracking with xy_centroids.py

_xy_centroids.py_ tracks the movement in x and y (given a spectral line) throughout all science FITS files within the current directory.

_ds9_finder_tool.py_ loads up either 2MASS or DSSSAO images around a certain given object and provides suitbale comparison objects in addtion to overlaying the slit.
NOTE: The manual entry of coordinates at the command line is not fully tested and may break down.

_interactive_column_check.py_ allows for testing of locations of spectral traces for bad columns, pixels and vignetting, given a flat field. **This has not been extensively tested.**

### Dependencies

Numpy
Scipy
Astropy
Matplotlib

DS9 (if using finding chart tool)
