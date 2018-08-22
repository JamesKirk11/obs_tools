### Example ACAM data for testing quick_look.py and xy_centroids.py

# quick_look.py

In this case, running:

> python quick_look.py r2706800.fit

will load in file r2706800.fit and automatically locate the spectra of the target (left-hand bright spectrum) and comparison (right-hand bright spectrum) and perform simple tracing and aperture photometry, with the aperture and background regions given by the green lines in plot window #1.

If the script is struggling to locate the desired spectra, approximate locations can be parsed using the "-x ..." option.

Following this, plot window #2 shows the extracted spectra for both stars. In this example, the line at an x-coordinate (Y pixel) of 1108 looks like a nice line to perform guiding on.

Plot window #3 shows trimmed spectra, which have been normalised for easier comparison between the 2 spectra. This is not necessary in this example but is useful when the 2 stars are of differing brightness.


# xy_centroids.py

Using the absorption line at a y-pixel of 1108 we can now run xy_centroids.py as:

> python xy_centroids.py WASP-93 -yleft 1108

This is equivalent to:

> python xy_centroids.py WASP-93 -yleft 1108 -ref r2706800.fit

as we have assumed that the first file within the directory is a spectrum of 'WASP-93' is the reference image. If this is untrue or you want to use a different reference image, "-ref" should be changed to the desired reference image.

The script then plots:
- a fit to the target's spectrum in x (plot window #1)
- a fit to the absorption line at y~1108 for the target (plot window #2)
- a fit to the comparison's spectrum in x (plot window #3)
- a fit to the absorption line at y~1108 for the comparison (plot window #4)
- the x positions of the target (top panel) and comparison (bottom panel) compared to the reference image (plot window #5)
- the y positions of the abosorption line at y~1108 for the target (top panel) and comparison (bottom panel) compared to the reference image (plot window #6). Note that in this case, there are some outliers in the top panel that show a delta y of ~40, this are due to failed fits and shouldn't be trusted, this can be seen by the fact that the comparison shows a delta y of < 0.5 pixels in these frames. Playing around with "-dyline" (the number of pixels over which to perform the fit) can sometimes bring these poor fits back into line. Also choosing a different absorption line to perform the fit on can be a good way to check discrepant points.
- the maximum counts recorded in the spectra (plot window #6), which is useful to make sure the counts are <40k.
- the FWHM of the traces (plot window #7), which is useful to monitor changes in the seeing

Note: it is useful to check the x & y positions in the red and blue ends along with the centre of the spectrum as sometimes a big shift in the red might not be replicated in the blue, in which case we might rethink a guiding correction. This will also reveal affects due to rotation.
