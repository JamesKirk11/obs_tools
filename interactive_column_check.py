from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Plot 1D spectra interactively to check for bad pixels/columns.')
parser.add_argument("fits_file",help="Which fits file to load?",type=str)
parser.add_argument("-n","--nclicks",help="Number of mouse clicks (and hence number of spectra) to be plotted, default = 1. If using this, NEEDS TO BE RUN INTERACTIVELY THROUGH IPYTHON!",type=int,default=1)
parser.add_argument("-pw","--pixel_width",help="Define width in pixels over which to perform the collapse, default = 50",default=50,type=int)
parser.add_argument("-d","--star_distance",help="Define distance between stars on chip, in arcminutes",type=float)
parser.add_argument("-y","--y_collapse",help="if wanting to collapse in y rather than x, default=x",action="store_true")
parser.add_argument("-px","--pixel",help="use this to input x pixel position, overriding use of mouse click",type=float)
parser.add_argument("-s","--slit_width",help='use to define which slit is used, default = 40"',type=int,default=40)
parser.add_argument("-w","--window",help='use to overplot window in red, default = 1 2148 1 2500"',nargs='+',default=[1,2148,1,2500],type=int)
parser.add_argument("-inst","--instrument",help="Which instrument? Needed for field of view & slit width and length. Default = ACAM",default='ACAM')
args = parser.parse_args()

f = fits.open(args.fits_file)
if len(f) == 2:
	data = f[1].data
else:
	data = f[0].data
	
if args.instrument == 'EFOSC':
	arcsec_per_pix = 0.12 # arcsec/pixel in manual
	print(len(data))
	if len(data) == 2060:
	    pix_scale = arcsec_per_pix
	else:
	    pix_conversion = 2060/1030. # This is needed because the calibration file used was done with 2x2 binning
	    pix_scale = pix_conversion*arcsec_per_pix
	y_range = [1030/2-7.5/pix_scale,1030/2+7.5/pix_scale]
	x_range = [0,1030] 
	pixel_width = args.pixel_width/2.
	
if args.instrument == 'ACAM':
	pix_scale = 0.253 # arcsec/pixel in x, in y it is 0.254
	pixel_width = args.pixel_width

colours = ['r','g','b','k','c','m']

def plot_slit(fits_file,slit_x_range=[279,1752],slit_y_range=[1247+800,1132+800]):
	# Check what window is enabled in fits file to allow for accurate plotting of slit
	hdr = f[0].header
	for i in range(1,5):
		window, enabled = hdr['WINSEC%d'%i].split()
		if enabled == 'enabled':
			break

	# window is currently a string, need to convert to useable list
	newstr = window.replace("]","").replace("[","").replace(":",",")
	window = list(map(int,newstr[:-1].split(",")))
	
	data_shape = np.shape(f[1].data)
	y_range = [slit_y_range[0] - window[2],slit_y_range[1] - window[2]] # use this to line up slit on fits image
	x_range = [slit_x_range[0] - window[0],slit_x_range[1] - window[0]]
	
	return window, x_range,y_range

if args.slit_width == 27 and args.instrument == 'ACAM':
	window,x_range,y_range = plot_slit(f)
	
elif args.slit_width == 40 and args.instrument == 'ACAM': # Using a 40" slit, which we don't have empirical data for so I'm assuming this is well centred
	x_centre = np.shape(data)[1]/2.
	y_centre = np.shape(data)[0]/2.
	hdr = f[0].header
	for i in range(1,5):
		window, enabled = hdr['WINSEC%d'%i].split()
		if enabled == 'enabled':
			break
	
	newstr = window.replace("]","").replace("[","").replace(":",",")
	window = list(map(int,newstr[:-1].split(",")))
	
	# unwindowed slit ranges
	x_range_unwin = [120,1904]
	y_range_unwin = [1855,2014]
	
	x_range = [x_range_unwin[0] - window[0],x_range_unwin[1] - window[0]]
	y_range = [y_range_unwin[0] - window[2],y_range_unwin[1] -window[2]]
	
else:
	pass
	
if args.instrument == 'ACAM':
    x_win = [args.window[0],args.window[1]]
    y_win = [args.window[2]-800,args.window[3]-800]

if args.instrument == 'EFOSC': # We don't window EFOSC as there's no advantage in the readout times
	x_win = [0,1030]
	y_win = [0,1030]
	

if args.pixel is None: # Only used if position of target is chosen with mouse click
    plt.ion()

# PLOT CHIP / FITS FILE
fig1 = plt.figure(figsize=(16,14))
if args.instrument == 'ACAM':
    plt.imshow(np.log10(data),cmap='gray')
if args.instrument == 'EFOSC':
    plt.imshow(data,cmap='gray')

plt.xlim(0,np.shape(data)[1])
plt.ylim(0,np.shape(data)[0])

# OVERLAY SLIT
plt.plot(x_range,[y_range[0],y_range[0]],color='g',lw=3)
plt.plot(x_range,[y_range[1],y_range[1]],color='g',lw=3)
plt.plot([x_range[0],x_range[0]],y_range,color='g',lw=3)
plt.plot([x_range[1],x_range[1]],y_range,color='g',lw=3)

# Plotting window
plt.plot(x_win,[y_win[0],y_win[0]],color='r',lw=3)
plt.plot(x_win,[y_win[1],y_win[1]],color='r',lw=3)
plt.plot([x_win[0],x_win[0]],y_win,color='r',lw=3)
plt.plot([x_win[1],x_win[1]],y_win,color='r',lw=3)

if args.pixel != None:
	x = np.array([[args.pixel,np.mean(y_range)]],dtype=int)
	args.nclicks = 1

else:
	# COLLECT POSITIONS OF MOUSE CLICKS
	x = plt.ginput(args.nclicks)

# Find the position of the comparison relative to the first click/chosen pixel
if args.star_distance != None:
	x_comp = [x[i][0] + args.star_distance*60/pix_scale for i in range(args.nclicks)]
	print("x position of comparison = %.2f"%x_comp[0])
	if args.instrument == "EFOSC":
		print("x position of centre of slit = ",(x_comp+x[i][0])/2.)
		if len(data) != 2060:
		    print("x position of centre of slit (unbinned) = ",(x_comp+x[i][0]))


# Plot frame with selected regions overlaid
for i in range(args.nclicks):
	x_pos = x[i][0]
	y_pos = x[i][1]
	
	plt.scatter(x_pos,y_pos,color=colours[i],marker='x')

	if args.star_distance != None:
		plt.scatter(x_comp[i],y_pos,color=colours[i],marker='+')
	
	if args.y_collapse:
		plt.axhline(y_pos + pixel_width,color=colours[i])
		plt.axhline(y_pos - pixel_width,color=colours[i])
	
	else:	# Collapsing in x and plotting in y
		plt.axvline(x_pos + pixel_width,color=colours[i])
		plt.axvline(x_pos - pixel_width,color=colours[i])

		if args.star_distance != None:
			plt.axvline(x_comp[i] + pixel_width,color=colours[i],ls='--')
			plt.axvline(x_comp[i] - pixel_width,color=colours[i],ls='--')

plt.xlim(0,np.shape(data)[1])
plt.ylim(0,np.shape(data)[0])
plt.show()

spec_ratio = []
fig2 = plt.figure()
for i in range(args.nclicks):

	if args.y_collapse:
		spectra = data[x[i][1]-pixel_width:x[i][1]+pixel_width]
		plt.plot(spectra.sum(axis=0),color=colours[i])

	else: # Collapsing in x and plotting in y
		spectra = data[:,int(x[i][0]-pixel_width):int(x[i][0]+pixel_width)]
		plt.plot(spectra,color=colours[i])
		
		if args.star_distance != None:
			spectra_comp = data[:,int(x_comp[i]-pixel_width):int(x_comp[i]+pixel_width)]
			spec_ratio.append(spectra.sum(axis=1)/np.array(spectra_comp,dtype=float).sum(axis=1))
			plt.plot(spectra_comp,color=colours[i+1])

plt.ylabel('Counts')
slit_pos_x = [279,1722]
slit_pos_y = [1247,1132]
if args.y_collapse:
	plt.xlabel('X position (pixels)')
	plt.title("Collapsing in y, plotting in x")
else:
	plt.title("Collapsing in x, plotting in y")
	plt.xlabel('Y position (pixels)')

plt.show()

# Collapsing in x and plotting in x
fig3 = plt.figure()
for i in range(args.nclicks):

	if args.y_collapse:
		spectra = data[x[i][1]-pixel_width:x[i][1]+pixel_width]
		print(spectra.sum(axis=0))
		plt.plot(spectra.sum(axis=0),color=colours[i])

	else:
		plt.plot(np.arange(x[i][0]-pixel_width,x[i][0]+pixel_width),spectra.sum(axis=0),color=colours[i])
		
		if args.star_distance != None:
			plt.plot(np.arange(x_comp[i]-pixel_width,x_comp[i]+pixel_width),spectra_comp.sum(axis=0),color=colours[i],ls='--')

plt.ylabel('Counts')
plt.xlabel('X position (pixels)')
plt.title("Plotting in x")
plt.show()

# Plot the ratio of counts from each column
if args.star_distance != None:
	plt.figure()
	for i in range(args.nclicks):
		# Plotting in y
		plt.plot(spec_ratio[i],color=colours[i])
plt.xlabel('Y position (pixels)')
plt.ylabel('Ratio of counts')
plt.title('Ratio of counts')
plt.show()

