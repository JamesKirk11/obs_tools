from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Plot 1D spectra interactively as a first look. Very rough, assumes the spectra are the birghtest things on the chip and does not account for curvature of the traces.')
parser.add_argument("fits_file",help="Which fits file to load?",type=str)
parser.add_argument("-sw","--spectral_width",help="Define width in pixels over which to sum the spectra, default = 50",default=50,type=int)
parser.add_argument("-bw","--background_width",help="Define width in pixels over which to estimate the background, default = 60",default=60,type=int)
parser.add_argument("-bo","--background_offset",help="Define offset in pixels from the spectra over which to estimate the background, default = 80",default=80,type=int)
parser.add_argument("-ylow","--lower_y_cutoff",help="Define lower cutoff in y pixels at which not to plot the ratio of spectra, default = 500",default=500,type=int)
parser.add_argument("-yup","--upper_y_cutoff",help="Define upper cutoff in y pixels at which not to plot the ratio of spectra, default = -150",default=-150,type=int)
args = parser.parse_args()

f = fits.open(args.fits_file)
print len(f)
nwindows = len(f) - 1

if nwindows == 0: # Using EFOSC data:
    data1 = f[0].data
if nwindows == 1:
    data1 = f[1].data
if nwindows == 2:
    data2 = f[2].data 

nrows,ncols = np.shape(data1)

if nwindows <= 1:
    trace1 = np.argmax(data1[nrows/2][20:ncols/2])
    trace2 = np.argmax(data1[nrows/2][ncols/2:-20])

    #x_pos1 = 337
    #x_pos2 = 900
    x_pos1 = trace1 + 20
    x_pos2 = trace2 + ncols/2

if nwindows == 2:
    trace1 = np.argmax(data1[nrows/2][20:-20])
    trace2 = np.argmax(data2[nrows/2][20:-20])
    
    x_pos1 = trace1 + 20
    x_pos2 = trace2 + 20
    
y_pos1 = y_pos2 = nrows/2

fig1 = plt.figure()
plt.imshow(data1,cmap='hot')

# Collapsing in x and plotting in y
plt.axvline(x_pos1 + args.spectral_width,color='g')
plt.axvline(x_pos1 - args.spectral_width,color='g')
plt.axvline(x_pos1 + args.background_offset+args.background_width//2,color='g',ls='--')
plt.axvline(x_pos1 + args.background_offset-args.background_width//2,color='g',ls='--')
plt.axvline(x_pos1 - args.background_offset-args.background_width//2,color='g',ls='--')
plt.axvline(x_pos1 - args.background_offset+args.background_width//2,color='g',ls='--')

if nwindows <= 1:
    plt.axvline(x_pos2 + args.spectral_width,color='g')
    plt.axvline(x_pos2 - args.spectral_width,color='g')
    plt.axvline(x_pos2 + args.background_offset+args.background_width//2,color='g',ls='--')
    plt.axvline(x_pos2 + args.background_offset-args.background_width//2,color='g',ls='--')
    plt.axvline(x_pos2 - args.background_offset-args.background_width//2,color='g',ls='--')
    plt.axvline(x_pos2 - args.background_offset+args.background_width//2,color='g',ls='--')


if nwindows == 2:
    plt.xlim(0,ncols)
    plt.ylim(0,nrows)
    plt.show()
    
    fig1a = plt.figure()
    plt.imshow(data2,cmap='hot')
    
    # Plot frame with selected regions overlaid
        
    # Collapsing in x and plotting in y
    plt.axvline(x_pos2 + args.spectral_width,color='g')
    plt.axvline(x_pos2 - args.spectral_width,color='g')
    plt.axvline(x_pos2 + args.background_offset+args.background_width//2,color='g',ls='--')
    plt.axvline(x_pos2 + args.background_offset-args.background_width//2,color='g',ls='--')
    plt.axvline(x_pos2 - args.background_offset-args.background_width//2,color='g',ls='--')
    plt.axvline(x_pos2 - args.background_offset+args.background_width//2,color='g',ls='--')
    
plt.xlim(0,ncols)
plt.ylim(0,nrows)
plt.show()


def get_spectrum(imdata,selected_region,spec_width,bkg_width,bkg_off):

    spectra = imdata[:,selected_region-spec_width:selected_region+spec_width]
    background_1 = imdata[:,selected_region-bkg_off-bkg_width//2:selected_region-bkg_off+bkg_width//2]
    background_2 = imdata[:,selected_region+bkg_off-bkg_width//2:selected_region+bkg_off+bkg_width//2]
    mean_bkg = np.column_stack((background_1,background_2)).mean(axis=1)
    spec_corrected = spectra - mean_bkg[:,None]
    
    return spectra,mean_bkg,spec_corrected.sum(axis=1)
    
def normalise(data):
    return data/data.sum()
    return (data - data.min())/(data.max()-data.min())

spec1,back1,corr1 = get_spectrum(data1,x_pos1,args.spectral_width,args.background_width,args.background_offset)
if nwindows <= 1:
    spec2,back2,corr2 = get_spectrum(data1,x_pos2,args.spectral_width,args.background_width,args.background_offset)
else:
    spec2,back2,corr2 = get_spectrum(data2,x_pos2,args.spectral_width,args.background_width,args.background_offset)

print "Peak counts (line 1) = ",np.max(corr1),"Y position = ",range(len(corr1))[np.argmax(corr1)]
print "Peak counts (line 2) = ",np.max(corr2),"Y position = ",range(len(corr2))[np.argmax(corr2)]

plt.figure()
plt.plot(corr1,color='r',label='Left')
plt.plot(corr2,color='b',label='Right')
plt.xlabel('Y pixel')
plt.ylabel('Counts')
plt.legend()
plt.show()

xaxis = range(args.lower_y_cutoff,nrows+args.upper_y_cutoff)

plt.figure()
plt.plot(xaxis,normalise(corr1)[args.lower_y_cutoff:args.upper_y_cutoff],color='r',label='Left')
plt.plot(xaxis,normalise(corr2)[args.lower_y_cutoff:args.upper_y_cutoff],color='b',label='Right')
plt.xlabel('Y pixel + %d'%args.lower_y_cutoff)
plt.ylabel('Normalised counts')
plt.legend()
plt.show()


plt.figure()
plt.plot(corr1/corr2)
plt.xlabel('Y pixel')
plt.ylabel('Left/Right')
plt.title('Ratio of spectra')
plt.show()
