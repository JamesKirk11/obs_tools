from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.ndimage import rotate

parser = argparse.ArgumentParser(description='Rotate images')
parser.add_argument("fits_files",help="Which fits files to load?",type=str)
parser.add_argument("-deg","--degrees",help="Define the required rotation, in degrees. Note this is a clockwise rotation!",type=int)
args = parser.parse_args()

file_names = np.loadtxt(args.fits_files,dtype=str)

for i,f in enumerate(file_names):
    
    print(f)
    
    fits_file = fits.open(f)

    image = fits_file[0].data
    
    # rotated_image = image.transpose()
    rotated_image = rotate(image, angle=args.degrees)
    
    rotated_image = np.rot90(image,k=360-(args.degrees//90))
    
    # numpy.rot90
    
    if i == 0:
	    plt.figure()
	    plt.imshow(rotated_image,vmin=np.median(rotated_image)*0.95,vmax=np.median(rotated_image)*1.05)
	    plt.show()
    # raise SystemExit
    
    header = fits_file[0].header
    
    hdu = fits.PrimaryHDU(data=rotated_image,header=header)
    
    new_filename = f[:-5]+'_rot.fits'
    
    hdu.writeto(new_filename,overwrite=True)

    fits_file.close()
