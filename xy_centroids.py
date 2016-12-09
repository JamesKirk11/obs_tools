import numpy as npimport scipy.optimize as soptfrom astropy.io import fitsimport matplotlib.pyplot as pltimport globimport pickleimport argparseparser = argparse.ArgumentParser(description='Track the x and y positions, FWHM and maximum counts as a function of time for 2 spectral traces. Any argument with a default value does not need to be supplied at the command line.')parser.add_argument('object_name',help='Name of the object as it appears in the FITS headers')parser.add_argument('-ref','--reference_image',help='Define the reference image for which subsequent frames will be compared to. If left blank this will use the first science image by default. Quicker to define this so code does not have to check header of each file.')parser.add_argument('-n','--n_images',type=int,help="Number of images to calculate the shift for, counting from the end, e.g. '-n 10' will find shifts for last 10 images. Default is all images.")parser.add_argument('-ns','--no_save',action='store_true',help='Use this argument if not wanting to save the outputted figures. Default is to save.')parser.add_argument('-ycont','--y_cont',type=int,help='y value at which x centroid will be evaluated (for both objects), default = nrows/2.')parser.add_argument('-yleft','--y_left',type=int,help='Rough y value of spectral line that you want to use for the y centroid, for the left hand trace')parser.add_argument('-yright','--y_right',type=int,help='Rough y value of spectral line that you want to use for the y centroid, for the right hand trace. If left blank, the y_left will be used here also.')parser.add_argument('-specw','--spec_w',type=int,help='Number of pixels in x to perform x profile fit over. Default = 40',default=40)parser.add_argument('-dyline','--dy_line',type=int,help='Number of pixels in y to perform line profile fit over. Default = 100',default=100)parser.add_argument('-dycont','--dy_cont',type=int,help='Number of pixels in y to perform continuum fit over. Default = 5',default=5)parser.add_argument('-bkgoff','--bkg_off',type=int,help='Offset in x between extraction aperture and background aperture. Default = 30',default=30)parser.add_argument('-bkgw','--bkg_w',type=int,help='Number of background pixels in x to extract spectra. Default = 20',default=20)args = parser.parse_args()y_cont = args.y_conty_line_left = args.y_leftspec_w = args.spec_wdy_line = args.dy_linedy_cont = args.dy_contbkg_w = args.bkg_wbkg_off = args.bkg_offif y_line_left is None:    raise SystemError('y_left is not defined')y_line_right = args.y_rightif y_line_right is None:    y_line_right = y_line_left# obj_name must match name of the object as in FITS headersobj_name = args.object_name
# Read names of all local fits files for time being
file_names = sorted(glob.glob("r*.fit"))if args.reference_image is not None:    # We've defined a reference image, and need to find its location within the local directory        ref_frame = args.reference_image    for i,j in enumerate(file_names):        if j == ref_frame:            file_names = file_names[i:]
else: # Find first science frame    for i,f in enumerate(file_names):         hdu = fits.open(f)         header = hdu[0].header         if header['OBJECT'] == obj_name and header['ACAMMODE'] == 'SPECTROSCOPY':             ref_frame = file_names[i]             print 'Reference frame taken to be %s'%ref_frame             file_names = file_names[i:]             hdu.close()             break         hdu.close()        # Use reference image to get number of windows and window sizeref_fits = fits.open(ref_frame)nrows,ncols = np.shape(ref_fits[1].data)if y_cont is None:    y_cont = nrows/2nwindows = len(ref_fits) - 1ref_fits.close()if nwindows == 1:    # Search range for each of the spectral traces    x_beg_left = 20    x_end_left = ncols/2 - 20        x_beg_right = ncols/2 + 20    x_end_right = -20if nwindows == 2:    x_beg_left = 20    x_end_left = -20    x_beg_right = 20    x_end_right = -20    if args.n_images is None:    # Loading in all files    l = np.loadtxt(file_names,str)else:    # Load in last n files plus the reference image    nfiles = args.n_images     fits_files = [ref_frame]+file_names[-nfiles:]    l = np.loadtxt(fits_files,str)    

def extract_flux_x(frame,x_beg,x_end,y_cont,dy_cont,spec_w,bkg_w,bkg_off):
    
    frame_subset = frame[y_cont-dy_cont//2:y_cont+dy_cont//2+1,x_beg:x_end]    
    max_ind = np.median(np.argmax(frame_subset,axis=1))    
    spec_2d = frame_subset[:,max_ind-spec_w//2:max_ind+spec_w//2+1]
    
    bkg_1_2d = frame_subset[:,max_ind-bkg_off-bkg_w//2:max_ind-bkg_off+bkg_w//2+1]    bkg_2_2d = frame_subset[:,max_ind+bkg_off-bkg_w//2:max_ind+bkg_off+bkg_w//2+1]    
    bkg = np.column_stack((bkg_1_2d,bkg_2_2d)).mean(axis=1)    
    spec = (spec_2d - bkg[:,None]).sum(axis=0)    
    return x_beg+max_ind, spec
def extract_flux_y(frame,x_beg,x_end,y_line_h32,dy_line,spec_w,bkg_w,bkg_off):        frame_subset = frame[y_line_h32-dy_line//2:y_line_h32+dy_line//2+1,x_beg:x_end]    
    y =  np.arange(frame_subset.shape[0])    
    spec_2d, bkg_1_2d, bkg_2_2d = spec_regions(frame_subset,y,spec_w,bkg_w,bkg_off)    
    bkg = np.column_stack((bkg_1_2d,bkg_2_2d)).mean(axis=1)    
    spec = (spec_2d - bkg[:,None]).sum(axis=1)    
    return spec

def spec_regions(frame_subset,y,spec_w,bkg_w,bkg_off):
    max_ind = model_track(frame_subset,y)
    spec_cut = np.array([np.arange(im-spec_w//2,im+spec_w//2+1).tolist() for im in max_ind])
    spec_cut[spec_cut<0] = 0
    spec_cut[spec_cut>=frame_subset.shape[1]] = frame_subset.shape[1] - 1
    bkg_1_cut = np.array([np.arange(im-bkg_off-bkg_w//2,im-bkg_off+bkg_w//2+1).tolist() for im in max_ind])
    bkg_1_cut[bkg_1_cut<0] = 0
    bkg_1_cut[bkg_1_cut>=frame_subset.shape[1]] = frame_subset.shape[1] - 1
    bkg_2_cut = np.array([np.arange(im+bkg_off-bkg_w//2,im+bkg_off+bkg_w//2+1).tolist() for im in max_ind]) 
    bkg_2_cut[bkg_2_cut<0] = 0
    bkg_2_cut[bkg_2_cut>=frame_subset.shape[1]] = frame_subset.shape[1] - 1
    b_spec_ind = np.broadcast_arrays(y[:,None],spec_cut)
    spec_ind = [b_spec_ind[0].tolist(),b_spec_ind[1].tolist()]
    b_bkg_1_ind = np.broadcast_arrays(y[:,None],bkg_1_cut)
    bkg_1_ind = [b_bkg_1_ind[0].tolist(),b_bkg_1_ind[1].tolist()]
    b_bkg_2_ind = np.broadcast_arrays(y[:,None],bkg_2_cut)
    bkg_2_ind = [b_bkg_2_ind[0].tolist(),b_bkg_2_ind[1].tolist()]

    return frame_subset[spec_ind], frame_subset[bkg_1_ind], frame_subset[bkg_2_ind]
  def model_track(data,y):    """Function that takes in some subset of the data frame and outputs a list of the modelled maxima """    
    # Cubic fit to the     p = np.polyfit(y,np.argmax(data,axis=1),3)    mod = np.around(np.polyval(p,y)).astype(int)    return mod
def moffat_x_mod(a,r):    x = 2. * (a[1]-1.) / (a[0] * a[0])    y = 1. + ((r-a[2])/a[0])**2.    return a[3] * (x * (y ** -a[1])) + a[4]
def moffat_x(a,r,y):    mod = moffat_x_mod(a,r)    return y - mod
def moffat_y_mod(a,r):    x = 2. * (a[1]-1.) / (a[0] * a[0])    y = 1. + ((r-a[2])/a[0])**2.    return a[3] * (x * (y ** -a[1])) + a[4] * r + a[5]
def moffat_y(a,r,y):    mod = moffat_y_mod(a,r)    return y - mod

x_left = []y_left = []x_right = []y_right = []mjd = []flist = []
max_counts = []y_pos_max_counts = []fwhm_x = []fwhm_y = []
for i,f in enumerate(l):
    hdu = fits.open(f)    hdr = hdu[0].header
    if np.logical_and("ACAMMODE" in hdr.keys(),"OBJECT" in hdr.keys()):
        print f,hdr["ACAMMODE"],hdr["OBJECT"]
        if hdr["ACAMMODE"] == "SPECTROSCOPY" and hdr["OBJECT"] == obj_name:
            frame1 = hdu[1].data                        if nwindows == 2:                frame2 = hdu[2].data                max_counts.append(np.max((frame1,frame2)))
            mjd.append(hdr["MJD-OBS"])            flist.append(f)
            hdu.close()
            # Extract flux in x collapsed in y            peak_x, flux_x = extract_flux_x(frame1,x_beg_left,x_end_left,y_cont,dy_cont,spec_w,bkg_w,bkg_off)
            x = peak_x+np.arange(-len(flux_x)/2,len(flux_x)/2)
            # Model flux in x with a gaussian.
            parms_x_left = sopt.leastsq(moffat_x,np.array([6.,5.,peak_x,200000.,0.]),args=(x,flux_x))                                                   fwhm_x.append(0.25*2*parms_x_left[0][0]*np.sqrt(2**(1/parms_x_left[0][1])-1)) # 0.25 needed for arcsec conversion
            flux_y = extract_flux_y(frame1,x_beg_left,x_end_left,y_line_left,dy_line,spec_w,bkg_w,bkg_off)
            y = y_line_left - dy_line//2 + np.arange(len(flux_y)) 
            parms_y_left = sopt.leastsq(moffat_y,np.array([4.0,1.5,y_line_left,-150000.,0.,50000.]),args=(y,flux_y))              fwhm_y.append(0.25*2*parms_y_left[0][0]*np.sqrt(2**(1/parms_y_left[0][1])-1)) # 0.25 needed for arcsec conversion
            if f == ref_frame:
                xfine = np.arange(x[0],x[-1],0.001)    
                plt.figure()                plt.plot(x,flux_x,'.')                plt.plot(xfine,moffat_x_mod(parms_x_left[0],xfine))                plt.title("x centroiding for left spectrum")                plt.ylabel("flux")                plt.xlabel("x pixel")                plt.show()                              
                yfine = np.arange(y[0],y[-1],0.001)    
                plt.figure()                plt.plot(y,flux_y,'.')                plt.plot(yfine,moffat_y_mod(parms_y_left[0],yfine))                plt.title("y centroiding for left spectrum")                plt.ylabel("flux")                plt.xlabel("y pixel")                plt.show()
            if nwindows == 2:                peak_x, flux_x = extract_flux_x(frame2,x_beg_right,x_end_right,y_cont,dy_cont,spec_w,bkg_w,bkg_off)                flux_y = extract_flux_y(frame2,x_beg_right,x_end_right,y_line_right,dy_line,spec_w,bkg_w,bkg_off)            else:                peak_x, flux_x = extract_flux_x(frame1,x_beg_right,x_end_right,y_cont,dy_cont,spec_w,bkg_w,bkg_off)                flux_y = extract_flux_y(frame1,x_beg_right,x_end_right,y_line_right,dy_line,spec_w,bkg_w,bkg_off)            x = peak_x+np.arange(-len(flux_x)/2,len(flux_x)/2)            # Model flux in x with a gaussian.            parms_x_right = sopt.leastsq(moffat_x,np.array([6.0,4.0,peak_x,650000.,0.]),args=(x,flux_x))

            y = y_line_right - dy_line//2 + np.arange(len(flux_y)) 
            parms_y_right = sopt.leastsq(moffat_y,np.array([4.,1.0,y_line_right,-5000000.,0.,200000.]),args=(y,flux_y)) 
            x_left.append(parms_x_left[0][2])
            y_left.append(parms_y_left[0][2])
            x_right.append(parms_x_right[0][2])
            y_right.append(parms_y_right[0][2])
            if f == ref_frame:            
                xfine = np.arange(x[0],x[-1],0.001)
                plt.figure()                plt.plot(x,flux_x,'.')                plt.plot(xfine,moffat_x_mod(parms_x_right[0],xfine))                plt.title("x centroiding for right spectrum")                plt.ylabel("flux")                plt.xlabel("x pixel")                plt.show()                          
                yfine = np.arange(y[0],y[-1],0.001)
                plt.figure()                plt.plot(y,flux_y,'.')                plt.plot(yfine,moffat_y_mod(parms_y_right[0],yfine))                plt.title("y centroiding for right spectrum")                plt.ylabel("flux")                plt.xlabel("y pixel")                plt.show()    
        else:            print "ACAMMODE OR OBJECT NAME DO NOT MATCH"
            hdu.close()            break
    else:
        hdu.close()
x_left = np.array(x_left)y_left = np.array(y_left)x_right = np.array(x_right)y_right = np.array(y_right)mjd = np.array(mjd)
ref_i = np.where(np.array(flist)==ref_frame)
print "X_left now = ",x_left[-1],"; X_left ref = ",x_left[ref_i]print "Y_left now = ",y_left[-1],"; Y_left ref = ",y_left[ref_i]print "X_right now = ",x_right[-1],"; X_right ref = ",x_right[ref_i]print "Y_right now = ",y_right[-1],"; Y_right ref = ",y_right[ref_i]
plt.figure(figsize=(11.69,8.27))
plt.subplot(211)if args.n_images is None:
    plt.plot(mjd-mjd[0],x_left-x_left[ref_i],'k.')else:    plt.plot((mjd-mjd[0])[1:],(x_left-x_left[ref_i])[1:],'k.')plt.title("x centroid at y~"+str(y_cont)+" (Target)")plt.xlabel("JD from start of observations")plt.ylabel("x centroid - reference")
plt.subplot(212)if args.n_images is None:    plt.plot(mjd-mjd[0],x_right-x_right[ref_i],'k.')else:    plt.plot((mjd-mjd[0])[1:],(x_right-x_right[ref_i])[1:],'k.')
plt.title("x centroid at y~"+str(y_cont)+" (Comparison)")plt.xlabel("JD from start of observations")plt.ylabel("x centroid - reference")if args.no_save:    plt.show()else:    plt.savefig('xshift_y%d.png'%y_line_left)
    plt.show()
plt.figure(figsize=(11.69,8.27))
plt.subplot(211)if args.n_images is None:
    plt.plot(mjd-mjd[0],y_left-y_left[ref_i],'k.')else:    plt.plot((mjd-mjd[0])[1:],(y_left-y_left[ref_i])[1:],'k.')plt.title("y centroid of spectral line at y~"+str(y_line_left)+" (Target)")plt.xlabel("JD from start of observations")plt.ylabel("y centroid - reference")
plt.subplot(212)if args.n_images is None:
    plt.plot(mjd-mjd[0],y_right-y_right[ref_i],'k.')else:    plt.plot((mjd-mjd[0])[1:],(y_right-y_right[ref_i])[1:],'k.')plt.title("y centroid of spectral line at y~"+str(y_line_right)+" (Comparison)")plt.xlabel("JD from start of observations")plt.ylabel("y centroid - reference")if args.no_save:    plt.show()else:    plt.savefig('yshift_y%d.png'%y_line_left)    plt.show()
plt.figure()plt.scatter(range(len(max_counts)),max_counts)plt.xlabel("Frames from start")plt.ylabel("Max counts")if not args.no_save:    plt.savefig('max_counts.png')plt.show()
plt.figure()plt.scatter(range(len(fwhm_x)),fwhm_x,color='b')plt.xlabel("Frames from start")plt.ylabel("FWHM (arcsec)")if not args.no_save:    plt.savefig('fwhm_y%d.png'%y_cont)plt.show()if not args.no_save:    pickle.dump(fwhm_x,open('fwhm_y%d.pickle'%y_cont,'w'))    pickle.dump(x_right-x_right[ref_i],open('x2_yline_%d.pickle'%y_line_right,'w'))    pickle.dump(x_left-x_left[ref_i],open('x1_yline_%d.pickle'%y_line_left,'w'))    pickle.dump(y_right-y_right[ref_i],open('y2_yline_%d.pickle'%y_line_right,'w'))    pickle.dump(y_left-y_left[ref_i],open('y1_yline_%d.pickle'%y_line_left,'w'))    pickle.dump(max_counts,open('max_counts.pickle','w'))    