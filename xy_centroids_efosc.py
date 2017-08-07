import numpy as np

file_names = sorted(glob.glob("EFOSC_Spectrum*.fits"))


def extract_flux_x(frame,x_beg,x_end,y_cont,dy_cont,spec_w,bkg_w,bkg_off):
    
    frame_subset = frame[y_cont-dy_cont//2:y_cont+dy_cont//2+1,x_beg:x_end]
    max_ind = np.median(np.argmax(frame_subset,axis=1)).astype(int)
    spec_2d = frame_subset[:,max_ind-spec_w//2:max_ind+spec_w//2+1]
    
    bkg_1_2d = frame_subset[:,max_ind-bkg_off-bkg_w//2:max_ind-bkg_off+bkg_w//2+1]
    bkg = np.column_stack((bkg_1_2d,bkg_2_2d)).mean(axis=1)
    spec = (spec_2d - bkg[:,None]).sum(axis=0)
    return x_beg+max_ind, spec

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
  
    # Cubic fit to the 
def moffat_x_mod(a,r):
def moffat_x(a,r,y):
def moffat_y_mod(a,r):
def moffat_y(a,r,y):
max_counts = []
for i,f in enumerate(l):
    hdu = fits.open(f)
    if hdr["OBJECT"] == obj_name:
        frame1 = hdu[0].data
        mjd.append(hdr["MJD-OBS"])
        hdu.close()
        # Extract flux in x collapsed in y
        x = peak_x+np.arange(-len(flux_x)/2,len(flux_x)/2)
        # Model flux in x with a gaussian.
        parms_x_left = sopt.leastsq(moffat_x,np.array([6.,5.,peak_x,200000.,0.]),args=(x,flux_x))
        flux_y = extract_flux_y(frame1,x_beg_left,x_end_left,y_line_left,dy_line,spec_w,bkg_w,bkg_off)
        y = y_line_left - dy_line//2 + np.arange(len(flux_y)) 
        parms_y_left = sopt.leastsq(moffat_y,np.array([4.0,1.5,y_line_left,-150000.,0.,50000.]),args=(y,flux_y))  
        if f == ref_frame:
            xfine = np.arange(x[0],x[-1],0.001)
            plt.figure()
            yfine = np.arange(y[0],y[-1],0.001)
            plt.figure()
        if nwindows == 2:

        y = y_line_right - dy_line//2 + np.arange(len(flux_y)) 
        parms_y_right = sopt.leastsq(moffat_y,np.array([4.,1.0,y_line_right,-5000000.,0.,200000.]),args=(y,flux_y)) 
        x_left.append(parms_x_left[0][2])
        y_left.append(parms_y_left[0][2])
        x_right.append(parms_x_right[0][2])
        y_right.append(parms_y_right[0][2])
        if f == ref_frame:
            xfine = np.arange(x[0],x[-1],0.001)
            plt.figure()
            yfine = np.arange(y[0],y[-1],0.001)
            plt.figure()
    else:
        hdu.close()
x_left = np.array(x_left)
ref_i = np.where(np.array(flist)==ref_frame)
print "X_left now = ",x_left[-1],"; X_left ref = ",x_left[ref_i]
plt.figure(figsize=(11.69,8.27))
plt.subplot(211)
    plt.plot(mjd-mjd[0],x_left-x_left[ref_i],'k.')
plt.subplot(212)
plt.title("x centroid at y~"+str(y_cont)+" (Comparison)")
    plt.show()

plt.subplot(211)
    plt.plot(mjd-mjd[0],y_left-y_left[ref_i],'k.')
plt.subplot(212)
    plt.plot(mjd-mjd[0],y_right-y_right[ref_i],'k.')

