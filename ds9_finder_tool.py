import os
from astroquery.vizier import Vizier
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
import argparse

parser = argparse.ArgumentParser(description='Make finding chart with DS9')
parser.add_argument("target",help="Name of object to be resolved by Vizier",type=str)
parser.add_argument("-inst","--instrument",help="Which instrument? Needed for field of view & slit width and length. Default = ACAM",default='ACAM')
parser.add_argument('-cat',"--catalog",help="Name of catalog to query [APASS9/UCAC4], default = UCAC4",default='UCAC4')

parser.add_argument("-s","--slit",help="Use this to define which ACAM slit is used, default = 40",type=int,default=40)
parser.add_argument("-v","--magnitude_cut",help="If wanting to cut V magnitude at a partciular value. Default = 14.",type=float,default=14)

parser.add_argument("-2mass","--twomass",help='If wanting to display a 2MASS bkg image.',action='store_true')
parser.add_argument("-apass","--apass",help="If wanting to use APASS9, rather than UCAC4 catalog",action="store_true")

parser.add_argument("--qr",help="Use this argument if the database cannot resolve the object's name. Must be used like --qr '00h00m00s 00d00m00s'",nargs='+',type=str)

parser.add_argument("--tc",help="Use if wanting to manually define target coords. Must be used like --tc '00h00m00s 00d00m00s'",nargs='+',type=str)
parser.add_argument("--cc",help="Use if wanting to manually define comparison coords. Must be used like --cc '00h00m00s 00d00m00s'",nargs='+',type=str)
parser.add_argument("--mc",help="Use if wanting to manually define midpoint coords. Must be used like --mc '00h00m00s 00d00m00s'",nargs='+',type=str)
parser.add_argument("-pa","--position_angle",help="Use if wanting to manually define position angle",type=float)

parser.add_argument("-p2s","--print_to_screen",help='Use this function to print the DS9 commands to the terminal for copying and pasting if having issues with automatic opening',action='store_true')

args = parser.parse_args()

if args.instrument == 'ACAM':
    fov = 8.1 # Field of view in arcmin
    search_radius = 6.8
    slit_width = args.slit
    
    if slit_width == 40:
        slit_length = 7.6*60
        
    elif slit_width == 27:
        slit_length = 6.8*60
    
    else:
        raise ValueError('Slit must be either 40 or 27')
        
if args.instrument == 'EFOSC':
    fov = 4.1
    search_radius = 4.1
    slit_width = 15
    slit_length = 4.1*60
    
    
if args.tc:
    targ_ra,targ_dec = args.tc[0].split()
    comp_ra,comp_dec = args.cc[0].split()
    mid_ra,mid_dec = args.mc[0].split()
    pa = args.position_angle

if args.qr:
    targ_ra_hms,targ_dec_dms = args.qr[0].split()
    c = SkyCoord(targ_ra_hms,targ_dec_dms)
    targ_ra = c.ra.deg
    targ_dec = c.dec.deg
    
v = Vizier(columns=['Full','+_r','_RAJ2000', '_DEJ2000','B-V', 'Vmag','Bmag'],column_filters={"Vmag":"<%f"%args.magnitude_cut})

if args.apass:
    catalog = 'APASS9'
else:
    catalog = 'UCAC4'

if not args.tc and not args.qr:

    Vizier.ROW_LIMIT = -1

    if args.apass:
        target = v.query_object(args.target,catalog=['APASS9'])
        viz_result = target[u'II/336/apass9']
    
    else:
        target = v.query_object(args.target,catalog=['UCAC4'])
        try:
            viz_result = target[u'I/322A/out']
        except:
            raise KeyError('Problem in querying UCAC4, try APASS9 instead')
    
    target_dict = {}
    
    if len(viz_result) != 1:
        print viz_result
        selection = int(raw_input("More than one candidate found, which is the correct object?  "))
        viz_result = viz_result[selection]
        if args.apass:
            keys = ['_r','_RAJ2000', '_DEJ2000','B-V', 'Vmag','Bmag']
        else:
            keys = ['_r','_RAJ2000', '_DEJ2000','Vmag','Bmag']
        for key in keys:
            target_dict[key] = np.array(viz_result[key])
    
    else:
        for key in viz_result.keys():
            target_dict[key] = np.array(viz_result[key].filled(999).tolist())
    
    targ_ra = target_dict['_RAJ2000']
    targ_dec = target_dict['_DEJ2000']
    
if not args.tc:  
    search_radius = str(search_radius)+'m'
    
    comparisons = v.query_region(SkyCoord(targ_ra,targ_dec,unit=(u.deg,u.deg),frame='icrs'),radius=search_radius,catalog=[catalog]) # Query UCAC4 for comparison objects within search radius
    
    if args.apass:
        table =  comparisons[u'II/336/apass9']
    else:
        table =  comparisons[u'I/322A/out']
        
    table.sort('_r')
    print table
    print "==="
    which_comp = int(raw_input("Which comparison to use?     ")) # Need to make a selection for which comparison object
    comp_ra = table[which_comp]['_RAJ2000']
    comp_dec = table[which_comp]['_DEJ2000']
    
    # find rough midpoint between target and comparison to load image into DS9
    mid_ra = (targ_ra + comp_ra)/2.
    mid_dec = (targ_dec + comp_dec)/2.
    mid = SkyCoord(ra=mid_ra*u.degree,dec=mid_dec*u.degree)
    
    # CALCULATE POSITION ANGLE USING FOUR PART FORMULA SEE Smart 1949 textbook on spherical astronomy, p12
    
    #np.cos(a)*np.cos(C) = np.sin(a)*np.cot(b) - np.sin(C)*np.cot(B)
    
    targ_coords = SkyCoord(ra=targ_ra*u.degree,dec=targ_dec*u.degree)
    comp_coords = SkyCoord(ra=comp_ra*u.degree,dec=comp_dec*u.degree)
    
    rarad1 = np.array([targ_ra*np.pi/12.0])
    rarad2 = np.array([comp_ra*np.pi/12.0])
    dcrad1 = np.array([targ_dec*np.pi/180.0])
    dcrad2 = np.array([comp_dec*np.pi/180.0])
    
    radif  = rarad2-rarad1
    angle  = np.arctan(np.sin(radif),np.cos(dcrad1)*np.tan(dcrad2)-np.sin(dcrad1)*np.cos(radif))
    
    sep = 60*(np.arccos(np.sin(targ_dec*np.pi/180.0)*np.sin(comp_dec*np.pi/180) + np.cos(comp_dec*np.pi/180)*np.cos(targ_dec*np.pi/180)*np.cos(comp_ra*np.pi/180 - targ_ra*np.pi/180)))*180.0/np.pi
    
    c1 = SkyCoord(targ_ra*u.deg,targ_dec*u.deg)
    c2 = SkyCoord(comp_ra*u.deg,comp_dec*u.deg)
    
    pa = c1.position_angle(c2).degree

    print "==="
    print "Target coords = ",targ_ra,targ_dec
    print "Target coords = ",targ_coords.to_string('hmsdms')
    print "Comparison coords = ",comp_ra,comp_dec
    print "Comparison coords = ",comp_coords.to_string('hmsdms')
    print "Position angle = ",pa
    print "Separation = ",sep," arcmin"
    print "Midpoint coords = ",mid.to_string('hmsdms')
    print "Comparison Vmag = ",table[which_comp]['Vmag']
    print "Comparison B-V = ",table[which_comp]['Bmag'] - table[which_comp]['Vmag']
    print "==="

print 'PLEASE CHECK BOX / SLIT IN DS9, IT MAY BE OF THE INCORRECT DIMENSIONS!!!!!'
print 'TARGET IS RED CROSSHAIR, COMPARISON IS BLUE'


if args.twomass:
    image = '-2mass'
    slit_correction_factor = 1
else:
    image = '-dsssao'
    slit_correction_factor = 0.5887445887445888 # This is needed as the scale of the dssao image is slightly off, not sure why.


MID_RA = mid_ra
MID_DEC = mid_dec

TARG_RA = targ_ra
TARG_DEC = targ_dec

COMP_RA = comp_ra
COMP_DEC = comp_dec

if not args.print_to_screen:
    os.system('ds9 %s size %.1f %.1f arcmin %s coord %f %f degrees  \
                      -regions system wcs -regions command "point %fd %fd # point=cross 40 color=red" \
                      -regions command "point %fd %fd # point=cross 40 color=blue" \
                      -regions command "box(%fd,%fd,%f",%f",%f)" \
                      -grid yes -grid grid no -grid type publication -grid axes type exterior \
                      -grid title no -grid numerics color black -grid labels yes -grid labels color black \
                       -grid numerics type exterior \
                       -grid title yes -grid title text %s -grid title def no \
                        -wcs fk5 &' %(image,fov,fov,image,MID_RA,MID_DEC,TARG_RA,TARG_DEC,COMP_RA,COMP_DEC,MID_RA,MID_DEC,slit_length*slit_correction_factor,slit_width*slit_correction_factor,pa-90,args.target))

else:
    print '\n ------- COPY AND PASTE THE FOLLOWING INTO A TERMINAL ---------'
    print 'ds9 %s size %.1f %.1f arcmin %s coord %f %f degrees  \
                      -regions system wcs -regions command "point %fd %fd # point=cross 40 color=red" \
                      -regions command "point %fd %fd # point=cross 40 color=blue" \
                      -regions command "box(%fd,%fd,%f",%f",%f)" \
                      -grid yes -grid grid no -grid type publication -grid axes type exterior \
                      -grid title no -grid numerics color black -grid labels yes -grid labels color black \
                       -grid numerics type exterior \
                       -grid title yes -grid title text %s -grid title def no \
                        -wcs fk5 &' %(image,fov,fov,image,MID_RA,MID_DEC,TARG_RA,TARG_DEC,COMP_RA,COMP_DEC,MID_RA,MID_DEC,slit_length*slit_correction_factor,slit_width*slit_correction_factor,pa-90,args.target)
print '\n'
