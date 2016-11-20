#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
#
# Copyright 2015
# Authors: Hector Vazquez Ramio (hvr@cefca.es)
#
# This file is part of JPAS pipeline jype
#
# jype is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# jype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with JPAS Pipeline. If not, see <http://www.gnu.org/licenses/>.
#
######################################################################
#
# Description:
#
#   Script to compute SNR for J-PAS filter system.
#
# Known issues:
#
# TODO:
#    * 
#
# Notes:
#
#
# Examples:
#
#
######################################################################

import sys
import os
import glob
import numpy as np
import pylab as plt
import datetime as dt
import tempfile

import pyfits
from scipy import ndimage
import ephem

import javadimm_utils as util
import dimm_conf 

try:
    import ipdb
    use_ipdb = True
except:
    use_ipdb = False
    
def go(conf):

    dimm = dimm_conf.DimmConf(
        conf['LAMBDA'], conf['DIAMETER'], conf['BASELINE'], 
        conf['ARCSEC_PIX'])
    
    #Make an observer
    oaj = ephem.Observer()
    #Location of OAJ
    oaj.long = str(conf['SITE_LONGITUDE']) # in string format
    oaj.lat = str(conf['SITE_LATITUDE'])  # in string format

    #Elevation of OAJ
    oaj.elev = conf['SITE_ELEVATION'] # (m)

    # horizon altitud (rad)
    hza = -np.arccos(ephem.earth_radius/(oaj.elev+ephem.earth_radius))


    # XXX Reading RoboDIMM data
    rdimm_file = os.path.join('./data/robodimm', 'Seeing161104.dat')
    pre = np.loadtxt(rdimm_file, unpack=True, dtype=str)
    rdimm = {}
    rdimm['DATE'] = np.array(pre[0])
    rdimm['UT'] = np.array(pre[1])
    rdimm['STAR'] = np.array(pre[2])
    rdimm['LONGITUDINAL'] = np.array(pre[3], float)
    rdimm['TRANSVERSE'] = np.array(pre[4], float)
    rdimm['AIRMASS'] = np.array(pre[5], float)
    
    for movie in conf['MOVIES']:

        # Root name for naming outputs
        ROOT_NAME = os.path.basename(movie).split(os.path.extsep)[0]
        
        # Outputs directory
        OUT_DIR = os.path.join(CONFIG['SELECTED_DIR'], ROOT_NAME+'files')
        if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)
        out_file = os.path.join(OUT_DIR, '%s.fits' % ROOT_NAME) 
        
        # Output seeing file
        seeing_file = out_file.replace(os.extsep+'fits', os.extsep+'dat')
        with open(seeing_file, 'w') as s_file:
            pass
            
        # Temporal directory where the frames will be extracted from movie
        frames_dir = os.path.join(OUT_DIR, 'frames')
        if not os.path.exists(frames_dir):
            os.makedirs(frames_dir)

        if conf['EXTRACT_FRAMES']:
            # Command used to extract frames from avi file
            if conf['USE_MPLAYER']:
            
                # Temporary copy of avi file in frames_dir
                file_descriptor, tmpavi = tempfile.mkstemp(
                    prefix=os.path.basename(movie).
                    split(os.extsep)[0]+'_', suffix=os.extsep+'avi',
                    dir=frames_dir, text=False)
    
                try:
                    
                    if os.path.isfile(tmpavi):
                        os.system('rm %s' % tmpavi)
                    cmd = 'ln -s %s %s' % (movie, tmpavi)
                    print cmd
                    res = os.system(cmd)

                    if res == 0:

                        cmd = "cd %s && " \
                            "mplayer -ao null -vo png %s && " \
                            "cd -" % (frames_dir, tmpavi)
                        print cmd
                        
                        res = os.system(cmd)
              
                        if res != 0:
                            print 'Frames extraction failed'
                            sys.exit()

                    else:
                        print \
                            'Something went wrong creating the ' \
                            'symbolic link %s to %s' % (movie, tmpavi)

                finally:

                    if os.path.isfile(tmpavi):
                        print 'Removing temporal file %s' % (tmpavi)
                        os.system('rm %s' % tmpavi)

            else:
                # Command used to extract frames from avi file
                cmd = \
                    "avconv -i %s %s" % (
                    movie, os.path.join(frames_dir, 'frame_%6d.png'))
            
                res = os.system(cmd)
            
                if res != 0:
                    print 'Frames extraction failed'
                    sys.exit()

            # Command to see avi metadata
            # avprobe ../T80_luckyimaging_20151015_focus-3_00000005.avi 
    
        # Length of 1arcsec in pixels
        one_arcsec_pix = 1./conf['ARCSEC_PIX']
    
        list_ = glob.glob(os.path.join(frames_dir, '*.png'))
        list_ = sorted(list_)
        num_frames = len(list_)

        if num_frames == 0:
            print 'No frames found %s' % os.path.join(
                frames_dir, 'frame_*.png')
            print 'Check paths or re-run with -x option to extract frames'
            sys.exit()
    
        num_analysis = num_frames
    
        # Checking the frames size
        dummy, xsize, ysize = util.read_img(
           list_[0], 'png', return_size=True)
    
        # Combined frames without re-centering
        combined_frames = np.zeros((conf['WINSIZE'], conf['WINSIZE']))
        # Combined frames subset without re-centering
        combined_frames_subset = {}
    
        Strehl_ratio = np.array([])
        Strehl_ratio_all = np.zeros(num_analysis)
        Strehl_ratio_sat = np.array([])
    
    
        # List of accepted frames
        list_ok = np.array([])
        ind_ok = np.array([])
        ind_sat = np.array([])
    
        # Maximum values of each frame
        max_values = np.zeros(num_analysis)
    
        new_hdu = pyfits.HDUList()

        EXPTIME = 8e-3
        
        not_sat = 0
        for k, file_ in enumerate(list_):
    
            file_ = list_[k]
            print file_
            
            img = util.read_img(file_, 'png')
            
            # Add to list_ok
            list_ok = np.append(list_ok, file_)
            ind_ok = np.append(ind_ok, k+1)
              
            # Appending extension
            new_hdu.append(pyfits.ImageHDU(np.array(img, np.short)))

            if (np.mod(k+1, conf['NFRAMES']) == 0) or (
                (k+1) == len(list_)):
#                Strehl_ratio = np.append(Strehl_ratio, sr_)
    
                num_ok = len(list_ok)
                num_sat = len(ind_sat)
    
                # Writing fits file
                cat_file = out_file.replace('.fits', '.cat')
                print 'Writing %s image...' % out_file
                new_hdu.writeto(out_file, clobber=True)
    
                # Run sextractor
                pre_cat = util.run_sextractor(
                    out_file, conf['CONFIG_SEX'], 
                    satur_level=conf['SATUR_LEVEL'],
                    detect_thresh=conf['DETECT_THRESH'])
                
                if isinstance(pre_cat, dict):   
                    # Checking the results (saturation, number of 
                    # detections, etc.)
                    dict_cat = util.reshape_catalog(pre_cat)
                    valid, rejected = util.check_catalog(dict_cat)
                
                    nvalid = np.sum(valid)
                    nrejected = np.sum(~valid)
                    nrej_ndet = np.sum(rejected['NDET'])
                    nrej_flags = np.sum(rejected['FLAGS'])
                else:
                    # No detections
                    nvalid = -1
                
                if nvalid < 5:
                    check = False
                else:
                    cat = util.flatten_valid_catalog(dict_cat, valid)

                    ndata = len(cat['X_IMAGE'])
                    XX = np.concatenate(
                        (cat['X_IMAGE'].reshape(1, ndata), 
                        cat['Y_IMAGE'].reshape(1, ndata))).T
        
                    from sklearn.cluster import KMeans
       
#                    kmeans = KMeans(n_clusters=4, n_init=10).fit(XX)
#                    YY = kmeans.predict(XX)

                    from sklearn.cluster import dbscan
                    YY = dbscan(XX, eps=5)[1]
                    
                    # Checking the result
                    check = (np.sum(YY == 0) == nvalid) and \
                        (np.sum(YY == 1) == nvalid) and \
                        (np.sum(YY == 2) == nvalid) and \
                        (np.sum(YY == 3) == nvalid)
  
                    if conf['PLOT_IT']:
                        for i in range(4):
                            mask = YY == i
                            plt.plot(
                                cat['X_IMAGE'][mask], 
                                cat['Y_IMAGE'][mask], 'o')
                        plt.grid(True)
                        plt.xlabel('X (pixels)')
                        plt.ylabel('Y (pixels)')
                        plt.savefig('./out/png/spots.png')
                        plt.close()
                
                if check:

                    nspot = {}
                    xmedians = np.array([
                        np.median(cat['X_IMAGE'][YY == 0]),
                        np.median(cat['X_IMAGE'][YY == 1]),
                        np.median(cat['X_IMAGE'][YY == 2]),
                        np.median(cat['X_IMAGE'][YY == 3])])
                    ymedians = np.array([
                        np.median(cat['Y_IMAGE'][YY == 0]),
                        np.median(cat['Y_IMAGE'][YY == 1]),
                        np.median(cat['Y_IMAGE'][YY == 2]),
                        np.median(cat['Y_IMAGE'][YY == 3])])
            

                    four = np.array(range(4))
                    nspot[0] = four[np.argmin(xmedians)]
                    remaining = (four != nspot[0])
    
                    nspot[1] = four[remaining][np.argmax(
                        xmedians[remaining])]
                
                    remaining = (four != nspot[0]) & (four != nspot[1])
                    nspot[2] = four[remaining][np.argmin(
                        ymedians[remaining])]
                
                    nspot[3] = four[
                        (four != nspot[0]) & (four != nspot[1]) & 
                        (four != nspot[2])][0]
    
    
                    if conf['PLOT_IT']:
                        for i in range(4):
                            mask = YY == nspot[i]
                            plt.plot(
                                cat['X_IMAGE'][mask], 
                                cat['Y_IMAGE'][mask], 'o',
                                label='%i' % i)
                        plt.legend(loc='best', numpoints=1)
                        plt.grid(True)
                        plt.xlabel('X (pixels)')
                        plt.ylabel('Y (pixels)')
                        plt.savefig('./out/png/spots_numbered.png')
                        plt.close()
    
    
                    spot = {}
                    for i in four:
                        spot[i] = YY == nspot[i]
        
                    # Find angles
                    angles = {}
                    angles[0] = np.arctan(
                        (cat['Y_IMAGE'][spot[1]]-
                        cat['Y_IMAGE'][spot[0]])/
                        (cat['X_IMAGE'][spot[1]]-
                        cat['X_IMAGE'][spot[0]]))
                    angles[1] = np.arctan(
                        (cat['Y_IMAGE'][spot[3]]-
                        cat['Y_IMAGE'][spot[2]])/
                        (cat['X_IMAGE'][spot[3]]-
                        cat['X_IMAGE'][spot[2]]))
                        
                    deltas = {}
                    for i in range(2):
                        deltas[i] = angles[i]-np.median(angles[i])
            
                    # Distances
                    lon_dist = {}
                    tra_dist = {}
            
                    lon_dist[0] = np.cos(deltas[0])*np.sqrt(
                        (cat['X_IMAGE'][spot[1]]-
                        cat['X_IMAGE'][spot[0]])**2+
                        (cat['Y_IMAGE'][spot[1]]-
                        cat['Y_IMAGE'][spot[0]])**2)
                    tra_dist[0] = np.sin(deltas[0])*np.sqrt(
                        (cat['X_IMAGE'][spot[1]]-
                        cat['X_IMAGE'][spot[0]])**2+
                        (cat['Y_IMAGE'][spot[1]]-
                        cat['Y_IMAGE'][spot[0]])**2)
            
                    lon_dist[1] = np.cos(deltas[1])*np.sqrt(
                        (cat['X_IMAGE'][spot[3]]-
                        cat['X_IMAGE'][spot[2]])**2+
                        (cat['Y_IMAGE'][spot[3]]-
                        cat['Y_IMAGE'][spot[2]])**2)
                    tra_dist[1] = np.sin(deltas[1])*np.sqrt(
                        (cat['X_IMAGE'][spot[3]]-
                        cat['X_IMAGE'][spot[2]])**2+
                        (cat['Y_IMAGE'][spot[3]]-
                        cat['Y_IMAGE'][spot[2]])**2)
                   
                    lon_pix = {}
                    tra_pix = {}
                    for i in range(2):
                        lon_pix[i] = np.var(lon_dist[i])
                        tra_pix[i] = np.var(tra_dist[i])
            
                    seeing_lon = {}
                    seeing_tra = {}
                    for i in range(2):
                        seeing_lon[i] = dimm.var2seeing(
                            np.var(dimm.pix2rad(lon_dist[i])[0]), 'L')[0]
                        seeing_tra[i] = dimm.var2seeing(
                            np.var(dimm.pix2rad(tra_dist[i])[0]), 'T')[0]
            
                    print seeing_lon
                    print seeing_tra
                    
                    seeings = np.array(
                        [seeing_lon[0], seeing_lon[1], 
                        seeing_tra[0], seeing_tra[1]])
                    print movie
                    print 'uncorrected seeing: %6.3f +/- %5.3f ' \
                        '(arcsec)' % (np.mean(seeings), 
                        np.std(seeings, ddof=1))
                    

                    # XXX Estimating time
                    # XXX Duration of the video
                    print movie
                    cmd = 'exiftool %s > /tmp/test.yaml' % movie
                    res = os.system(cmd)
                    movie_info = yaml.load(open('/tmp/test.yaml', 'r'))
                    str_duration = movie_info['Duration']
                    bit_depth = movie_info['Bit Depth']
#                    cmd = 'avprobe %s 2> /tmp/dummy.txt' % movie
#                    res = os.system(cmd)
#                    cmd = "cat /tmp/dummy.txt | grep Duration " \
#                        "| awk '{print $2}' > /tmp/dummy1.txt" 
#                    res = os.system(cmd)
#                    str_duration = str(np.loadtxt(
#                        '/tmp/dummy1.txt', dtype=str)).replace(',','')
                    dt_duration = dt.timedelta(
                        seconds=np.sum(np.array(
                        str_duration.split(':'), float)*
                        [3600., 60., 1.]))
                        
                    # Frames per second
                    fps = float(len(list_))/dt_duration.total_seconds()
                    print movie, str_duration, dt_duration, fps
                    
                    mask_movie = \
                        conf['METADATA']['FILENAME'] == \
                        os.path.basename(movie)
                    end_dt = conf['METADATA']['END_DT'][mask_movie][0]
                    ini_dt = end_dt-dt_duration
                    
                    # Time elapsed since the beginning of the video
                    dt_elapsed = dt.timedelta(
                        seconds=(k-float(conf['NFRAMES'])/2.)/
                        float(fps))
                    
                    # Estimated mid-observation time of the current
                    # measurement UT    
                    dt_estimated = ini_dt+dt_elapsed
                    
                    print 't0=%s est_time=%s' % (
                        ini_dt.isoformat(), dt_estimated.isoformat())
                    
                    
                    # date "_dt" in ephem format
                    ephem_date = ephem.Date(dt_estimated) 
                    oaj.date = ephem_date
                    
                    obs_star = conf['METADATA']['STAR'][mask_movie][0]
                    mask_star = conf['STAR']['NAME'] == obs_star
                    
                    star = ephem.FixedBody()
                    star._ra = conf['STAR']['RA'][mask_star][0]
                    star._dec = conf['STAR']['DEC'][mask_star][0]
                    
                    star.compute(oaj)
                    star_alt = star.alt
                    airmass = 1./np.cos(np.pi/2.-star.alt)
                    
                    print '%s alt=%s  airmass=%6.4f' % (
                        obs_star, star_alt, airmass)
                    
        
#                    # XXX Taking RoboDIMM data
#                    sUT = os.path.basename(movie).split('UT')[0]
#                    t0 = float(sUT[:2])+float(sUT[2:4])/60.+ \
#                        float(sUT[4:])/3600.

#                    est_time = t0+k*(EXPTIME+1./float(FPS))/3600.
#                    print 't0=%s est_time=%s' % (t0, est_time)
#                    
#                    est_UT = '%02i%02i%02i' % (
#                        int(est_time),
#                        int((est_time-int(est_time))*60.), 
#                        int(est_time*3600.-int(est_time)*3600.-
#                        int((est_time-int(est_time))*60.)*60.))
#                    print est_UT   
#             
#                    t1s = np.array([
#                        float(time[:2])+float(time[2:4])/60.+
#                        float(time[4:])/3600 
#                        for time in rdimm['UT']])
#                    
#                    airmass = np.interp(est_time, t1s, rdimm['AIRMASS'])
#                    print 'airmass = %6.4f' % airmass
            
#                    if 'airmass' not in vars().keys():
#                        ans = raw_input('Aprox. airmass: ')
#                        airmass = float(ans)    
                    
                    corr_seeings = dimm.seeing_corr_airmass(
                        seeings, airmass)[0]
                    seeing = np.mean(corr_seeings)
                    seeing_err = np.std(corr_seeings, ddof=1) 
                    print '  corrected seeing: %6.3f +/- %5.3f ' \
                        '(arcsec)' % (seeing, seeing_err)

                    with open(seeing_file, 'a') as s_file:
                        s_file.write(
                            '%s %6i %6.3f %5.3f %6.4f %10s %5.2f '
                            '%4i %4i %4i %4i '
                            '%6.3f %6.3f %6.3f %6.3f %6.1f %6.1f '
                            '%5.1f %5.1f %5.1f %5.1f\n' % (
                            util.round2sec(dt_estimated),
                            k+1, 
                            seeing,
                            seeing_err,
                            airmass,
                            obs_star,
                            np.rad2deg(star.alt.real),
                            nvalid,
                            nrejected,
                            nrej_ndet,
                            nrej_flags,
                            corr_seeings[0],
                            corr_seeings[1],
                            corr_seeings[2],
                            corr_seeings[3],
                            np.median(np.rad2deg(angles[0])),
                            np.median(np.rad2deg(angles[1])),
                            np.median(lon_dist[0]),
                            np.median(tra_dist[0]),
                            np.median(lon_dist[1]),
                            np.median(tra_dist[1])
                            ))

                
                
#                if use_ipdb:
#                    ipdb.set_trace()

                new_hdu = pyfits.HDUList()

        if 'airmass' in vars().keys():
            del airmass

    
if __name__ == '__main__':

    import argparse
    import yaml
    
    parser = argparse.ArgumentParser()

#    parser.add_argument(
#        'input', type=str, nargs='*', help='Input movie/s filename/s')

    parser.add_argument(
        'conf', type=str, default='javadimm.yaml',
        help='yaml configuration file for the instrument'
        'file will be read.')

    parser.add_argument(
        '-x', action='store_true', default=False,
        help='Extract frames from input movie/s.')

    parser.add_argument(
        '-p', '--plot', default=False, action='store_true',
        help='Show plots on screen')
        
    parser.add_argument(
        '--mplayer', action='store_true', default=False,
        help='Use mplayer to extract frames from input movie/s.')

    parser.add_argument(
        '-S', '--Strehl', type=float, default=0.95,
        help='Fraction of Strehl ratio to define threshold for lucky ' 
        ' imaging (only applies if BEST_FRACTION == -1)')

#    parser.add_argument('--sat', type=float, default=253,
#                        help="Satuartion level (ADUs)")

    parser.add_argument(
        '-w', '--wsize', type=int, default=24,
        help='Size of the square subframe window side (pixels)')

                        
    args = parser.parse_args()
    
    CONFIG = yaml.load(open(args.conf, 'r'))
    
    # yaml does not understand scientific notation (take is as a str)
    CONFIG['LAMBDA'] = float(CONFIG['LAMBDA'])
    
    # Directory selection
    list_dir = glob.glob(os.path.join(CONFIG['INPUT_DIR'], '*'))
    list_dir = sorted(list_dir)
    for i, dir_ in zip(range(len(list_dir)), list_dir):
        print '%3i: %s' % (i+1, dir_)

    dir_num = int(raw_input('Select directory number: '))

    CONFIG['SELECTED_DIR'] = list_dir[dir_num-1]

    # Movie/s selection
    list_avi = glob.glob(os.path.join(CONFIG['SELECTED_DIR'], '*.avi'))
    list_avi = sorted(list_avi)
    for i, avi_ in zip(range(len(list_avi)), list_avi):
        print '%3i: %s' % (i+1, avi_)

    pre_avi_num = raw_input(
        'Select avi/s numbers (e.g. "0" or "0 4 5" or "2-5"): ')

    if '-' in pre_avi_num:
        pre_ = pre_avi_num.split('-')
        avi_num = np.arange(int(pre_[0]), int(pre_[1])+1)
    else:
        avi_num = pre_avi_num.split(' ')
    CONFIG['MOVIES'] = [list_avi[int(x)-1] for x in avi_num]

    for files in CONFIG['MOVIES']:
        print files

    # XXX Read info from info.txt
    data_info = np.loadtxt(
        os.path.join(CONFIG['SELECTED_DIR'], 'info.txt'), dtype=str).T
    CONFIG['METADATA'] = {}
    CONFIG['METADATA']['END_DATE'] = data_info[0]
    CONFIG['METADATA']['END_UT'] = data_info[1]
    CONFIG['METADATA']['FILENAME'] = data_info[3]
    CONFIG['METADATA']['END_DT'] = np.array([
        dt.datetime.strptime(date_+'T'+ut_[:-3], 
        '%Y-%m-%dT%H:%M:%S.%f')
        for date_, ut_ in zip(
        CONFIG['METADATA']['END_DATE'], CONFIG['METADATA']['END_UT'])])

    # XXX Read stars
    data_stars = np.loadtxt('./data/stars_list.dat', delimiter=';',
        dtype=str).T
    CONFIG['STAR'] = {}
    CONFIG['STAR']['NAME'] = data_stars[0]
    CONFIG['STAR']['CONST'] = data_stars[1]
    CONFIG['STAR']['RA'] = np.array(data_stars[2])
    CONFIG['STAR']['DEC'] = np.array(data_stars[3])
    CONFIG['STAR']['V'] = np.array(data_stars[4], float)

    CONFIG['METADATA']['STAR'] = np.array([], str)
    for name_ in CONFIG['METADATA']['FILENAME']:
        for star_ in CONFIG['STAR']['NAME']:
            # XXX The [:-1] is for a Capell_16bits... name
            if star_.upper()[:-1] in name_.upper():
                CONFIG['METADATA']['STAR'] = np.append(
                    CONFIG['METADATA']['STAR'], star_.upper())
    
    
    # Extract frames or not
    CONFIG['EXTRACT_FRAMES'] = args.x
    
#    # Input movie/s filenames
#    CONFIG['MOVIES'] = args.input

    # Plot results on screen
    CONFIG['PLOT_IT'] = args.plot

    # Plot results on screen
    CONFIG['USE_MPLAYER'] = args.mplayer

    # Fraction of Strehl ratio to define threshold for luccky imaging 
    # (Only applies if BEST_FRACTION == -1)
#    STREHL_FACTOR = 0.95
    CONFIG['STREHL_FACTOR'] = args.Strehl
    
    # Saturation level (counts)
#    CONFIG['SATUR_LEVEL'] = args.sat
    
    # Window size of the cropped image (WARNING: even value required)
#    WINSIZE = 44
    CONFIG['WINSIZE'] = args.wsize
    
    go(CONFIG)    
    
    
