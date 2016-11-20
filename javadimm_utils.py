import os
import numpy as np
from PIL import Image

def read_img(sfile, sformat, return_size=False):
    r''' Reads a bmp or png image.
    '''
    
    if sformat in ('bmp', 'png'):
        im = Image.open(sfile)
        im1 = im.resize((im.size[0], im.size[1]))

    img = np.array(im1, float)
    if len(img.shape) == 3:
        img = img[:, :, 0]
        
    if return_size:
        return img, im.size[0], im.size[1]
    else:
        return img


def eval_strehl(image, lambda_, diameter, arcsec_pix):
    r''' Computes Strehl ratio:
    
       lambda_: effective wavelength (m)
       diameter: aperture (m)
       arcsec_pix: size of the pixel in arcsec
    '''
    
    pre = lambda_/(diameter*np.deg2rad(arcsec_pix/3600.))
    background = np.median(image)
    max_flux = np.max(image-background)
    total_flux = np.sum(image-background)
    
    return (max_flux/total_flux)*(4./np.pi)*pre*pre
    
def read_sex_cat(catalog):
    r''' Reads sextractor catalag.
    '''
    
    # Dictionary definition
    output = {}
    
    with open(catalog, 'r') as fcat:
        count = 0
        keys = []
        for line in fcat.readlines():
            if line != '\n':
               if (line[0] == '#'):
                   keys.append(line.replace('  ',' ').split()[2]) 
                   output[keys[-1]] = np.array([])
               else:
                   pre = line.replace('  ',' ').split()
                   count += 1
                   for key, value in zip(keys, pre):
#                       print catalog, key, value
                       if (key == 'NUMBER') or (key == 'FLAGS'):
                           output[key] = np.append(output[key], int(value))
                       else:
                           output[key] = np.append(output[key], np.float(value))
    if count == 0:
       return 1
    else:
       return output

def read_and_crop(filename, window_size, xsize, ysize,
        xmax=None, ymax=None, return_max=False):
    r''' Reads and crops the image acording to window_size. If
    xmax and ymax are provided, then the image is cropped using
    them as the coordinates of the center. Otherwise the position
    of the maximum is computed and the image cropped accordantly.
    '''

    frame_ = read_img(filename, 'png')
    
    if (xmax is None) or (ymax is None):
        # Maximum search
        argmax = np.argmax(frame_)
        xmax = np.divide(argmax, xsize)
        ymax = np.mod(argmax, ysize)
    
    window = frame_[
        xmax-window_size/2:xmax+window_size/2,
        ymax-window_size/2:ymax+window_size/2]
    
    if return_max:
        return window, xmax, ymax
    else:
        return window


def run_sextractor(fits_file, config_file, satur_level=251, 
    detect_thresh=3.):
    r''' Runs sextractor and returns the catalogue in a dictionary.
    '''

    # Temporal file
    out_cat = '/dev/shm/tmp.cat'
    
#    command = \
#        'sex -c %s %s -DETECT_THRESHOLD %5.1f -ANALYSIS_THRESH 1000 -CATALOG_NAME %s' % \
#        (config_file, fits_file, detect_thresh, out_cat)

    command = \
        'sex {fits_file} -DETECT_THRESH {det_thr} -FILTER N ' \
        '-SATUR_LEVEL {satur_level} -CATALOG_NAME {cat_name}'.format(
        fits_file=fits_file, det_thr=detect_thresh, 
        satur_level=satur_level, cat_name=out_cat)
    
    print command
    
    result = os.system(command)

    if result == 0:
        # Reads the resulting catalogue
        dict_cat = read_sex_cat(out_cat)
        return dict_cat
    else:
        print 'Sextractor failed'
        return result


def reshape_catalog(catalog):
    r''' Reshapes a catalag obtained from a multi-extension fits to 
    a dictionary with an element per extension.
    '''
    
    out = {}
    keys = catalog.keys()
    
    argsw = np.where(catalog['NUMBER'] == 1)[0]
    argsw = np.append(argsw, len(catalog['NUMBER']))

    for k, (i, j) in enumerate(zip(argsw[:-1], argsw[1:])):
        out[k] = {}
        for key in keys:
            out[k][key] = catalog[key][i:j]
    
    return out

def flatten_valid_catalog(catalog, valid):
    r''' Returns a flattened dict catalog with the valid elements only.
    '''
    
    out = {}
    nkeys = catalog.keys()
    keys = catalog[0].keys()
    
    for k, bol in zip(nkeys, valid):
        if bol:
            for key in keys:
                if key not in out.keys():
                    out[key] = np.array([])
                out[key] = np.append(out[key], catalog[k][key])
    
    return out
    
def check_catalog(catalog, nspots=4):
    r''' Check the catalog: saturation, number of detections, etc.
    '''
    
    valid = np.array([], bool)
    rejected = {}
    rejected['FLAGS'] = np.array([], bool)
    rejected['NDET']  = np.array([], bool)
    
    for k in catalog.keys():

        cat = catalog[k]

        rej = False
        
        if np.max(cat['NUMBER']) != nspots:
            rejected['NDET'] = np.append(
                rejected['NDET'], True)
        else:
            rejected['NDET'] = np.append(
                rejected['NDET'], False)
            
        if np.sum(cat['FLAGS'] != 0) > 0:
            rejected['FLAGS'] = np.append(
                rejected['FLAGS'], True)
        else:
            rejected['FLAGS'] = np.append(
                rejected['FLAGS'], False)
        
        valid = np.append(
            valid, 
            np.logical_not(rejected['NDET'][-1] or 
            rejected['FLAGS'][-1]))
        
    return valid, rejected
        
        
def round2sec(dt_):
    ''' Returns the date (str) with the time rounded according to the 
    number of seconds.
    '''
    
    import datetime as dt
    
    if dt_.microsecond >= 500000:
        # Add 1 second
        dt_ += dt.timedelta(0, 1)
    
    return dt_.isoformat().split('.')[0]


