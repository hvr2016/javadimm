import numpy as np
import pdb

class DimmConf(object):

    def __init__(self, lambda_, diameter, baseline, pixscale):
	
	self.lambda_ = lambda_
	self.diameter = diameter
	self.baseline = baseline
	self.b = self.baseline/self.diameter
        self.pixscale = pixscale
    
        self.kl, self.kt = self.klt()
        self.R = self.kl/self.kt
        
        
    def klt(self):
        ''' K_l and Kl coeffitients calculation for the given 
            configuration.
        '''
        b_13 = np.power(self.b, -1./3.)
        b_73 = np.power(self.b, -7./3.)
        
        kl = 0.340*(1.-0.570*b_13-0.040*b_73)
        kt = 0.340*(1.-0.855*b_13+0.030*b_73)
        
        return kl, kt
        
    def arcsec2pix(self, value):
        ''' Transforms arcsec to pixels '''
        
        if np.isscalar(value):
            out = np.float(value)/self.pixscale
        else:
            out = np.array(value, float)/self.pixscale
        

        return (out, "pixels")

    def arcsec2rad(self, value):
        ''' Transforms arcsec to radians. '''
        
        if np.isscalar(value):
            out = np.deg2rad(np.float(value)/3600.)
        else:
            out = np.deg2rad(np.array(value, float)/3600.)
        
        return (out, "rad")
        

    def rad2arcsec(self, value):
        ''' Transforms radians to arcsec. '''
        
        return (np.rad2deg(value)*3600., "arcsec")


    def rad2pix(self, value):
        ''' Transforms radians to pixels '''
        
        return (self.rad2arcsec(value)[0]/self.pixscale, "pixels")


    def pix2arcsec(self, value):
        ''' Transforms pixels to arcsec. '''
        
        if np.isscalar(value):
            out = np.float(value)*self.pixscale
        else:
            out = np.array(value, float)*self.pixscale

        return (out, "arcsec")


    def pix2rad(self, value):
        ''' Transforms pixels to rad. '''
        
        if np.isscalar(value):
            out = self.arcsec2rad(np.float(value)*self.pixscale)[0]
        else:
            out = self.arcsec2rad(np.array(value, float)* \
                  self.pixscale)[0]


        return (out, "rad")

    def seeing_corr_airmass(self, seeing, airmass):
        ''' Corrects the measured seeing in arcsec using the airmass.
        '''
        
        corr_seeing = seeing*np.power(airmass, -3./5.)
            
        return (corr_seeing, "arcsec")


    def seeing_uncorr_airmass(self, seeing, airmass):
        ''' Uncorrects the already corrected seeing in arcsec using the 
        airmass.
        '''
        
        uncorr_seeing = seeing*np.power(airmass, 3./5.)
            
        return (uncorr_seeing, "arcsec")

    def var2seeing(self, var_rad, slt):
        ''' Transforms variance in radians square to seeing in
        arcsec.
        '''
        
        if slt.upper() == 'L':
            klt = self.kl
        elif slt.upper() == 'T':
            klt = self.kt
        else:
            sys.exit('Longitudinal or transverse component not defined')

        seeing_rad = 0.98*np.power(var_rad/klt, 3./5.)* \
                          np.power(self.lambda_/self.diameter, 
                                   -1./5.)
            
        return self.rad2arcsec(seeing_rad)
                    
    def seeing2var(self, seeing, slt):
        ''' Transforms seeing in arcsec to variance in radians square.
        '''
        
        if slt.upper() == 'L':
            klt = self.kl
        elif slt.upper() == 'T':
            klt = self.kt
        else:
            sys.exit('Longitudinal or transverse component not defined')

        seeing_rad = self.arcsec2rad(seeing)[0]
        var_rad = klt*np.power(self.lambda_/self.diameter, 2.)*\
                      np.power(self.diameter*seeing_rad/ \
                               (0.98*self.lambda_), 5./3.)

        return (var_rad, "rad*rad")   
                    
    def seeing2std(self, seeing, slt):
        ''' Transforms seeing in arcsec to standard deviation 
        - sqrt(variance) - in radians.
        '''
        
        var_rad = self.seeing2var(seeing, slt)[0]

        return (np.sqrt(var_rad), "rad")


    def std2seeing(self, std_rad, slt):
        ''' Transforms standard deviation - sqrt(variance) - in 
        radians to seeing in arcsec.
        '''
        
        var_rad = std_rad*std_rad
        
        return self.var2seeing(var_rad, slt)


    def r02seeing(self, r0):
        ''' Transforms Fried parameter r0 in cm to seeing in arcsec. 
        '''
        
        r0_m = np.float(r0)*1e-2
        seeing_rad = 0.98*self.lambda_/r0_m

        return self.rad2arcsec(seeing_rad)
           
    def seeing2r0(self, seeing):
        ''' Transforms seeing in arcsec to Fried parameter r0 in cm.
        '''
        
        seeing_rad = self.arcsec2rad(seeing)[0]
        r0 = 0.98*self.lambda_/seeing_rad

        return (r0*100, "cm")
    

    def Airy(self):
        ''' Returns the system's Airy disk diameter in arcsec.
        Lambda and diameter in the same units.
        '''
        disk_rad = 1.22*self.lambda_/self.diameter
        
        return self.rad2arcsec(disk_rad)

    def seeing_noise_corr(self, fwhm_l, fwhm_t):
        ''' Estimates the var_noise from the ratio fwhm_l/fwhm_t and 
        corrects fwhm_l and fwhm_t for var_noise. fwhm_l and fwhm_t in 
        arcsec.
        ''' 
        # From seeing to variance
        var_l = self.seeing2var(fwhm_l, 'l')
        var_t = self.seeing2var(fwhm_t, 't')

        # The observed ratio
        ratio_obs = var_l[0]/var_t[0]

        # Theoretical ratio
        ratio_theo = self.R

        var_n_var_t = (ratio_theo-ratio_obs)/(ratio_obs-1)
        
        # Inferred variance caused by noise
        var_n = var_t[0]*(var_n_var_t/(1+var_n_var_t))

        # Correction of var_l and var_t for var_n
        corr_var_l = var_l[0]-var_n
        corr_var_t = var_t[0]-var_n

        # From variance to seeing: corrected fwhm_l and fwhm_t
        corr_fwhm_l = self.var2seeing(corr_var_l, 'l')
        corr_fwhm_t = self.var2seeing(corr_var_t, 't')

        return corr_fwhm_l, corr_fwhm_t
        
    def seeing_noise_corr_from_ratio(self, fwhm_l, fwhm_t, ratio_obs):
        ''' Estimates the var_noise from a given ratio of variances
        ratio_obs and corrects fwhm_l and fwhm_t for var_noise. fwhm_l 
        and fwhm_t in arcsec.
        ''' 
        # From seeing to variance
        var_l = self.seeing2var(fwhm_l, 'l')
        var_t = self.seeing2var(fwhm_t, 't')

        # Theoretical ratio
        ratio_theo = self.R

        var_n_var_t = (ratio_theo-ratio_obs)/(ratio_obs-1)
        
        # Inferred variance caused by noise
        var_n = var_t[0]*(var_n_var_t/(1+var_n_var_t))

        # TODO Improve this
        # The correction will be done only if the following condition
        # is satisfied
        ind_corr = ratio_obs < ratio_theo
        corr_var_l = var_l[0]
        corr_var_t = var_t[0]
            
        # Correction of var_l and var_t for var_n
        corr_var_l[ind_corr] = var_l[0][ind_corr]-var_n[ind_corr]
        corr_var_t[ind_corr] = var_t[0][ind_corr]-var_n[ind_corr]

        # From variance to seeing: corrected fwhm_l and fwhm_t
        corr_fwhm_l = self.var2seeing(corr_var_l, 'l')
        corr_fwhm_t = self.var2seeing(corr_var_t, 't')

        return corr_fwhm_l, corr_fwhm_t
        
              
