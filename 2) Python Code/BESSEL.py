"""
=====================================
BESSEL THERMAL DIFFUSIVITY CALCULATOR
=====================================

Hello fellow physicist! 

This script was designed to compute thermal diffusivity from the values of 
transmission factor and phase lag that you have obtained from the Fourier 
analysis of the data.

It uses Bessel functions of first kind, J0, which can be treated as a series 
expansion. More about this can be found in the Script's Appendix.

In this program, we can choose the number of terms that we wish to truncate J0 
to - this is controlled by the variable 'no_terms'. By default, we set 
'no_terms' = 100, which should work in most cases. If you get an error message 
saying that "the method failed to converge", lower the variable 'no_terms', 
as suggested by the message outputed. Sometimes working with too many terms 
can give rise to convergence problems!

You may use this script as a calculator, but you are very welcome to have a 
look at what it is actually doing by scrolling down a bit :)

To start, just input your data and press run!


Created on Wed Aug 19 18:46:50 2020
@author: Carolina_Rossi
"""

# =============================================================================
# YOUR DATA:
# =============================================================================





# #fourA
# trans_data = [0.1472444997171886, 0.013626831060655238, 0.0028642566550574106, 0.0004718937373537813]#input your transmission data here
# trans_periods = [240, 240/3, 240/5, 240/7] #input the corresponding periods here (s)
# phase_data = [3.150437387016244, 5.6703702430308445, 0.7796292901114679, 2.295166807313928]#input your phase lag data here
# phase_periods = [240, 240/3, 240/5, 240/7]#input the corresponding periods here (s)
# r_out = 0.01012 #outer radius (m)
# r_in = 0.0025 #inner  radius (m)
# no_terms = 100 #number of terms you want to truncate the J0 Bessel Series to

# # oneA
# trans_data = [0.017615699606359638, 0.0002533180592878783, 0.0008329137124024899, 0.0006109599195345082]#input your transmission data here
# trans_periods = [60/ 1, 60/3, 60/5, 60/7] #input the corresponding periods here (s)
# phase_data = [5.358726887041305, 0.9477543334237679, 0.7158429269928867, 2.5853489681912682]#input your phase lag data here
# phase_periods = [60/ 1, 60/3, 60/5, 60/7]#input the corresponding periods here (s)
# r_out = 0.01012 #outer radius (m)
# r_in = 0.0025 #inner  radius (m)
# no_terms = 100 #number of terms you want to truncate the J0 Bessel Series to

# # twoA
# trans_data = [0.03845372601109291, 0.0021360894974493303, 0.0014936274537170275, 0.0015300197592436422]#input your transmission data here
# trans_periods = [120/1, 120/3, 120/5, 120/7] #input the corresponding periods here (s)
# phase_data = [4.53991191723654, 0.4545742737989054, 6.274496130964069, 6.23966060601856]#input your phase lag data here
# phase_periods = [120/1, 120/3, 120/5, 120/7]#input the corresponding periods here (s)
# r_out = 0.01012 #outer radius (m)
# r_in = 0.0025 #inner  radius (m)
# no_terms = 100 #number of terms you want to truncate the J0 Bessel Series to

# # six
# trans_data = [0.27661764589088883, 0.04166375470563922, 0.012048289828457195, 0.005561523583587996]#input your transmission data here
# trans_periods = [360/1, 360/3, 360/5, 360/7] #input the corresponding periods here (s)
# phase_data = [2.4785852663689494, 4.677025559738848, 6.081069692252148, 0.5755948054482881]#input your phase lag data here
# phase_periods = [360/1, 360/3, 360/5, 360/7]#input the corresponding periods here (s)
# r_out = 0.01012 #outer radius (m)
# r_in = 0.0025 #inner  radius (m)
# no_terms = 100 #number of terms you want to truncate the J0 Bessel Series to

# # eight
# trans_data = [0.4432412338094975, 0.11312832803755474, 0.03287289664907653, 0.012294316998891958]#input your transmission data here
# trans_periods = [480/1, 480/3, 480/5, 480/7] #input the corresponding periods here (s)
# phase_data = [1.8480068762963306, 3.638720018483986, 4.843820096244032, 6.03419487110334]#input your phase lag data here
# phase_periods = [480/1, 480/3, 480/5, 480/7]#input the corresponding periods here (s)
# r_out = 0.01012 #outer radius (m)
# r_in = 0.0025 #inner  radius (m)
# no_terms = 100 #number of terms you want to truncate the J0 Bessel Series to

# sixteen
trans_data = [0.6016240363068186, 0.13928434185643934, 0.04153834995788276, 0.0750101076711088]#input your transmission data here
trans_periods = [960/1, 960/3, 960/5, 960/7] #input the corresponding periods here (s)
phase_data = [1.1742599987696252, 2.4762311997472244, 4.786709655775473, 5.765026917903765]#input your phase lag data here
phase_periods = [960/1, 960/3, 960/5, 960/7]#input the corresponding periods here (s)
r_out = 0.01012 #outer radius (m)
r_in = 0.0025 #inner  radius (m)
no_terms = 100 #number of terms you want to truncate the J0 Bessel Series to




# =============================================================================   
# BESSEL TRANSMISSION AND BESSEL PHASE LAG FUNCTIONS
# =============================================================================
# General comments:
# -----------------
#   Terms in the Bessel Series are complex. They take the form:
#       -  J0 = re_in + i(im_in) ,  at r_in;
#       -  J0 = re_out + i(im_out) ,  at r_out.
#   We consider J0 as a function of alpha, where alpha = sqrt(w/D)*r_in.
#
#   The functions are designed to be used for fitting the transmission factors
#   and phase lags obtained in the Fourier analysis of the experimental data.
#    
#   A Secant root-finding method is used to find the alpha than matches
#   the Bessel phase lag / transmision values to the Fourier ones. 
#        
#   Finally, the diffusivity is obtained from alpha.
#     
# About the iterator I:
# --------------------
#   The Bessel series is obtained via an iterative process. 
#   Each iteration corresponds to 2 terms in the series, so upper limit is set 
#   to be I < (no_terms/2).
# 
# =============================================================================


from scipy.optimize import newton
import numpy as np

def bessel_trans_fit(alpha): 
    """
    FUNCTION THAT CALCULATES TRANSMISSION FACTORS VIA BESSEL ANALYSIS 
    
    It computes the J0 at the inner radius (r_in), the outer radius (r_out), 
    and the transmission factor J0_in/J0_out (bessel_trans).
    
    Parameters:
    ----------
    alpha = sqrt(w/D)*r_in      
    
    Returns:
    -------
    For a given Fourier transmission value (fourier_trans), bessel_trans_fit 
    returns:
            bessel_trans - fourier_trans   
    Later, a secant root-finding method is used to find alpha when
    bessel_trans - fourier_trans = 0. 
    """
    I = 0 
    ratio = r_out/r_in
    re_in = re_out = im_in = im_out = 0
    
    #TO UPDATE VARIABLE PART OF TERMS IN SERIES:
    term_in = 1  #r_in series
    term_out = 1 #r_out series
    
    #FOR UPDATING CONSTANT COEFFICIENT OF TERMS IN SERIES:
    re_coeff = 1  #real part
    im_coeff = -re_coeff / 4 #imaginary part

    while I<(no_terms/2):
        #updating r_in series:
        re_in = re_in + re_coeff * term_in
        term_in = term_in * alpha**2
        im_in = im_in + im_coeff * term_in
        term_in = term_in * alpha**2 
        
        #updating r_out series:
        re_out = re_out + re_coeff * term_out
        term_out = term_out * (alpha * ratio)**2 
        im_out = im_out + im_coeff * term_out
        term_out = term_out * (alpha * ratio)**2
        
        #updating coeffs and iterator:
        re_coeff = -re_coeff / (4*(2*I + 1)*(2*I + 2))**2
        I = I + 1 
        im_coeff = -re_coeff / (4*(2*I + 1)**2)
           
    #transmission factor:
    bessel_trans = np.sqrt((re_in**2 + im_in**2) / (re_out**2 + im_out**2))
    
    return bessel_trans - fourier_trans


def bessel_phase_fit(alpha):
    """
    FUNCTION THAT CALCULATES PHASE LAGS VIA BESSEL ANALYSIS 
    
    It computes the J0 at the inner radius (r_in), the outer radius (r_out), 
    and the phase lag, phase_in - phase_out (bessel_phase_lag).
    
    Parameters:
    ----------
    alpha = sqrt(w/D)*r_in      
    
    Returns:
    -------
    For a given Fourier phase lag value (fourier_phase_lag),bessel_phase_fit 
    returns:
            bessel_phase_lag - fourier_phase_lag 
    Later, a secant root-finding method is used to find alpha when
    bessel_phase_lag - fourier_phase_lag = 0. 
    """
    I = 0 
    ratio = r_out/r_in
    re_in = re_out = im_in = im_out = 0
    
    #TO UPDATE VARIABLE PART OF TERMS IN SERIES:
    term_in = 1  #r_in series
    term_out = 1 #r_out series
    
    #FOR UPDATING CONSTANT COEFFICIENT OF TERMS IN SERIES:
    re_coeff = 1  #real part
    im_coeff = -re_coeff / 4 #imaginary part

    while I<(no_terms/2):
        #updating r_in series:
        re_in = re_in + re_coeff * term_in
        term_in = term_in * alpha**2
        im_in = im_in + im_coeff * term_in
        term_in = term_in * alpha**2 
        
        #updating r_out series:
        re_out = re_out + re_coeff * term_out
        term_out = term_out * (alpha * ratio)**2 
        im_out = im_out + im_coeff * term_out
        term_out = term_out * (alpha * ratio)**2
        
        #updating coeffs and iterator:
        re_coeff = -re_coeff / (4*(2*I + 1)*(2*I + 2))**2
        I = I + 1 
        im_coeff = -re_coeff / (4*(2*I + 1)**2)
        
    phase_out = np.arctan2(im_out,re_out)
    phase_in = np.arctan2(im_in,re_in)
    
    if phase_in > phase_out:
         bessel_phase_lag = phase_in - phase_out 
    else:
         bessel_phase_lag = phase_in - phase_out + 2*np.pi 
    
    return bessel_phase_lag - fourier_phase_lag



# =============================================================================
# FITTING FUNCTIONS TO EXPERIMENTAL DATA AND CALCULATING DIFFUSIVITY:
# =============================================================================

#FITTING bessel_trans_fit TO TRANSMISSION DATA 
trans_alphas = []
for trans in trans_data:
    fourier_trans = trans
    b = newton(bessel_trans_fit,x0=1) #secant method 
    trans_alphas.append(b)

#FITTING bessel_phase_fit TO PHASE DATA 
phase_alphas = []
for phase in phase_data:
    fourier_phase_lag = phase
    b = newton(bessel_phase_fit,x0=1) #secant method
    phase_alphas.append(b)

#CALCULATING PHASE DIFFUSIVITY FROM ALPHA
D_phase = []
for i in range(len(phase_alphas)):
    w = (2*np.pi/phase_periods[i])
    alpha = phase_alphas[i]
    D_p = (w*r_in**2)/(alpha**2)
    D_phase.append(D_p)
    
#CALCULATING TRANSMISSION DIFFUSIVITY FROM ALPHA
D_trans = []
for i in range(len(trans_alphas)):
    w = (2*np.pi/trans_periods[i])
    alpha = trans_alphas[i]
    D_t = (w*r_in**2)/(alpha**2)
    D_trans.append(D_t)


print("\nDIFFUSIVITY - BESSEL ANALYSIS")
print('\nDiffusivity - Transmission Factors:')
print(D_trans)
print('\nDiffusivity - Phase Lags:')
print(D_phase)
