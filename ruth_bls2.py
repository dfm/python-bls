# BLS wrapper

import pyfits
import numpy
import scipy
import scipy.ndimage.filters
import pylab as p
from math import *
import bls
import filter
from numpy.random import normal
import atpy
import medium_filter

print 'Dan rocks! x2'

gap_days = 0.02043365  # long cadence
#ascii_DIR = '/Users/angusr/angusr/data2/SAMSI/ascii_inj'
ascii_DIR = '/Users/angusr/angusr/data2/SAMSI'
ascii_file = 'kplr000892203-2009131105131_llc.txt' 
#ascii_file = 'KIC_005383248_long.dat'
#ascii_file = 'k7372635.dat'
ascii_DIR = '/Users/angusr/angusr/data2/SAMSI/ascii_inj/detrended'
ascii_file = 'k2161400.dat'
#ascii_file = 'nopl_k2161400.dat'

fits_DIR = '/Users/angusr/.kplr/data/old'
fits_file = 'kplr006448890-2009259160929_llc.fits'

def run_BLS(Baines = False, type = 'ascii'):

    print 'Loading lc...'
    time, lc = load_data(type)

    p.close(6)
    p.figure(6)
    p.plot(time,lc)
    
    print 'Median normalise...'
    lc = lc/numpy.median(lc)

    print 'Median filter...'
    lc = medium_filter.detrend_using_mf(lc)

    ''' Define parameters '''
    min_period = 300
    max_period = 450
    freq_step = 0.000001  #0.0001
    min_duration_hours = 10
    max_duration_hours = 15
    
    print '1st pass bls...'
    ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
        duration, f_1, convolved_bls, approx_duration = compute_bls(time, \
                lc, df = freq_step, \
                nf =  (1./min_period)/freq_step, nb = 1400, \
                qmi = float(min_duration_hours)/24./450., \
                qma = float(max_duration_hours)/24./300., \
            fmin = (1./(float(max_period)*1.1)))

   
    p.close(1)
    p.figure(1)
    p.subplot(2,1,1)
    p.plot(time,lc,".k")
    p.ylabel('Flux')
    p.xlabel('Time')
    for i in range(len(ingresses)):
        p.axvline(ingresses[i], color = 'c')
        p.axvline(egresses[i], color = 'c')
    p.subplot(2,1,2)
    p.plot(f_1, convolved_bls)
    p.axvline(bper, color = 'r')
    for i in range(2, 10):
        p.axvline(i*bper, color = 'y')
    p.axvline(bper/2., color = 'y')
    p.axvline(3*bper/2., color = 'y')
    print 'period = ', bper
    p.xlim(min(f_1), max(f_1))
    p.xlabel('Period')

    if Baines == False:

        p.close(2)

        '''Folding lc ...'''
        Fold(time, lc, bper, ingresses, egresses, plot_time = numpy.zeros(len(time)), \
             plot_lc = numpy.zeros(len(time)), figure = 2)

        print 'Cutting out 1st planet transits...'
        time, lc, plot_time, plot_lc = cut_out_transits(time, lc, ingresses, egresses)

        '''Folding lc ...'''
        Fold(time, lc, bper, ingresses, egresses, plot_time, plot_lc, figure = 3)

        print 'Interpolating...'
        time, lc = interpolate(time, lc, approx_duration)

        '''Folding lc ...'''
        Fold(time, lc, bper, ingresses, egresses, plot_time, plot_lc, figure = 4)

        freq_step =  0.000001

        print '2nd pass bls...'
        # ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
        #     duration, f_1, convolved_bls = compute_bls(time, lc)
        ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
            duration, f_1, convolved_bls, approx_duration = compute_bls(time, \
                lc, df = freq_step, nf =  (1./min_period)/freq_step,\
                nb = 1400, qmi = float(min_duration_hours)/24./450., \
                qma = float(max_duration_hours)/24./300., \
                fmin = (1./(float(max_period)*1.1)))

        p.close(5)
        p.figure(5)
        p.subplot(2,1,1)
        p.plot(time,lc,".k")
        p.ylabel('Flux')
        p.xlabel('Time')
        for i in range(len(ingresses)):
            p.axvline(ingresses[i], color = 'c')
            p.axvline(egresses[i], color = 'c')
        p.subplot(2,1,2)
        p.plot(f_1, convolved_bls)
        p.axvline(bper, color = '0.5')
        print 'period measurements = ', periods, bper
        p.xlim(min(f_1), max(f_1))
        p.xlabel('Period')

        p.close(3)

        '''Folding lc ...'''
        Fold(time, lc, bper, ingresses, egresses, plot_time = numpy.zeros(len(time)), \
             plot_lc = numpy.zeros(len(time)), figure = 2)

        print 'Cutting out 2nd planet transits...'
        time, lc, plot_time, plot_lc = cut_out_transits(time, lc, ingresses, egresses)

        '''Folding lc ...'''
        Fold(time, lc, bper, ingresses, egresses, plot_time, plot_lc, figure = 3)

        print 'Interpolating...'
        time, lc = interpolate(time, lc, approx_duration)

        '''Folding lc ...'''
        Fold(time, lc, bper, ingresses, egresses, plot_time, plot_lc, figure = 4)
    
    return 


def load_data(type):
    
    if type == 'fits':
        file = '%s/%s' %(fits_DIR, fits_file)
        hdulist = pyfits.open(file)
        tbdata = hdulist[1].data 

        '''Remove NANs'''
        x = numpy.where(numpy.isfinite(tbdata['TIME']))
        time = tbdata['TIME'][x]; lc = tbdata['PDCSAP_FLUX'][x]
        x = numpy.where(numpy.isfinite(tbdata['PDCSAP_FLUX']))
        time = tbdata['TIME'][x]; lc = tbdata['PDCSAP_FLUX'][x]

    elif type == 'ascii':
        f = open('%s/%s' %(ascii_DIR, ascii_file), 'r')
        time = []
        lc = []
        for line in f:
            line = line.strip()
            columns = line.split()
            time.append(float(columns[0][0:-1]))
            lc.append(float(columns[1][0:-1]))
        time = numpy.array(time); lc = numpy.array(lc)

        x = numpy.where(numpy.isfinite(time))
        time = time[x[0]]; lc = lc[x[0]]
        x = numpy.where(numpy.isfinite(lc))
        time = time[x[0]]; lc = lc[x[0]]

    return numpy.array(time), numpy.array(lc)

# def medfilt (x, k):
#     """Apply a length-k median filter to a 1D array x.
#     Boundaries are extended by repeating endpoints.
#     """
#     assert k % 2 == 1, "Median filter length must be odd."
#     assert x.ndim == 1, "Input must be one-dimensional."
#     k2 = (k - 1) // 2
#     y = numpy.zeros ((len (x), k), dtype=x.dtype)
#     y[:,k2] = x
#     for i in range (k2):
#         j = k2 - i
#         y[j:,i] = x[:-j]
#         y[:j,i] = x[0]
#         y[:-j,-(i+1)] = x[j:]
#         y[-j:,-(i+1)] = x[-1]
#     return numpy.median (y, axis=1)


def compute_bls(time, lc, df, nf, nb, qmi, qma, fmin):
    '''Calculate BLS'''
    bls, f_1, nb = BLS(time, lc, df, nf, nb, qmi, qma, fmin)
    print 'Complete'

    bper = bls[1]
    bpow = bls[2]
    depth = bls[3]
    qtran = bls[4]
    duration = bper*qtran
    in1 = bls[5]
    in2 = bls[6]
    phase1 = in1/float(nb)
    phase2 = in2/float(nb)
    
    '''Convolve with gaussian'''
    convolved_bls = scipy.ndimage.filters.gaussian_filter(bls[0], 2.0) # maybe adjust this?\
                                                                       #related to df
                                                                       #Play with this
    '''Locate all peaks'''
    peak = numpy.r_[True, convolved_bls[1:] > convolved_bls[:-1]] & \
       numpy.r_[convolved_bls[:-1] > convolved_bls[1:], True]

    '''Sort peaks'''
    sel_peaks = numpy.sort(convolved_bls[peak])
    sel_peaks = sel_peaks[-1:]

    '''locate highest peak'''
    periods = f_1[numpy.where(convolved_bls == sel_peaks)]

    '''calculate number of transits, epoch, ingress and egress times'''
    t_number = int((max(time) - min(time)) / bper)
    epoch = time[0] + phase1*bper
    ingresses = numpy.zeros(t_number)
    egresses = numpy.zeros(t_number)
    for n in range(0,t_number):
        ingresses[n] = (epoch + bper*n) - 0.2
        egresses[n] = epoch + bper*n + duration + 0.2 # add a margin each side

    approx_duration = egresses[0] - ingresses[0]

    return ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
        duration, f_1, convolved_bls, approx_duration

#---------------------------------------------------------------------------------------

''' df = frequency step, \
    nf = number of frequencies, \
    nb = number of bins, \
    qmi = minimum fractional transit duration, \
    qma = maximum transit duration, \
    fmin = minimum frequency '''

def BLS(time, lc, df = 0.0001, nf = 500,  nb = 200, qmi = 0.01,\
        qma = 0.8, fmin = (1./(400.0*1.1))):
    diffs = time[1:] - time[:-1]

    u = numpy.ones(len(time))
    v = numpy.ones(len(time))
    
    BLS = bls.eebls(time, lc, u, v, nf, fmin, df, nb, qmi, qma)
    f = fmin + (numpy.arange(len(BLS[0])))*df
    
    return BLS, 1/f, nb

#---------------------------------------------------------------------------------------

def cut_out_transits(time, lc, ingresses, egresses):
    
    new_time = list(time); new_lc = list(lc)
    plot_time = []; plot_lc = []
    for j in range(0, len(ingresses)):
        for i in range(0,len(time)):
            if ingresses[j] < time[i] < egresses[j]:
                new_time.remove(time[i])
                new_lc.remove(lc[i])
                plot_time.append(time[i])
                plot_lc.append(lc[i])

    time = numpy.array(new_time)
    lc = numpy.array(new_lc)

    # p.close(2)
    # p.figure(2)
    # p.plot(time,lc, 'k.')
    # for i in range(len(ingresses)):
    #     p.axvline(ingresses[i], color = 'c')
    #     p.axvline(egresses[i], color = 'c')

    return time, lc, plot_time, plot_lc

#---------------------------------------------------------------------------------------
def Fold(time, lc, bper, ingresses, \
         egresses, plot_time, plot_lc, figure ):
    ''' Fold '''
    #phases = (time-time[0]) % bper
    ingress = (ingresses[0] - time[0]) / bper
    egress = (egresses[0] - time[0]) / bper
    plot_phase = (plot_time - plot_time[0]) / bper
    for i in range(0,len(plot_time)):
        plot_phase[i] -= int(plot_phase[i])
    phases =  (time - time[0]) / bper
    for i in range(0, len(phases)):
        phases[i] -= int(phases[i])

    # p.close(figure)
    # p.figure(figure)
    #p.close(2)
    p.figure(2)
    #p.subplot(2,1,1)
    p.subplot(3,1,figure-1)
    p.plot(phases, lc, 'k.')
    if figure == 3:
        p.plot(plot_phase, plot_lc, 'r.')
    p.axvline(ingress, color = 'c')
    p.axvline(egress, color = 'c')
    #p.xlabel('Phase')
    p.ylabel('Flux')
    p.xlim(ingress - (egress-ingress)*2, egress + (egress - ingress)*2)
    # #p.subplot(2,1,2)
    # p.plot(time, lc, 'k.')
    # p.xlabel('Time (days)')
    # p.ylabel('Flux')
    # for i in range(len(ingresses)):
    #     p.axvline(ingresses[i], color = 'c')
    #     p.axvline(egresses[i], color = 'c')
    # if figure == 3:
    #     p.plot(plot_time, plot_lc, 'r.')
    # p.xlim(min(time), max(time))

    return

#---------------------------------------------------------------------------------------

def interpolate(time, flux, approx_duration):

    min_gap = approx_duration/gap_days

    # p.close(5)
    # p.figure(5)
    # p.subplot(2,1,1)
    # p.plot(time, flux, 'k.')
    
    ''' Calculate noise properties '''
    flux_filt = filter.filt1d(flux, 10, 5)
    med_noise, sig_noise = filter.medsig(flux-flux_filt)


    ''' find gaps greater than 1 (1.1) data points and linear interp with noise'''
    diff1 = time[1:] - time[:-1]
    diff1 = scipy.append(diff1, gap_days)
    gap_find = diff1 > min_gap*gap_days*0.9 #1.1*gap_days
    gap_times = time[gap_find]
    time_index = scipy.r_[0:len(time):1]

    fill_arr_t = scipy.zeros(1)
    fill_arr_f = scipy.zeros(1)
    #fill_arr_nan = scipy.zeros(1)

    print 'Filling gaps...'
    for m in scipy.arange(len(gap_times)):
        time_start = time[time_index[gap_find][m]]
        flux_start = flux[time_index[gap_find][m]]
        time_end = time[time_index[gap_find][m] + 1]
        flux_end = flux[time_index[gap_find][m] + 1]
        span =  time_end - time_start
        if span < 2.0*gap_days:
            fill = scipy.array([time_start, time_start+gap_days, time_end])
        else: fill = scipy.r_[time_start: time_end: gap_days]
        fill = fill[1:-1]
        if time_end - fill.max() > 1.1*gap_days:  #1.1*gap_days*min_gap: 
            fill = scipy.append(fill, fill.max()+gap_days)

        f = scipy.interpolate.interp1d([time_start, time_end], [flux_start, flux_end])
        gap_new = f(fill)

       
        if sig_noise > 0: gap_new += normal(0,sig_noise,len(fill))

        fill_arr_t = scipy.append(fill_arr_t, fill)
        fill_arr_f = scipy.append(fill_arr_f, gap_new)
        #fill_arr_nan = scipy.append(fill_arr_nan, scipy.ones(len(fill))*scipy.nan)

    # combine time and flux arrays with their filled sections
    fill_arr_t = fill_arr_t[1:]
    fill_arr_f = fill_arr_f[1:]
    #fill_arr_nan = fill_arr_nan[1:]
    fill_arr_t = scipy.append(fill_arr_t, time)
    fill_arr_f = scipy.append(fill_arr_f, flux)
    #fill_arr_nan = scipy.append(fill_arr_nan, flux_base)

    if len(fill_arr_t) == 0:
        print '*************** empty time array ***************'
        return False, atpy.Table(), 0, 0
    
    # put in table and sort
    tf = atpy.Table()
    tf.add_column('time', fill_arr_t)
    tf.add_column('flux', fill_arr_f)
    # tf.add_column('flux_pdc', fill_arr_nan)
    tf.sort('time')

    # for i in range(0, len(gap_times)):
    #     p.axvline(gap_times[i], color = 'y')
    # p.subplot(2,1,2)
    # p.plot(fill_arr_t, fill_arr_f, 'k.')
        

    return tf.time, tf.flux

if __name__ == '__main__':
    BLS()
