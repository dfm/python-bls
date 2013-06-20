# BLS wrapper

import BLS
import pyfits
import numpy
import scipy
import scipy.ndimage.filters
import pylab as p
from math import *
import bls
#import filter
from numpy.random import normal
import atpy

gap_days = 0.02043365  # long cadence

ascii_DIR = '/Users/angusr/angusr/data2/SAMSI'
ascii_file = 'KIC_005383248_long.dat'

fits_DIR = '/Users/angusr/.kplr/data/old'
fits_file = 'kplr006448890-2009259160929_llc.fits'

def run_BLS():

    print 'Loading lc...'
    time, lc = load_data(type = 'ascii')
    
    print 'Median normalise...'
    lc = lc/numpy.median(lc)

    print '1st pass bls...'
    ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
        duration, f_1, convolved_bls = compute_bls(time, lc)

   
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
    p.axvline(bper, color = '0.5')
    for i in range(2, 10):
        p.axvline(i*bper, color = 'y')
    p.axvline(bper/2., color = 'y')
    p.axvline(3*bper/2., color = 'y')
    print 'period = ', bper
    p.xlim(min(f_1), max(f_1))
    p.xlabel('Period')


    # print 'Cutting out 1st planet transits...'
    # time, lc = cut_out_transits(time, lc, ingresses, egresses)

    # print 'Interpolating...'
    # time, lc = interpolate(time, lc)

    # p.close(3)
    # p.figure(3)
    # p.plot(time,lc, 'k.')
    # for i in range(len(ingresses)):
    #     p.axvline(ingresses[i], color = 'c')
    #     p.axvline(egresses[i], color = 'c')

    # print '2nd pass bls...'
    # ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
    #     duration, f_1, convolved_bls = compute_bls(time, lc)

    # p.close(4)
    # p.figure(4)
    # p.subplot(2,1,1)
    # p.plot(time,lc,".k")
    # p.ylabel('Flux')
    # p.xlabel('Time')
    # for i in range(len(ingresses)):
    #     p.axvline(ingresses[i], color = 'c')
    #     p.axvline(egresses[i], color = 'c')
    # p.subplot(2,1,2)
    # p.plot(f_1, convolved_bls)
    # p.axvline(bper, color = '0.5')
    # print 'period measurements = ', periods, bper
    # p.xlim(min(f_1), max(f_1))
    # p.xlabel('Period')
    
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
            lc.append(float(columns[2][0:-1]))
        time = numpy.array(time); lc = numpy.array(lc)

        x = numpy.where(numpy.isfinite(time))
        time = time[x[0]]; lc = lc[x[0]]
        x = numpy.where(numpy.isfinite(lc))
        time = time[x[0]]; lc = lc[x[0]]

    return numpy.array(time), numpy.array(lc)


def compute_bls(time, lc):
    '''Calculate BLS'''
    bls, f_1, nb = BLS(time, lc)
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

    return ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
        duration, f_1, convolved_bls



def cut_out_transits(time, lc, ingresses, egresses):
    new_time = list(time); new_lc = list(lc)
    for j in range(0, len(ingresses)-1):
        for i in range(0,len(time)):
            if ingresses[j] < time[i] < egresses[j]:
                new_time.remove(time[i])
                new_lc.remove(lc[i])

    time = numpy.array(new_time)
    lc = numpy.array(new_lc)

    # p.close(2)
    # p.figure(2)
    # p.plot(time,lc, 'k.')
    # for i in range(len(ingresses)):
    #     p.axvline(ingresses[i], color = 'c')
    #     p.axvline(egresses[i], color = 'c')

    return time, lc


''' df = frequency step, \
    nf = number of frequencies, \
    nb = number of bins, \
    qmi = minimum fractional transit duration, \
    qma = maximum transit duration, \
    fmin = minimum frequency '''

def BLS(time, lc, df = 0.0001, nf = 500,  nb = 200, qmi = 0.01,\
        qma = 0.8, fmin = (1./(400.0*1.1))):

    u = numpy.ones(len(time))
    v = numpy.ones(len(time))

    BLS = bls.eebls(time, lc, u, v, nf, fmin, df, nb, qmi, qma)
    f = fmin + (numpy.arange(len(BLS[0])))*df
    
    return BLS, 1/f, nb


def interpolate(time, flux):
    
    ''' Calculate noise properties '''
    flux_filt = filter.filt1d(flux, 10, 5)
    med_noise, sig_noise = filter.medsig(flux-flux_filt)


    ''' find gaps greater than 1 (1.1) data points and linear interp with noise'''
    diff1 = time[1:] - time[:-1]
    diff1 = scipy.append(diff1, gap_days)
    gap_find = diff1 > 1.1*gap_days
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
        if time_end - fill.max() > 1.1*gap_days: \
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


    return time, flux

if __name__ == '__main__':
    BLS()
