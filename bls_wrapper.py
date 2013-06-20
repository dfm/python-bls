# BLS wrapper
# adjustable parameters: convolving window; median filter length
# Could try using a running mean filter

import BLS
import pyfits
import numpy
import scipy
import scipy.ndimage.filters
import pylab as p
from math import *


def run_BLS():

    print 'Loading lc...'
    time, lc = load_data(type = 'ascii')
    
    print 'Median normalise...'
    lc = lc/numpy.median(lc)

    #print 'Median filter...'
    #lc = medfilt(lc,97) # not shorter than the transit duration, maybe 2?

    print '1st pass bls...'
    ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
        qtran, duration, f_1, convolved_bls = compute_bls(time, lc)
    
    print 'Cutting out 1st planet transits...'
    time, lc = cut_out_transits(time, lc, ingresses, egresses)

   
    p.close(1)
    p.figure(1)
    p.subplot(2,1,1)
    p.plot(time,lc,".k")
    p.ylabel('Flux')
    p.xlabel('Time')
    #p.xlim(min(time), max(time))
    for i in range(len(ingresses)):
        p.axvline(ingresses[i], color = 'c')
        p.axvline(egresses[i], color = 'c')
    p.subplot(2,1,2)
    p.plot(f_1, convolved_bls)
    #p.axvline(periods, color = 'r')
    p.axvline(bper, color = '0.5')
    print 'period measurements = ', periods, bper
    p.xlim(min(f_1), max(f_1))
    p.xlabel('Period')
    # p.subplot(3,1,3)
    # p.plot(phases, lc, 'k.')
    # p.xlabel('Phase')
    # p.ylabel('Flux')

    print '2nd pass bls...'
    ingresses, egresses, t_number, epoch, periods, bper, bpow, depth, qtran, \
        qtran, duration, f_1, convolved_bls = compute_bls(time, lc)

    p.close(2)
    p.figure(2)
    p.subplot(2,1,1)
    p.plot(time,lc,".k")
    p.ylabel('Flux')
    p.xlabel('Time')
    #p.xlim(min(time), max(time))
    for i in range(len(ingresses)):
        p.axvline(ingresses[i], color = 'c')
        p.axvline(egresses[i], color = 'c')
    p.subplot(2,1,2)
    p.plot(f_1, convolved_bls)
    #p.axvline(periods, color = 'r')
    p.axvline(bper, color = '0.5')
    print 'period measurements = ', periods, bper
    p.xlim(min(f_1), max(f_1))
    p.xlabel('Period')
    # p.subplot(3,1,3)
    # p.plot(phases, lc, 'k.')
    # p.xlabel('Phase')
    # p.ylabel('Flux')

    return 


def load_data(type):
    
    if type == 'Fits':
        file = '/Users/angusr/.kplr/data/old/kplr006448890-2009259160929_llc.fits'
        hdulist = pyfits.open(file)
        tbdata = hdulist[1].data 

        '''Remove NANs'''
        x = numpy.where(numpy.isfinite(tbdata['TIME']))
        time = tbdata['TIME'][x]; lc = tbdata['PDCSAP_FLUX'][x]
        x = numpy.where(numpy.isfinite(tbdata['PDCSAP_FLUX']))
        time = tbdata['TIME'][x]; lc = tbdata['PDCSAP_FLUX'][x]

    elif type == 'ascii':
        
        #f = open('/Users/angusr/angusr/data2/SAMSI/ascii_inj/k7372635.dat', 'r')
        f = open('/Users/angusr/angusr/data2/SAMSI/KIC_005383248_long.dat', 'r')
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

    #return numpy.array(time[0:3000]), numpy.array(lc[0:3000])
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


# def smooth(x,window_len=11,window='hanning'):
#     if x.ndim != 1:
#         raise ValueError, "smooth only accepts 1 dimension arrays."
#     if x.size < window_len:
#         raise ValueError, "Input vector needs to be bigger than window size."
#     if window_len<3:
#         return x
#     if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
#         raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
#     s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
#     #print(len(s))
#     if window == 'flat': #moving average
#         w=numpy.ones(window_len,'d')
#     else:
#         w=eval('numpy.'+window+'(window_len)')
#     y=numpy.convolve(w/w.sum(),s,mode='valid')
#     return y


def compute_bls(time, lc):
    '''Calculate BLS'''
    bls, f_1, nb = BLS.BLS(time, lc, fmin = (1./(50.0*1.1)))
    print 'Complete'

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

    return ingresses, egresses, t_number, epoch, periods, bls[1], bls[2], bls[3], bls[4], \
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

    return time, lc


if __name__ == '__main__':
    test ()


    #lc[time % period < duration] *= (1 - depth)


    # period = 10.5
    # duration = 0.1
    # depth = 1e-3
    # time = numpy.linspace(0,200,5000)
    # time = (numpy.random.rand(5000))*200
    # time = numpy.sort(time)
    # lc = numpy.ones(len(time))
    # lc[time % period < duration] *= (1 - depth)
    # lc += .5e-3 * numpy.random.randn(len(lc))


     # p.close(2)
    # p.figure(2)
    # p.subplot(3,1,1)
    # p.plot(time, lc)
    # p.subplot(3,1,2)
    # p.plot (medfilt (lc,3))
    # p.subplot(3,1,3)
    # p.plot (medfilt (lc,5))


    # '''locate 3 highest peaks'''
    # periods = numpy.zeros(3)
    # for i in range(0, len(sel_peaks)):
    #     for j in range(0, len(convolved_bls)):
    #         if sel_peaks[i] == convolved_bls[j]:
    #             periods[i] = f_1[j]


    
    # ''' Fold '''
    # phases = (time-time[0]) % bper
    # phases = ( (time - time[0]) / bper)
    # for i in range(0, len(phases)):
    #     phases[i] -= int(phases[i])
