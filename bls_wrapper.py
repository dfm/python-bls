# BLS wrapper
# parameters to play with: convolving window

import BLS
import pyfits
import numpy
import scipy
import scipy.ndimage.filters
import pylab as p


def run_BLS():
    # period = 10.5
    # duration = 0.1
    # depth = 1e-3
    # time = numpy.linspace(0,200,5000)
    # time = (numpy.random.rand(5000))*200
    # time = numpy.sort(time)
    # lc = numpy.ones(len(time))
    # lc[time % period < duration] *= (1 - depth)
    # lc += .5e-3 * numpy.random.randn(len(lc))

    # Load data
    file = '/Users/angusr/.kplr/data/old/kplr006448890-2009259160929_llc.fits'
    hdulist = pyfits.open(file)
    tbdata = hdulist[1].data #"add the first extension to tbdata"
    x = numpy.where(numpy.isfinite(tbdata['TIME']))
    time = tbdata['TIME'][x]
    lc = tbdata['PDCSAP_FLUX'][x]
    x = numpy.where(numpy.isfinite(tbdata['PDCSAP_FLUX']))
    time = tbdata['TIME'][x]
    lc = tbdata['PDCSAP_FLUX'][x]
    lc = lc/numpy.median(lc)

    # Median filter
    lc = medfilt(lc,1)

    # Calculate BLS
    bls, f_1 = BLS.BLS(time, lc, fmin = (1./(50.0*1.1)))
    #print 'Max BLS = ', max(bls[0])
    #x = numpy.where(bls[0] == max(bls[0]))
    #period = f_1[x]

    # Convolve gaussian
    convolved_bls = scipy.ndimage.filters.gaussian_filter(bls[0], 1.0)
    peak = numpy.r_[True, convolved_bls[1:] < convolved_bls[:-1]] & \
        numpy.r_[convolved_bls[:-1] < convolved_bls[1:], True]

    # print 'peaks = ', convolved_bls[peak]
    #x = numpy.where(convolved_bls[peak] > 2.5e-5)
    #sel_peaks = numpy.sort(convolved_bls[peak][x])

    sel_peaks = numpy.sort(convolved_bls[peak])
    print sel_peaks[-3:]
    x = numpy.where(convolved_bls[peak] == sel_peaks[-1], sel_peaks[-2], \
                    sel_peaks[-3])
    print f_1[x]


    #print 'Period = ', period

    p.close(1)
    p.figure(1)
    p.subplot(2,1,1)
    p.plot(time,lc,".k")
    p.ylabel('Flux')
    p.xlabel('Time')
    p.xlim(min(time), max(time))
    p.subplot(2,1,2)
    #p.plot(f_1, bls[0])
    p.plot(f_1, convolved_bls)
    p.xlim(min(f_1), max(f_1))
    p.xlabel('Period')

    # p.close(2)
    # p.figure(2)
    # p.subplot(3,1,1)
    # p.plot(time, lc)
    # p.subplot(3,1,2)
    # p.plot (medfilt (lc,3))
    # p.subplot(3,1,3)
    # p.plot (medfilt (lc,5))


    return


def medfilt (x, k):
    """Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating endpoints.
    """
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = numpy.zeros ((len (x), k), dtype=x.dtype)
    y[:,k2] = x
    for i in range (k2):
        j = k2 - i
        y[j:,i] = x[:-j]
        y[:j,i] = x[0]
        y[:-j,-(i+1)] = x[j:]
        y[-j:,-(i+1)] = x[-1]
    return numpy.median (y, axis=1)
 
 
if __name__ == '__main__':
    test ()
