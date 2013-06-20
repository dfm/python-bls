import numpy
from scipy.signal import signaltools
from scipy import interpolate
import pylab
import scipy

def medsig(array):
    med = numpy.median(array)
    sig = 1.48 * numpy.median(abs(array - med))
    return med, sig

def filt1d(array, nmed, nlin, fill = False, circ = False):
    """Nonlinear (median+boxcar) filter with edge reflection."""
    nmed = 2 * int(nmed/2) + 1
    N = len(array)
    if N < (3*nmed): 
        return array
    # Check for NaNs
    lnan = numpy.isnan(array)
    nnan = sum(lnan)
    if (nnan != 0):
        lnotnan = numpy.where(lnan==0)
        s = numpy.shape(lnotnan)
        if (len(s)>1):
            s = numpy.size(lnotnan)
            lnotnan = numpy.reshape(lnotnan, s)
            med = numpy.median(array[lnotnan])
        # Fill in any data gaps by interpolating over them ...
        work = numpy.zeros(N)
        il = numpy.min(lnotnan)
        ih = numpy.max(lnotnan)
        xnew = numpy.arange(ih-il+1) + il
        f = interpolate.interp1d(lnotnan, array[lnotnan])  
        work[il:ih+1] = f(xnew)
        # ... or extending slope of nearest data points if at the edges
        if (il!=0):
            slope = work[il+1] - work[il]
            for i in range(il): work[i] = work[il] - slope*(il-i)
        if (ih!=N-1):
            slope = work[ih] - work[ih-1]
            for i in range(N-ih-1)+ih+1: work[i] = work[ih] + slope*(i-ih)
    else:
        work = numpy.copy(array)
    # Edge reflection
    nr = min(20, nlin)
    sz = max([nmed, nlin])
    if sz >= (N-1): 
        sz = N-2
    if circ != False:
        wl = work[N-sz:]
        wr = work[:sz]
    else:
        wl = array[0:sz]
        pivot = numpy.median(array[:nr])
        wl = 2 * pivot - wl
        wl = wl[::-1] # reverse array
        wr = array[N-sz:N]
        pivot = numpy.median(array[N-nr:N])
        wr = 2 * pivot - wr
        wr = wr[::-1] # reverse array
    work2 = numpy.zeros(N + 2 * sz)
    work2[0:sz] = wl
    work2[sz:sz+N] = work
    work2[sz+N:2*sz+N] = wr
    # Filter
    if nmed > 1: work = signaltools.medfilt(work2, nmed) 
    else: work = work2
    box = scipy.signal.boxcar(2*nlin+1)
    work = signaltools.convolve(work, box) / float(2*nlin+1)
    padd = (len(work) - N - 2 * sz) / 2
    work = work[padd:padd+N+2*sz]
    # return to orginal array bounds
    result = work[sz:N+sz]
    # replace bad data if present
    if (fill==False) * (nnan!=0): result[lnan] = numpy.nan
    return result

def smoothe(array, box, fill = False, circ = False):
    """Boxcar filter with edge reflection"""
    nlin = int(box / 2)
    return filt1d(array, 1, nlin, fill = fill, circ = circ)

def NIF(array, nmed, nlin, nsig = 3.0, prefilt = False, \
            fill = False, circ = False):
    """Non-linear iterative filter with k-sigma clipping (Aigrain &
    Irwin 2004)."""
    # Pre-filter if requested
    if prefilt != False:
        wd = filt1d(array[:], 7, 3) 
    else:
        wd = array[:]
    wd2 = wd
    # Start iteration loop
    irejo = 0
    irej = -1
    for iter in range(5):
        if iter > 0:
            irejo = irej
            out = numpy.isnan(wd) + (abs(wd2-wd) > (nsig*sigma))
            irej = sum(out)
            wd = array[:]
            print iter, irej, sigma
            if irej != 0: wd[out] = numpy.nan
        wd = filt1d(wd, nmed, nlin, fill = fill, circ = circ)
        sigma = 1.48 * numpy.median(abs(wd-array))
        if irej <= irejo: break
    return wd

def IRF(period, time, flux, timescale = 0.99, nbin = 3000):
    """Iterative reconstruction filter (Alapini & Aigrain 2009)."""
    nlin = 5
    fstar = numpy.copy(flux)
    ndata = len(time)
    tstep = numpy.median(time[1:]-time[:ndata-1])
    nmed = int(timescale / tstep)
    nsmooth = int(ndata / nbin)
    cvlim = 1.e-4

    # Normalise LC
    medstar, sigstar = medsig(fstar)
    fnorm = fstar / medstar
    signorm = sigstar / medstar

    # Fold light curve
    tmod = (time % period)
    t_pf = tmod / period
    spf = numpy.argsort(t_pf)
    st = numpy.argsort(spf)
  
    # Initial estimate of stellar variability
    fvar = fnorm - fnorm
    fvar_im1 = fvar

    # Initial estimate of transit signal
    trans = fnorm - fvar

    # Initial residuals
    fres = fnorm / trans - fvar
    sigres = signorm
  
    # Start iteration
    i = -1
    conv = False
    while conv == False:

        # Store results of previous iteration
        fvar_im2   = fvar_im1
        fvar_im1   = fvar
        trans_old  = trans
        fres_old   = fres
        sigres_old = sigres
        i += 1  

        print i
        pylab.clf()
        pylab.subplot(211)
        pylab.plot(t_pf, fnorm-fvar, 'k.')

        # Estimate transit signal
        transf = smoothe(fnorm[spf]-fvar[spf], nsmooth, circ = True)
        trans  = transf[st]
        pylab.plot(t_pf, trans, 'r.')

        # Estimate transit-less signal
        fnotrans = fnorm / trans
        mednotrans, signotrans = medsig(fnotrans)

        # Estimate stellar signal
        fvar = NIF(fnotrans, nmed, nlin)
        pylab.subplot(212)
        pylab.plot(time, fnorm, 'k-')
        pylab.plot(time, fvar, 'b-')
        fvar -= numpy.median(fvar)

        # Evaluate residuals and scatter
        fres = fnorm / trans - fvar
        medres, sigres = medsig(fres)

        # Test for convergence
        stcur = numpy.sum((fres-medres)**2) / float(ndata - 1)        
        print stcur
        raw_input('?')

        if i == 0: stat = stcur
        else: stat = numpy.append(stat, stcur)
        if i >= 2:
            dstat = abs(stat[1:] - stat[:i])
            if (dstat > cvlim).any() == False: conv = 1

    # Final arrays
    fvar = fvar_im2
    transf = smoothe(fnorm[spf] - fvar[spf], nsmooth, circ = True)
    trans  = transf[st]
    fres = fnorm / trans - fvar
    fvar *= medstar
    fres *= medstar

    return trans, fvar, fres

def bin_old(x, y, binsz = None, nbin = None):
    xmin = min(x)
    xmax = max(x)
    xrange = xmax - xmin
    if nbin == None:
        if binsz == None: 
            print 'bin: must supply either binsz or nbin'
        nbin = int(xrange / binsz)
    if nbin <= 1: return x, y
    binsz = xrange / float(nbin)
    xbin = scipy.r_[xmin:xmax-binsz:nbin*1j] + binsz / 2.
    ybin = scipy.zeros(nbin) + scipy.nan
    for i in scipy.arange(nbin):
        l = (abs(x-xbin[i]) <= binsz / 2.)
        n = scipy.sum(l)
        if n == 0: continue
        ybin[i] = scipy.mean(y[l])
    return xbin, ybin


def bin(x, y, binsz = None, nbin = None):
    xmin = min(x)
    xmax = max(x)
    xrange = xmax - xmin
    if nbin == None:
        if binsz == None: 
            print 'bin: must supply either binsz or nbin'
        nbin = int(xrange / binsz)
    if nbin <= 1: return x, y
    binsz = xrange / float(nbin)
    
    xbin = scipy.r_[xmin:xmax-binsz:nbin*1j] + binsz / 2.
    xbinedge = xbin - (binsz / 2.)
    xbinedge = scipy.append(xbinedge, xbinedge.max()+(binsz))
    ybin = scipy.zeros(nbin) + scipy.nan
    indices = numpy.digitize(x, xbinedge)
    for i in scipy.arange(len(xbin)):
        ybin[i] = y[indices == i+1].mean()
        
    return xbin, ybin

