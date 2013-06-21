#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import numpy
import scipy.ndimage

#this shit only works with DFM's injected LCs

def find_da_gap_bitch(lc,pidx):
    gaps = []
    gap = []

    #pr
    #print lc[:,pidx]
    #flux = lc[:,pidx]
    flux = lc
    
    for i in range(len(lc)-1):
        if numpy.isnan(flux[i]):
            if numpy.isnan(flux[i+1]):
                gap.append(i)
            else:
                gap.append(i)
                gaps.append(gap)
                gap = []
                
    return [(len(i),i) for i in gaps]

def fill_da_gap_bitch(v1,v2,n):
    
    #v1 - flux value before the gap starts
    #v2 - flux value after the gap ends
    #n  - number of points in the gap
    
    x = numpy.linspace(0,1,n+2)
    y = (v2-v1)*x + v1
    
    return y[1:-1]
    
def slice_da_bitch_up(lc,pidx,hsize=10):
    gaps = find_da_gap_bitch(lc,pidx)
    lc_chunks = []
    flux = lc#[:,pidx]
    
    for g in gaps:
        
        gs = g[1][0] #gap starting index
        ge = g[1][-1] #gap starting index  
        
        if g[0]<hsize: #fill the gap if it's smaller than hsize points
            flux[numpy.array(g[1])] = fill_da_gap_bitch(flux[gs-1],
                                                        flux[ge+1],g[0])
        else:
            lc_chunks.append(gs)
            lc_chunks.append(ge)

    lc_chunks.insert(0,0)
    lc_chunks.append(len(lc))

    return lc_chunks

def detrend_using_mf(lc,flavor='SAP',size=14,mode='constant'):
    
    # either use a gap (size is int):
    #
    #             size
    #<.......--------------.........................................>
    #_/▔﹀\_︿╱﹀╲/╲︿_/︺╲▁︹_/﹀\_︿╱▔︺\/\︹▁╱﹀▔╲︿_/︺▔╲▁︹_/﹀▔\⁄﹀\╱
    #
    #
    #
    # or a gap with with a hole (size is 2-element array):
    #
    #        size[0] size[1] size[0]
    #<.......--------_______--------................................>
    #_/▔﹀\_︿╱﹀╲/╲︿_/︺╲▁︹_/﹀\_︿╱▔︺\/\︹▁╱﹀▔╲︿_/︺▔╲▁︹_/﹀▔\⁄﹀\╱
    
    pidx={'SAP':1,'PDC':3}[flavor]
    
    lc_chunks = slice_da_bitch_up(lc,pidx)
    dtlc = numpy.array(lc)#[:,pidx])
    
    for i in range(0,len(lc_chunks),2):
        
        try:
            len(size)
            
            sidesize = size[0]
            gapsize = size[1]
            
            footprint = [True]*sidesize + [False]*gapsize + [True]*sidesize
            
            dt = scipy.ndimage.filters.median_filter(
            #lc[:,pidx][lc_chunks[i]:lc_chunks[i+1]],footprint=footprint,
            lc[lc_chunks[i]:lc_chunks[i+1]],footprint=footprint,
            mode=mode,cval=1.0)
            
        except TypeError:             
            dt = scipy.ndimage.filters.median_filter(
            #lc[:,pidx][lc_chunks[i]:lc_chunks[i+1]],size=size,
            lc[lc_chunks[i]:lc_chunks[i+1]],size=size,
            mode=mode,cval=1.0)

            
        #dtlc[lc_chunks[i]:lc_chunks[i+1]] = lc[:,pidx][lc_chunks[i]:lc_chunks[i+1]] / dt
        dtlc[lc_chunks[i]:lc_chunks[i+1]] = lc[lc_chunks[i]:lc_chunks[i+1]] / dt
    return dtlc


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    fn='/Users/angusr/angusr/data2/SAMSI/kplr000892203-2009131105131_llc.txt'
    lc = numpy.loadtxt(fn,delimiter=',')

    dt = detrend_using_mf(lc,size=15)
    dth = detrend_using_mf(lc,size=(10,10))

    plt.plot(lc[:,0],dt)
    plt.plot(lc[:,0],dth+0.005)

    dt = detrend_using_mf(lc,flavor='PDC',size=15)
    dth = detrend_using_mf(lc,flavor='PDC',size=(10,10))

    plt.plot(lc[:,0],dt+0.015)
    plt.plot(lc[:,0],dth+0.02)


    plt.show()
