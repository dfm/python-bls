# python BLS script

import numpy
import bls


def BLS(time, lc, df = 0.0001, nf = 50000,  nb = 200, qmi = 0.01,\
        qma = 0.8, fmin = (1./(20.0*1.1))):

    u = numpy.ones(len(time))
    v = numpy.ones(len(time))

    BLS = bls.eebls(time, lc, u, v, nf, fmin, df, nb, qmi, qma)
    f = fmin + (numpy.arange(len(BLS[0])))*df
    
    return BLS, 1/f


if __name__ == '__main__':
    BLS()
