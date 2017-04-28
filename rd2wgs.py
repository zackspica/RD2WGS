#!/usr/bin/env python
import numpy as np


def rd2bessel(x, y): 
    """ Convert xy to Bessel """

    x0      = 1.55e5 
    y0      = 4.63e5 
    k       = 0.9999079 
    bigr    = 6382644.571 
    m       = 0.003773953832 
    n       = 1.00047585668 
    e       = 0.08169683122

    lambda0 = np.pi*0.029931327161111111 
    b0      = np.pi*0.28956165138333334

    d1 = x-x0
    d2 = y-y0
    r = np.sqrt(d1**2 + d2**2)

    sa = np.zeros( x.shape )
    ca = np.zeros( x.shape )
    
    ix = (np.abs(r) > 1e-10).nonzero()[0]
    
    sa[ix] = ( x[ix]-x0 ) / r[ix]
    ca[ix] = ( y[ix]-y0 ) / r[ix]

    psi = np.arctan2(r, k*2*bigr) * 2

    cpsi = np.cos(psi) 
    spsi = np.sin(psi)

    sb = ca*np.cos(b0)*spsi + np.sin(b0)*cpsi 
    d1 = sb 
    cb = np.sqrt(1- d1**2) 
    b = np.arccos(cb) 
    sdl = sa * spsi / cb 
    dl = np.arcsin(sdl) 
    lambd = dl / n + lambda0 
    w = np.log( np.tan(b / 2 + np.pi / 4) ) 
    q = (w - m) / n 
    phiprime = np.arctan(np.exp(q)) * 2 - np.pi / 2
    
    for i in xrange(4):    
        dq = e / 2 * np.log( ( e * np.sin(phiprime) + 1 ) / \
             ( 1 - e * np.sin(phiprime) ) ) 
        phi = np.arctan( np.exp( q + dq ) ) * 2 - np.pi / 2
        phiprime = phi 
    
    lambd = lambd/np.pi * 180
    phi   = phi/np.pi * 180

    return phi, lambd


def bessel2wgs84(phibes, lambes):
    """   
    Convert Bessel2 WGS84
    """
    a = 52
    b = 5 
    c = -96.862 
    d = 11.714 
    e = 0.125
    f = 1e-5 
    g = 0.329 
    h = 37.902 
    i = 14.667
    
    dphi = phibes-a
    dlam = lambes-b
    phicor = (c - dphi * d - dlam * e) * f 
    lamcor = (dphi *g - h - dlam * i) * f
    phiwgs = phibes + phicor
    lamwgs = lambes + lamcor

    return np.round(phiwgs,5), np.round(lamwgs,5)


def rd2wgs(x, y): 
    """
    main fct
    """
    if isinstance(x, int):
        x = np.array([x]); y = np.array([y])
    phibes, lambes = rd2bessel(x, y)
    phiwgs, lamwgs = bessel2wgs84(phibes,lambes)
    
    for idx in xrange(0,len(phiwgs),1):
        print 'N %s, E %s'%(phiwgs[idx], lamwgs[idx]) 

    return phiwgs, lamwgs 

if __name__=='__main__':
    x = 24220824
    y = 59273422
    rd2wgs(x,y)
    
#EOF
