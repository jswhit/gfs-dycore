import shtns, numpy

class Spharmt(object):
    """
    wrapper class for commonly used spectral transform operations in
    atmospheric models.  Provides an interface to shtns compatible
    with pyspharm (pyspharm.googlecode.com).

    Jeffrey S. Whitaker <jeffrey.s.whitaker@noaa.gov>
    """
    def __init__(self,nlons,nlats,ntrunc,rsphere,gridtype='gaussian'):
        """initialize
        nlons:  number of longitudes
        nlats:  number of latitudes"""
        self._shtns = shtns.sht(ntrunc, ntrunc, 1, \
                                shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)
        if gridtype == 'gaussian':
            self._shtns.set_grid(nlats,nlons,shtns.sht_gauss_fly|shtns.SHT_PHI_CONTIGUOUS,1.e-10)
        elif gridtype == 'regular':
            self._shtns.set_grid(nlats,nlons,shtns.sht_reg_dct|shtns.SHT_PHI_CONTIGUOUS,1.e-10)
        self.lats = numpy.arcsin(self._shtns.cos_theta)
        self.lons = (2.*numpy.pi/nlons)*numpy.arange(nlons) 
        self.nlons = nlons
        self.nlats = nlats
        self.ntrunc = ntrunc
        self.nlm = self._shtns.nlm
        self.degree = self._shtns.l
        self.lap = -self.degree*(self.degree+1.0).astype(numpy.complex)
        self.invlap = numpy.zeros(self.lap.shape, self.lap.dtype)
        self.invlap[1:] = 1./self.lap[1:]
        self.rsphere = rsphere
        self.lap = self.lap/rsphere**2
        self.invlap = self.invlap*rsphere**2
    def grdtospec(self,data):
        """compute spectral coefficients from gridded data"""
        data = numpy.ascontiguousarray(data, dtype=numpy.float)
        dataspec = numpy.empty(self.nlm, dtype=numpy.complex)
        self._shtns.spat_to_SH(data, dataspec)
        return dataspec
    def spectogrd(self,dataspec):
        """compute gridded data from spectral coefficients"""
        dataspec = numpy.ascontiguousarray(dataspec, dtype=numpy.complex)
        data = numpy.empty((self.nlats,self.nlons), dtype=numpy.float)
        self._shtns.SH_to_spat(dataspec, data)
        return data
    def getuv(self,vrtspec,divspec):
        """compute wind vector from spectral coeffs of vorticity and divergence"""
        vrtspec = numpy.ascontiguousarray(vrtspec, dtype=numpy.complex)
        divspec = numpy.ascontiguousarray(divspec, dtype=numpy.complex)
        u = numpy.empty((self.nlats,self.nlons), dtype=numpy.float)
        v = numpy.empty((self.nlats,self.nlons), dtype=numpy.float)
        self._shtns.SHsphtor_to_spat((self.invlap/self.rsphere)*vrtspec,\
               (self.invlap/self.rsphere)*divspec, u, v)
        return u,v
    def getvrtdivspec(self,u,v):
        """compute spectral coeffs of vorticity and divergence from wind vector"""
        u = numpy.ascontiguousarray(u, dtype=numpy.float)
        v = numpy.ascontiguousarray(v, dtype=numpy.float)
        vrtspec = numpy.empty(self.nlm, dtype=numpy.complex)
        divspec = numpy.empty(self.nlm, dtype=numpy.complex)
        self._shtns.spat_to_SHsphtor(u, v, vrtspec, divspec)
        return self.lap*self.rsphere*vrtspec, self.lap*rsphere*divspec
    def getgrad(self,divspec):
        """compute gradient vector from spectral coeffs"""
        divspec = numpy.ascontiguousarray(divspec, dtype=numpy.complex)
        vrtspec = numpy.zeros(divspec.shape, dtype=numpy.complex)
        u = numpy.empty((self.nlats,self.nlons), dtype=numpy.float)
        v = numpy.empty((self.nlats,self.nlons), dtype=numpy.float)
        self._shtns.SHsphtor_to_spat(vrtspec,divspec)
        return u/rsphere,v/rsphere

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap, addcyclic
    import time

    # non-linear barotropically unstable shallow water test case
    # of Galewsky et al (2004, Tellus, 56A, 429-440).
    # "An initial-value problem for testing numerical models of the global
    # shallow-water equations"
    # http://www-vortex.mcs.st-and.ac.uk/~rks/reprints/galewsky_etal_tellus_2004.pdf

    # grid, time step info
    nlons = 256  # number of longitudes
    ntrunc = nlons/3  # spectral truncation (for alias-free computations)
    nlats = nlons/2 # for gaussian grid.
    dt = 150 # time step in seconds
    itmax = 6*(86400/dt) # integration length in days

    # parameters for test
    rsphere = 6.37122e6 # earth radius
    omega = 7.292e-5 # rotation rate
    grav = 9.80616 # gravity
    hbar = 10.e3 # resting depth
    umax = 80. # jet speed
    phi0 = numpy.pi/7.; phi1 = 0.5*numpy.pi - phi0; phi2 = 0.25*numpy.pi
    en = numpy.exp(-4.0/(phi1-phi0)**2)
    alpha = 1./3.; beta = 1./15.
    hamp = 120. # amplitude of height perturbation to zonal jet
    efold = 3.*3600. # efolding timescale at ntrunc for hyperdiffusion
    ndiss = 8 # order for hyperdiffusion

    # setup up spherical harmonic instance, set lats/lons of grid
    x = Spharmt(nlons,nlats,ntrunc,rsphere,gridtype='gaussian')
    lons,lats = numpy.meshgrid(x.lons, x.lats)
    f = 2.*omega*numpy.sin(lats) # coriolis

    # zonal jet.
    vg = numpy.zeros((nlats,nlons),numpy.float)
    u1 = (umax/en)*numpy.exp(1./((x.lats-phi0)*(x.lats-phi1)))
    ug = numpy.zeros((nlats),numpy.float)
    ug = numpy.where(numpy.logical_and(x.lats < phi1, x.lats > phi0), u1, ug)
    ug.shape = (nlats,1)
    ug = ug*numpy.ones((nlats,nlons),dtype=numpy.float) # broadcast to shape (nlats,nlonss)
    # height perturbation.
    hbump = hamp*numpy.cos(lats)*numpy.exp(-(lons/alpha)**2)*numpy.exp(-(phi2-lats)**2/beta)

    # initial vorticity, divergence in spectral space
    vrtspec, divspec =  x.getvrtdivspec(ug,vg)
    vrtg = x.spectogrd(vrtspec)
    divg = x.spectogrd(divspec)

    # create spectral indexing arrays, laplacian operator and its inverse.
    hyperdiff_fact = numpy.exp((-dt/efold)*(x.lap/x.lap[-1])**(ndiss/2))

    # solve nonlinear balance eqn to get initial zonal geopotential,
    # add localized bump (not balanced).
    vrtg = x.spectogrd(vrtspec)
    tmpg1 = ug*(vrtg+f); tmpg2 = vg*(vrtg+f)
    tmpspec1, tmpspec2 = x.getvrtdivspec(tmpg1,tmpg2)
    tmpspec2 = x.grdtospec(0.5*(ug**2+vg**2))
    phispec = x.invlap*tmpspec1 - tmpspec2
    phig = grav*(hbar + hbump) + x.spectogrd(phispec)
    phispec = x.grdtospec(phig)

    # initialize spectral tendency arrays
    ddivdtspec = numpy.zeros(vrtspec.shape+(3,), numpy.complex)
    dvrtdtspec = numpy.zeros(vrtspec.shape+(3,), numpy.complex)
    dphidtspec = numpy.zeros(vrtspec.shape+(3,), numpy.complex)
    nnew = 0; nnow = 1; nold = 2

    # time loop.
    time1 = time.clock()
    for ncycle in range(itmax+1):
        t = ncycle*dt
        # get vort,u,v,phi on grid
        vrtg = x.spectogrd(vrtspec)
        ug,vg = x.getuv(vrtspec,divspec)
        phig = x.spectogrd(phispec)
        print 't=%6.2f hours: min/max %6.2f, %6.2f' % (t/3600.,vg.min(), vg.max())
        # compute tendencies.
        tmpg1 = ug*(vrtg+f); tmpg2 = vg*(vrtg+f)
        ddivdtspec[:,nnew], dvrtdtspec[:,nnew] = x.getvrtdivspec(tmpg1,tmpg2)
        dvrtdtspec[:,nnew] *= -1
        tmpg = x.spectogrd(ddivdtspec[:,nnew])
        tmpg1 = ug*phig; tmpg2 = vg*phig
        tmpspec, dphidtspec[:,nnew] = x.getvrtdivspec(tmpg1,tmpg2)
        dphidtspec[:,nnew] *= -1
        tmpspec = x.grdtospec(phig+0.5*(ug**2+vg**2))
        ddivdtspec[:,nnew] += -x.lap*tmpspec
        # update vort,div,phiv with third-order adams-bashforth.
        # forward euler, then 2nd-order adams-bashforth time steps to start.
        if ncycle == 0:
            dvrtdtspec[:,nnow] = dvrtdtspec[:,nnew]
            dvrtdtspec[:,nold] = dvrtdtspec[:,nnew]
            ddivdtspec[:,nnow] = ddivdtspec[:,nnew]
            ddivdtspec[:,nold] = ddivdtspec[:,nnew]
            dphidtspec[:,nnow] = dphidtspec[:,nnew]
            dphidtspec[:,nold] = dphidtspec[:,nnew]
        elif ncycle == 1:
            dvrtdtspec[:,nold] = dvrtdtspec[:,nnew]
            ddivdtspec[:,nold] = ddivdtspec[:,nnew]
            dphidtspec[:,nold] = dphidtspec[:,nnew]
        vrtspec += dt*( \
        (23./12.)*dvrtdtspec[:,nnew] - (16./12.)*dvrtdtspec[:,nnow]+ \
        (5./12.)*dvrtdtspec[:,nold] )
        divspec += dt*( \
        (23./12.)*ddivdtspec[:,nnew] - (16./12.)*ddivdtspec[:,nnow]+ \
        (5./12.)*ddivdtspec[:,nold] )
        phispec += dt*( \
        (23./12.)*dphidtspec[:,nnew] - (16./12.)*dphidtspec[:,nnow]+ \
        (5./12.)*dphidtspec[:,nold] )
        # implicit hyperdiffusion for vort and div.
        vrtspec *= hyperdiff_fact
        divspec *= hyperdiff_fact
        # switch indices, do next time step.
        nsav1 = nnew; nsav2 = nnow
        nnew = nold; nnow = nsav1; nold = nsav2

    time2 = time.clock()
    print 'CPU time = ',time2-time1

    # make a orthographic plot of potential vorticity.
    m = Basemap(projection='ortho',lat_0=45,lon_0=0)
    # dimensionless PV
    pvg = (0.5*hbar*grav/omega)*(vrtg+f)/phig
    print 'max/min PV',pvg.min(), pvg.max()
    lons1d = (180./numpy.pi)*x.lons; lats1d = (180./numpy.pi)*x.lats
    pvg,lons1d = addcyclic(pvg,lons1d)
    lons, lats = numpy.meshgrid(lons1d,lats1d)
    x,y = m(lons,lats)
    levs = numpy.arange(-0.2,1.801,0.1)
    m.drawmeridians(numpy.arange(-180,180,60))
    m.drawparallels(numpy.arange(-80,81,20))
    CS=m.contourf(x,y,pvg,levs,cmap=plt.cm.spectral,extend='both')
    m.colorbar()
    plt.title('PV (T%s with hyperdiffusion, hour %6.2f)' % (ntrunc,t/3600.))
    plt.show()
