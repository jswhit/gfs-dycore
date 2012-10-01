from read_sigma import read_header, read_specdata, read_griddata
import spharm
import numpy as np
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.pyplot as plt
import sys

class ncepsigma(object):
    def __init__(self,filename):
        nlons,nlats,nlevs,ntrunc = read_header(filename)
        self.nlons = nlons; self.nlats = nlats
        self.ntrunc = ntrunc; self.nlevs = nlevs
        self.filename = filename
        self.sp = spharm.Spharmt(nlons,nlats,6.3712e6,gridtype='gaussian')
        lats,wts = spharm.gaussian_lats_wts(nlats)
        self.lats = lats
        self.lons = (360./nlons)*np.arange(nlons)
    def spectogrd(self,specdata):
        return self.sp.spectogrd(specdata)
    def getuv(self,vrtdata,divdata):
        return self.sp.getuv(vrtdata,divdata)
    def specdata(self):
        vrtspec, divspec,tempspec,zspec,lnpsspec,qspec =\
        read_specdata(self.filename,self.ntrunc,self.nlevs)
        return vrtspec.T,divspec.T,tempspec.T,zspec,lnpsspec,qspec.T

#filename = sys.argv[1]
#m = Basemap(projection='npstere',boundinglat=15,lon_0=90,round=True)
m = Basemap(llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution=None)
fig = plt.figure(figsize=(15,10))
npanel = 1
#for fhour in [168,216,288,360]:
for fhour in [96,144,192,240]:
    filename = 'SIG.F%s' % fhour

    sigfile = ncepsigma(filename)
    vrtspec,divspec,tempspec,zspec,lnpsspec,qspec = sigfile.specdata()
    
    lons,lats = np.meshgrid(sigfile.lons,sigfile.lats)
    psg = sigfile.spectogrd(lnpsspec)
    psg = 10.*np.exp(psg) # hPa
    print psg.min(), psg.max(), psg.shape
    
    # add wraparound (cyclic) points.
    lons1 = lons[0,:]
    lats1 = lats[:,0]
    psg, lons1 = addcyclic(psg, lons1)
    lons, lats = np.meshgrid(lons1,lats1)
    x,y = m(lons,lats)
    levs = np.arange(900,1050,10).tolist()
    #levs.remove(1000)
    ax = fig.add_subplot(2,2,npanel)
    #m.contour(x,y,psg,[1000.],linestyles='dotted',linewidths=1,colors='k')
    cs1 = m.contour(x,y,psg,levs,linewidths=1,colors='k')
    cs2 = m.contourf(x,y,psg,levs,cmap=plt.cm.jet,extend='both')
    m.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    plt.title('day %s' % int(fhour/24.), y=1.1)
    npanel = npanel + 1

# a single colorbar.
cax = plt.axes([0.1, 0.05, 0.8, 0.025])
plt.colorbar(cs2, cax=cax, orientation='horizontal')
plt.show()
