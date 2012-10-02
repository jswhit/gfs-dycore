from pyspharm import Spharmt, ncepsigma
import numpy as np
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.pyplot as plt
import sys

m = Basemap(llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution=None)

fig = plt.figure(figsize=(15,10))

npanel = 1
sigfile = None

for fhour in [168,216,288,360]:

    filename = 'SIG.F%s' % fhour

    if sigfile is None:
        sigfile = ncepsigma(filename)
    else:
        sigfile.filename = filename
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
