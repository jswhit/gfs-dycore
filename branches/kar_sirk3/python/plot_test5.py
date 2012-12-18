from pyspharm import Spharmt, ncepsigma
import numpy as np
from mpl_toolkits.basemap import Basemap, addcyclic, cm
import matplotlib.pyplot as plt
import sys

m = Basemap(llcrnrlat=-0,urcrnrlat=40,llcrnrlon=140,urcrnrlon=220,resolution=None)

fig = plt.figure(figsize=(15,5))

npanel = 1
sigfile = None
nlevsout = 30

for fhour in [72,120,240]:

    filename = 'SIG.F%s' % fhour

    if sigfile is None:
        sigfile = ncepsigma(filename)
    else:
        sigfile.filename = filename
    vrtspec,divspec,tempspec,zspec,lnpsspec,qspec = sigfile.specdata()

    lons,lats = np.meshgrid(sigfile.lons,sigfile.lats)
    nlons = lons.shape[1]; nlats = lons.shape[0]; nlevs = vrtspec.shape[0]
    ug = np.empty((nlevs,nlats,nlons),np.float)
    vg = np.empty((nlevs,nlats,nlons),np.float)
    for k in range(nlevs):
        ug[k,:,:],vg[k,:,:] = sigfile.getuv(vrtspec[k],divspec[k])
    spd = np.sqrt(ug**2+vg**2)
    psg = sigfile.spectogrd(lnpsspec)
    psg = 10.*np.exp(psg) # hPa
    jmin,imin = np.unravel_index(np.argmin(psg),(nlats,nlons))
    print fhour,spd.max(),psg.min(), psg[jmin,imin]
    spdxsxn = spd[:nlevsout,jmin,imin-10:imin+10].squeeze()
    lonsxsxn = lons[jmin,imin-10:imin+10].squeeze()

    x,y = m(lons,lats)
    levs = np.arange(5,71,5)
    ax = fig.add_subplot(2,3,npanel)
    cs = m.contourf(x,y,spd[12],levs,cmap=cm.GMT_haxby_r,extend='both')
    m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=9)
    m.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1],fontsize=9)
    plt.title('day %s' % int(fhour/24.), y=1.1)

    ax = fig.add_subplot(2,3,npanel+3,aspect=0.15)
    print spdxsxn.shape, nlevs, lonsxsxn.shape
    print spdxsxn.min(), spdxsxn.max()
    ax.contourf(lonsxsxn,np.arange(1,nlevsout+1),spdxsxn,levs,cmap=cm.GMT_haxby_r,extend='both')
    ax.set_xlabel('longitude')
    ax.set_ylabel('model level')
    npanel = npanel + 1

# a single colorbar.
cax = plt.axes([0.95, 0.1, 0.02, 0.8])
plt.colorbar(cs, cax=cax, orientation='vertical')
plt.show()
