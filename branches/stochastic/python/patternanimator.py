import numpy as np
from spharm import gaussian_lats_wts
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation

nlons = 384
nlats = 192
rsphere = 6.3712e6
lats1d, wts = gaussian_lats_wts(nlats)
lons1d = (360./nlons)*np.arange(nlons)
dt = 1200.
ntmax = 5000


filename = 'patterns.dat'
f = open(filename,'rb')
psig = np.fromstring(f.read(nlons*nlats*4),'<f4').reshape(nlats,nlons)

fig = plt.figure()
m = Basemap(projection='ortho',lon_0=-105,lat_0=40)
m.drawcoastlines()
lons, lats = np.meshgrid(lons1d,lats1d)
x, y = m(lons,lats)
levs = np.linspace(-1.,1,21)
txt = plt.title('t = 0')
CS = m.contourf(x,y,psig,levs)
f.seek(0)

nt = 0
psigmean = np.zeros(lons.shape, np.float32)
def updatefig(*args):
    global CS,nt,txt,psigmean
    t = nt*dt/3600.
    psig = np.fromstring(f.read(nlons*nlats*4),'<f4').reshape(nlats,nlons)
    psigmean = psigmean + psig
    print nt,(psigmean/(nt+1)).min(),(psigmean/(nt+1)).max(),(psigmean/(nt+1)).mean()
    print psig.min(), psig.max()
    for c in CS.collections: c.remove()
    CS=m.contourf(x,y,psig,levs)
    txt.set_text('t = %6.2f' % t)
    if nt == ntmax:
        f.seek(0)
        nt = -1
        psigmean = np.zeros(lons.shape, np.float32)
    nt += 1

ani = animation.FuncAnimation(fig, updatefig)
plt.show()
