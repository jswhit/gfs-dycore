from pyspharm import Spharmt, ncepsigma
import numpy as np
import matplotlib.pyplot as plt

ncount = 0
hrs = range(2400,9600,24); nhrs = len(hrs)
sigfile = None
for hr in hrs:
    filename = 'SIG.F%03i' % hr
    print filename
    if sigfile is None:
        sigfile = ncepsigma(filename)
    else:
        sigfile.filename = filename
    vrtspec,divspec,tempspec,zspec,lnpsspec,qspec = sigfile.specdata()
    if not ncount:
        vrtspectm = vrtspec/nhrs
        divspectm = divspec/nhrs
        tempspectm = tempspec/nhrs
    else:
        vrtspectm = vrtspectm + vrtspec/nhrs
        divspectm = vrtspectm + divspec/nhrs
        tempspectm = tempspectm + tempspec/nhrs
    ncount = ncount + 1

lons,lats = np.meshgrid(sigfile.lons,sigfile.lats)
psg = sigfile.spectogrd(lnpsspec)
psg = 10.*np.exp(psg) # hPa
print psg.min(), psg.max(), psg.shape
ug = np.empty((sigfile.nlevs,sigfile.nlats,sigfile.nlons),np.float32)
vg = np.empty((sigfile.nlevs,sigfile.nlats,sigfile.nlons),np.float32)
temp = np.empty((sigfile.nlevs,sigfile.nlats,sigfile.nlons),np.float32)
for k in range(sigfile.nlevs):
    ug[k],vg[k] = sigfile.getuv(vrtspectm[k],divspectm[k])
    temp[k] = sigfile.spectogrd(tempspectm[k])
    spd = np.sqrt(ug[k]**2 + vg[k]**2)
    print k,spd.max(),temp[k].min(),temp[k].max()

ugzm = ug.mean(axis=2)
vgzm = vg.mean(axis=2)
tempzm = temp.mean(axis=2)

clevs = np.arange(-50,51,5)
print ugzm.min(), ugzm.max()
print vgzm.min(), vgzm.max()
print tempzm.min(), tempzm.max()
CS=plt.contour(lats[:,0],np.arange(1,sigfile.nlevs+1),ugzm,clevs,colors='k',linewidths=0.5)
plt.contour(lats[:,0],np.arange(1,sigfile.nlevs+1),ugzm,[0],colors='k',linewidths=1.5)
plt.clabel(CS,fmt='%i')
plt.contourf(lats[:,0],np.arange(1,sigfile.nlevs+1),ugzm,clevs,cmap=plt.cm.RdBu_r,extend='both')
plt.colorbar()
plt.xlabel('latitude')
plt.ylabel('model level')
plt.title('zonal wind')

plt.figure()
clevs = np.arange(186,311,4)
CS=plt.contour(lats[:,0],np.arange(1,sigfile.nlevs+1),tempzm,clevs,colors='k',linewidths=0.5)
plt.clabel(CS,fmt='%i')
plt.contourf(lats[:,0],np.arange(1,sigfile.nlevs+1),tempzm,clevs,cmap=plt.cm.RdBu_r,extend='both')
plt.colorbar()
plt.xlabel('latitude')
plt.ylabel('model level')
plt.title('temperature')
plt.show()
