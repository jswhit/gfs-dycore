from read_sigma import read_header, read_specdata
import spharm
import numpy as np
import matplotlib.pyplot as plt

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
        #return (vrtspec.T).astype(np.complex64),(divspec.T).astype(np.complex64),(tempspec.T).astype(np.complex64),zspec.astype(np.complex64),lnpsspec.astype(np.complex64),(qspec.T).astype(np.complex64)

ncount = 0
hrs = range(2400,9997,12); nhrs = len(hrs)
for hr in hrs:
    filename = 'sig.f%04i' % hr
    print filename
    sigfile = ncepsigma(filename)
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
clevs = np.arange(200,311,4)
CS=plt.contour(lats[:,0],np.arange(1,sigfile.nlevs+1),tempzm,clevs,colors='k',linewidths=0.5)
plt.clabel(CS,fmt='%i')
plt.contourf(lats[:,0],np.arange(1,sigfile.nlevs+1),tempzm,clevs,cmap=plt.cm.RdBu_r,extend='both')
plt.colorbar()
plt.xlabel('latitude')
plt.ylabel('model level')
plt.title('temperature')
plt.show()
