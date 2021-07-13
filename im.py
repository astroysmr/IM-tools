import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from scipy import constants
from scipy.interpolate import interp1d
import sys
import os
args = sys.argv
aa = 1.0
bb = (0.8/0.8)**-1
'''
PCII_mean = np.loadtxt('PCII.z5.mean.csv',delimiter=',')
kk = PCII_mean[:,0]
PP = PCII_mean[:,1]
plt.plot(kk,PP*2.*np.pi**2/kk**2,lw=4)
plt.xscale('log')
plt.yscale('log')
plt.show()
'''

# params
nu_line = 1900.53690000 * u.GHz# GHz
cc   = constants.c*u.m/u.s # m
zz   = 4.6
DD   = 50.*u.m # m
beam = (1.2*cc/(nu_line/(1.+zz))/DD).to('')/np.pi*180.*3600. * u.arcsec # arcsec
RR   = 500.
Wscan = float(args[2])*u.deg # deg
#dnu  = 0.6*u.GHz # GHz
dnu = nu_line/(1.+zz)/RR #GHz
Wnu  = 24.5*u.GHz # GHz

dz = dnu/nu_line*(1.+zz)**2
Totalobstime = float(args[1]) *u.hour
case = args[3]
if case == 'opt':
    color = 'red'
    casename = 'optimistic'
elif case == 'mean':
    color='green'
    casename = 'mean'
else:
    color='blue'
    casename = 'pesimistic'

tobs = Totalobstime/(Wscan**2/beam**2).to('')# hour
print(tobs)

# main
cosmo = FlatLambdaCDM(H0=100.,Om0=0.3)
PCII_mean = np.loadtxt('PCII.z5.'+case+'.csv',delimiter=',')

Vpix = (cosmo.kpc_comoving_per_arcmin(zz) * beam)**2 * cc/cosmo.H(zz) * dz
Wz = ((Wnu/dnu)*dz).to('').value
kmax_s = 1./(cosmo.kpc_comoving_per_arcmin(zz) * beam)
kmax_f = 1./(cc/cosmo.H(zz) * dz)
Vs = Vpix * (Wscan/beam)**2 * (Wnu/dnu)

#sigma_pix = 8. * u.MJy*u.s**0.5/u.sr
sigma_pix = (1./1.37)**0.5 * u.mJy * u.h**0.5 / beam

Pnoise = Vpix*sigma_pix**2/tobs

'''
def NN(kk_,delta_k_):
    kk = kk_/u.Mpc
    delta_k = delta_k_/u.Mpc
    return (np.pi*kk*delta_k*Vs/(2.*np.pi)**2).to('').value
'''

def NN(kk_,delta_k_):
    kk = kk_/u.Mpc
    delta_k = delta_k_/u.Mpc
    return (2.*np.pi*kk**2*delta_k*Vs/(2.*np.pi)**3).to('').value


sigma_pix = bb*1.3*u.mJy * (Totalobstime/(1000.*u.hour))**-0.5 * (Wscan**2/(1.*u.deg**2))**0.5 / beam**2
#sigma_pix = (1./1.37)**0.5 * u.mJy * u.h**0.5 / beam

Pnoise = Vpix*sigma_pix**2
print((sigma_pix.to('Jy / sr'))/3600.**0.5)
Pnoise = Pnoise.to('Mpc3 Jy2 / sr2').value

kk = 10.**np.linspace(-2,1,15)
PP = interp1d(np.log10(PCII_mean[:,0]), np.log10(PCII_mean[:,1]), kind='cubic',fill_value='extrapolate')(np.log10(kk))
PP = 10.**PP

print(Pnoise/1.0e9)
print((Wscan/beam).to('')**2)

plt.close()
plt.plot(kk,PP*aa,lw=4,color=color)

kmax = 2.*np.pi*np.max(np.array([kmax_s.to('/Mpc').value,kmax_f.to('/Mpc').value]))
kk = kk[kk<kmax]
dk_ = np.diff(kk)
SN = []
for i in range(dk_.shape[0]-1):
    PP_k = aa*10**(interp1d(np.log10(PCII_mean[:,0]), np.log10(PCII_mean[:,1]), kind='cubic',fill_value='extrapolate')(np.log10(kk[i+1])))
    dk = (dk_[i]+dk_[i+1])/2.
    N = NN(kk[i+1],dk)
    print('n:',N)
    if N>1:
        err = (PP_k+kk[i+1]**3*Pnoise/2./np.pi**3)/np.sqrt(N)
        #err = (PP_k+Pnoise)
        plt.errorbar([kk[i+1]],[PP_k],yerr=err,fmt='ko')

        SN.append((PP_k/err)**2)
        print('n:',N, ' err:', N*(PP_k/err))

print(kmax_s.to('/Mpc'))
print(kmax_f.to('/Mpc'))
print((np.array(SN).sum())**0.5)

plt.xscale('log')
plt.yscale('log')
plt.xlim(10**-2,10)
plt.ylim(1.0e2,1.0e11)
titlename = r'(t$_{on}$='+str(int(Totalobstime.value))+r'h, $S_{area}=$'+str((Wscan**2).value)+r'deg$^2$, $z$='+'{0:1.1f}'.format(zz)+r', $\Delta z$='+'{0:1.2f}'.format(Wz)+') ==> SNR='+'{0:1.2f}'.format((np.array(SN[0:8]).sum())**0.5)
plt.title('[CII] intensity mapping w/ LMT PWV=2mm: '+casename+' case\n'+titlename)
plt.xlabel(r'$k$ [$h$/Mpc]')
plt.ylabel(r'$k^3P_\mathrm{line} (k)/2\pi^3$ [(Jy/sr)$^2$]')
os.system('mkdir -p figs')
plt.savefig('figs/fig.'+str(int(Totalobstime.value))+'h.'+str(Wscan.value)+'degscan.'+casename+'.png')
plt.show()

'''
###

Vpix = (cosmo.kpc_comoving_per_arcmin(zz) * beam)**2 * cc/cosmo.H(zz) * dz
Vs = Vpix * (Wscan/beam)**2 * (Wnu/dnu)

sigma_pix = 2.5e6 * u.MJy*u.s**0.5/u.sr

Pnoise = Vpix*sigma_pix**2/tobs

def NN(kk,delta_k):
    return np.pi*kk*delta_k*Vs/(2.*np.pi)**2

print(Pnoise.to('MJy2 Mpc3/ sr2')/1.0e9)
print((Wscan/beam).to(''))

##
'''
