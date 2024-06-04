# -*- coding: utf-8 -*-A
import numpy as np
import matplotlib.pyplot as plt
from Acheson_equation import Acheson_coeffs
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc('font', family='sans-serif')
mpl.rc('font', serif='Helvetica')
plt.rc('font', size=12)          # controls default text sizes
plt.rcParams["axes.formatter.limits"] = [-3,4] #1x10^3 et non 100
plt.rcParams["axes.labelsize"] = 'large'#large
plt.rcParams["axes.titlesize"] = 'large'
plt.rcParams["xtick.labelsize"] = 'large'
plt.rcParams["ytick.labelsize"] = 'large'
plt.rcParams["legend.fontsize"] ='small'


theta=0.775*np.pi/2 #1.0 equator
g= 6.0e13 #6.5-7.5e13
Omega =5.468e3# in s^-1
N=4967#3900#Omega/2
q=1.12#1.6
rho = 3.7e14 # in g/cm^3
r=7e5 #in cm
dOmega=q*Omega/r
dOmegar2=(q+2)*Omega*r
Bp=3.45e17 #in G
V=Bp/np.sqrt(4*np.pi*rho)
print("Omega",Omega,"w_a",V/r)
cs2 = (1.26e10)**2 #From EOS
gamma = 2.19#2.19  #test for now, EOS dependence later on
#Hypothesis Bp \propto r^p
p=2.0
s=-2.0 # magnetic field dependence on z^s
l=2*np.pi/1e5
m=1
I=1j
min_wave = l
Npoints=2000
n = np.logspace(np.log10(min_wave),3,Npoints,dtype='float64') #in cm^-1
#k = np.arange(1,350,dtype='float64')
dE= gamma*N**2/g*(np.sin(theta)-l/n*np.cos(theta)) #Spruit 1999 #=0 for Polytropic EOS, EOS dependence later on
coeffs = np.arange(7,dtype='complex128')
plot_growth = np.zeros((Npoints,6))
plot_growth2 = np.zeros((Npoints,6))
plot_growth_diss= np.zeros_like(dissipation)
wavelength= np.zeros_like(dissipation)
z=r*np.cos(theta)
for k in range(len(dissipation)):
    eta=dissipation[k]
    nu = dissipation[k]
    kappa= dissipation[k]
    print(k,eta)
    for j in range(n.shape[0]):
        dF=(p-1)/r - l/n[j]*s/z #d np.log(B/rho/r)/dh #neglecting rho dependence
        dQ=(p+1)/r - l/n[j]*s/z#d np.log(B r )/dh
        G= g*np.sin(theta) - g*np.cos(theta)*(l/n[j])
        s2 = n[j]**2+l**2
        coeffs_acheson=np.array(Acheson_coeffs(Omega,q,N,rho,g,r,theta,Bp,eta,nu,kappa,cs2,gamma,p,s,l,n[j],s2,m),dtype=np.complex128)
        coeffs_acheson_l_large=np.array(Acheson_coeffs(Omega,q,N,rho,g,r,theta,Bp,eta,nu,kappa,cs2,gamma,p,s,n[j],l,s2,m),dtype=np.complex128)
        sigma2=np.roots(coeffs_acheson_l_large)
        sigma = np.roots(coeffs_acheson)
        plot_growth[j,0]=(np.imag(sigma).max())
        plot_growth[j,1]=np.sqrt(s2)
        plot_growth2[j,0]=(np.imag(sigma2).max())
        plot_growth2[j,1]=np.sqrt(s2)
    max1=(plot_growth[:,0].max())
    max2=(plot_growth2[:,0].max())
    case=np.argmax((max1,max2))
    print("case=",case)
    if case == 0:
        plot_growth_diss[k]=max1
        wavelength[k]=plot_growth[np.argmax(plot_growth[:,0]),1]
    else:
        plot_growth_diss[k]=max2
        wavelength[k]=plot_growth2[np.argmax(plot_growth2[:,0]),1]
    #print(plot_growth[j])
    #if(plot_growth[j] <= 0.0):
     #   print('stable')

plt.figure()
plt.loglog(dissipation,plot_growth_diss,label="Growth rate")
plt.hlines(V/r,dissipation[0],dissipation[-1],label=r"$\omega_A$",color='r')
plt.vlines(V/r/l**2,1e-4,1e7,color='k',label="dissipation limit",ls="--")
plt.legend()
plt.xlabel(r"Dissipation [cm$^2$ s$^-1$]")
plt.ylabel(r"Growth rate [s$^-1$]")
plt.tight_layout()
plt.xlim(dissipation[0],dissipation[-1])
plt.ylim(1e3,3e4)
plt.savefig("/Users/asalze/Documents/Articles/MyArticles/2024_TS_one_zone_model/Figure_8.pdf")
