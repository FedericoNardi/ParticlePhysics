import numpy as np
import matplotlib.pyplot as plt
"""
#Differential Cross section for s=(14GeV)^2
E=7
ChargeSquare = 4*np.pi/137
K= ChargeSquare**2/(32*np.pi*4*E**2)
x= np.linspace(-0.9,0.9,100)
differentialCS = K*8*( (1+x)**2/16 + (1-x)**2/16 + (4 + (1+x)**2)/(1-x)**2 - (1+x)**2/(8*(1-x)) )
differentialCS_annihilation = K*8*( (1+x)**2/16 + (1-x)**2/16 )
differentialCS_scattering = K*8*( (4 + (1+x)**2)/(1-x)**2 )
differentialCS_crossTerm = K*8*( - (1+x)**2/(8*(1-x)) )

plt.figure()
plt.plot(x, differentialCS)
#plt.plot(x, differentialCS_annihilation)
#plt.plot(x, differentialCS_scattering)
#plt.plot(x, differentialCS_crossTerm)
plt.title("Bhabha scattering - differential cross section")
plt.xlabel(r"$\cos\theta$")
plt.ylabel(r"$d\sigma / d\cos\theta\quad[GeV^{-2}]$")
plt.grid(linestyle='--')
#plt.show()
plt.savefig("figures/differentialCS")

plt.figure()
plt.plot(x, differentialCS_annihilation)
plt.title(r"Differential cross section - $s$ channel")
plt.xlabel(r"$\cos\theta$")
plt.ylabel(r"$d\sigma / d\cos\theta\quad[GeV^{-2}]$")
plt.grid(linestyle='--')
#plt.show()
plt.savefig("figures/annihilationCS")

plt.figure()
plt.plot(x, differentialCS_scattering)
plt.title(r"Bhabha scattering - $t$ channel")
plt.xlabel(r"$\cos\theta$")
plt.ylabel(r"$d\sigma / d\cos\theta\quad[GeV^{-2}]$")
plt.grid(linestyle='--')
#plt.show()
plt.savefig("figures/scatteringCS")

plt.figure()
plt.plot(x, differentialCS_crossTerm)
plt.title(r"Bhabha scattering - Cross term")
plt.xlabel(r"$\cos\theta$")
plt.ylabel(r"$d\sigma / d\cos\theta\quad[GeV^{-2}]$")
plt.grid(linestyle='--')
#plt.show()
plt.savefig("figures/crossCS")


"""
#Total Cross section
"""
from scipy import integrate

ChargeSquare = 4*np.pi/137

CrossSection = np.zeros(100)
E = np.linspace(1,15,100)
cutoff=0.001
x= np.linspace(-1,1-cutoff,100)
for i in range(100):
	K= ChargeSquare**2/(32*np.pi*4*E[i]**2)
	differentialCS = K*8*( (1+x)**2/16 + (1-x)**2/16 + (4 + (1+x)**2)/(1-x)**2 - (1+x)**2/(8*(1-x)) )
	CrossSection[i] = integrate.simps(differentialCS,x)
plt.figure()
plt.plot(2*E, CrossSection)
plt.title("Bhabha scattering - total cross section")
plt.xlabel(r"$\sqrt{s}$ $[GeV]$")
plt.ylabel(r"$\sigma$ $[GeV^{-2}]$")
plt.grid(linestyle='--')
#plt.show()
plt.savefig("figures/totalCS")

# Muon pair production
x= np.linspace(-1,1,100)
CrossSection_muon = np.zeros(100)
for i in range(100):
	K= ChargeSquare**2/(32*np.pi*4*E[i]**2)
	differentialCS = K*8*( (1+x)**2/16 + (1-x)**2/16 )
	CrossSection_muon[i] = integrate.simps(differentialCS,x)
plt.figure()	
plt.plot(2*E, CrossSection_muon)
plt.title("Muon pair production - total cross section")
plt.xlabel(r"$\sqrt{s}$ $[GeV]$")
plt.ylabel(r"$\sigma$ $[GeV^{-2}]$")
plt.grid(linestyle='--')
#plt.show()
plt.savefig("figures/muontotalCS")



# Number of events at sqrt(s)=14GeV
CMEnergy = 7
K= ChargeSquare**2/(32*np.pi*4*CMEnergy**2)
#differentialCS = K*8*( (1+x)**2/16 + (1-x)**2/16 + (4 + (1+x)**2)/(1-x)**2 - (1+x)**2/(8*(1-x)) )
differentialCS = K*8*( (1+x)**2/16 + (1-x)**2/16 ) # For muon production
Cross_Section = integrate.simps(differentialCS,x)
# 1 GeV^-2 = 0.3894 mb
Luminosity = 10 #pb^-1
Efficiency = 0.5
NumEvents = Cross_Section*0.3894e9*Luminosity*Efficiency
print("-----CM Energy: 14 GeV-----")
print(NumEvents, "events")
print(Cross_Section)
"""

# Comparison with data
from scipy import integrate
ChargeSquare = 4*np.pi/137
CMEnergy = 29/2

x1=np.linspace(-0.55,0.55,50)
x2=np.linspace(0.75,0.85,50)
x3=np.linspace(-0.85,-0.75,50)


K= ChargeSquare**2/(32*np.pi*4*CMEnergy**2)
#differentialCS = K*8*( (1+x)**2/16 + (1-x)**2/16 + (4 + (1+x)**2)/(1-x)**2 - (1+x)**2/(8*(1-x)) )
differentialCS1 = K*8*( (1+x1)**2/16 + (1-x1)**2/16 )#+ (4 + (1+x1)**2)/(1-x1)**2 - (1+x1)**2/(8*(1-x1)) )
differentialCS2 = K*8*( (1+x2)**2/16 + (1-x2)**2/16 )#+ (4 + (1+x2)**2)/(1-x2)**2 - (1+x2)**2/(8*(1-x2)) )
differentialCS3 = K*8*( (1+x3)**2/16 + (1-x3)**2/16 )#+ (4 + (1+x3)**2)/(1-x3)**2 - (1+x3)**2/(8*(1-x3)) )
 # For muon production
Cross_Section = integrate.simps(differentialCS1,x1)
Cross_Section = integrate.simps(differentialCS2,x2)
Cross_Section += integrate.simps(differentialCS3,x3)
# 1 GeV^-2 = 0.3894 mb
Luminosity = 19.6 #pb^-1
Efficiency = 0.88
NumEvents = Cross_Section*0.3894e9*Luminosity*Efficiency
print("-----CM Energy: ",2*CMEnergy," GeV-----")
print(NumEvents, "events")
print(Cross_Section)



