# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 13:34:59 2018

@author: angelo en inge
"""

import matplotlib.pyplot as plt
import numpy as np
plt.close("all")

yeartosec=1.#365.*24.*3600.
L=1.5*10.**6. #length domain in meter
xstep=(1.5/50.)*10.**6. #x step in meter
Tend=20000.*yeartosec #final time in seconds, years if yeartosec=1.
tstep=2.*yeartosec #t step in seconds
Nx=int(L/xstep) #number of x steps
Nt=int(Tend/tstep) #number of t steps
b=np.empty([Nx,Nt]) #bottom topography empty list
for k in range (0,Nt):
    b[:,k]=[i for i in np.arange(0.,L*10.**(-3.),xstep*10.**(-3.))]
timeplot=[i for i in np.arange(0.,(Tend/yeartosec),(tstep/yeartosec))]#time list in years
#, can be used for plotting, no calculations (or with conversion)!
xplot=[i for i in np.arange(0.,L*10.**(-3.),xstep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
h=np.empty([Nx,Nt])

A=1.*10.**(-16.)
rho=910.
g=9.81
n=3.
a=0.3
dconst=((rho*g)**n)*2.*A/(n+2.)

for l in range (0,Nt-1):
    print(l)
    for y in range (1,Nx-1):
        dgrad1=(np.abs((h[y+1,l]+b[y+1,l]-h[y,l]-b[y,l])/xstep))**(n-1.)
        dgrad2=((h[y,l]+h[y+1,l])/2.)**(n+2.)
        dplus=dconst*dgrad1*dgrad2
        dgrad1min=(np.abs((h[y,l]+b[y,l]-h[y-1,l]-b[y-1,l])/xstep))**(n-1.)
        dgrad2min=((h[y-1,l]+h[y,l])/2.)**(n+2.)
        dmin=dconst*dgrad1min*dgrad2min
        fac1=h[y+1,l]+b[y+1,l]-h[y,l]-b[y,l]
        fac2=h[y,l]+b[y,l]-h[y-1,l]-b[y-1,l]
        h[y,l+1]=h[y,l]+(tstep/(xstep**2.))*(dplus*fac1-dmin*fac2)+a*tstep

#%%        
plt1=plt.figure()
plt.plot(xplot,h[:,int(0.7*Nt)])
      