# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 14:50:54 2018

@author: angelo en inge
"""

import matplotlib.pyplot as plt
import numpy as np
plt.close("all")

yeartosec=1.#365.*24.*3600.
L=1.5*10.**6. #length domain in meter
xstep=(1.5/50.)*10.**6. #x step in meter
ystep=xstep
Tend=20000.*yeartosec #final time in seconds, years if yeartosec=1.
tstep=2.*yeartosec #t step in seconds
Nx=int(L/xstep) #number of x steps
Ny=Nx
Nt=int(Tend/tstep) #number of t steps
b=np.empty([Nx,Ny,Nt]) #bottom topography empty list
#for k in range (0,Nt):
#    b[:,k]=[i for i in np.arange(0.,L*10.**(-3.),xstep*10.**(-3.))]
timeplot=[i for i in np.arange(0.,(Tend/yeartosec),(tstep/yeartosec))]#time list in years
#, can be used for plotting, no calculations (or with conversion)!
xplot=[i for i in np.arange(0.,L*10.**(-3.),xstep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
yplot=[i for i in np.arange(0.,L*10.**(-3.),ystep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
s=np.empty([Nx,Ny,Nt])

A=1.*10.**(-16.)
rho=910.
g=9.81
n=3.
a=0.3

for l in range (0,Nt-1):
    print(l)
    for y in range (1,Nx-1):
        for k in range (1,Ny-1):
            xgradplusx=((s[y+1,k,l]-s[y,k,l])/xstep)**2
            ygradplusx=((s[y,k+1,l]-s[y,k,l])/ystep)**2
            avsplusx=((s[y,k,l]+s[y+1,k,l])/2.)**5
            dplusx=(xgradplusx+ygradplusx)*avsplusx
            xgradminx=((s[y,k,l]-s[y-1,k,l])/xstep)**2
            ygradminx=((s[y,k,l]-s[y,k-1,l])/ystep)**2
            avsminx=((s[y,k,l]+s[y-1,k,l])/2.)**5
            dminx=(xgradminx+ygradminx)*avsminx
            x2grad=(tstep/(xstep**2.))*(dplusx*(s[y+1,k,l]-s[y,k,l])-dminx*(s[y,k,l]-s[y-1,k,l]))
            
            xgradplusy=((s[y+1,k,l]-s[y,k,l])/xstep)**2
            ygradplusy=((s[y,k+1,l]-s[y,k,l])/ystep)**2
            avsplusy=((s[y,k,l]+s[y,k+1,l])/2.)**5
            dplusy=(xgradplusy+ygradplusy)*avsplusy
            xgradminy=((s[y,k,l]-s[y-1,k,l])/xstep)**2
            ygradminy=((s[y,k,l]-s[y,k-1,l])/ystep)**2
            avsminy=((s[y,k,l]+s[y,k-1,l])/2.)**5
            dminy=(xgradminy+ygradminy)*avsminy
            y2grad=(tstep/(ystep**2.))*(dplusy*(s[y,k+1,l]-s[y,k,l])-dminy*(s[y,k,l]-s[y,k-1,l]))

            s[y,k,l+1]=s[y,k,l]+tstep+(tstep/(xstep**2.))*x2grad+tstep+(tstep/(ystep**2.))*y2grad
            

#%%      
X, Y = np.meshgrid(xplot, yplot) 
Z=s[:,:,int(0.7*Nt)]           
plt1=plt.figure()
plt.contour(X,Y,Z)#, [levels], **kwargs)
#plt.plot(s[:,:,int(0.7*Nt)])
      