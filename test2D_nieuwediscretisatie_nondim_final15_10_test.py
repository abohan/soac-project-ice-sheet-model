# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 14:50:54 2018

@author: angelo en inge
"""

import matplotlib.pyplot as plt
import numpy as np
plt.close("all")

yeartosec=365.*24.*3600.
L=1.5*10.**6./10. #length domain in meter
xstep=(1.5/50.)*10.**6. #x step in meter, de 100. was eerst 50.
ystep=xstep
Tend=20000.*yeartosec #final time in seconds, years if yeartosec=1.
tstep=2.*yeartosec/4. #t step in seconds
Nx=int(L/xstep) #number of x steps
Ny=Nx
Nt=200#int(Tend/tstep) #number of t steps
b=np.zeros([Nx,Ny,Nt]) #bottom topography empty list
H=np.zeros([Nx,Ny,Nt])
H[1:-1,1:-1,0]=1000. #bottom topography empty list
#h=np.empty([Nx,Ny,Nt]) #bottom topography empty list
#for k in range (0,Nt):
#    b[:,:,k]=[i for i in np.arange(0.,L*10.**(-3.),xstep*10.**(-3.))]
timeplot=[i for i in np.arange(0.,(Tend/yeartosec),(tstep/yeartosec))]#time list in years
#, can be used for plotting, no calculations (or with conversion)!
xplot=[i for i in np.arange(0.,L*10.**(-3.),xstep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
yplot=[i for i in np.arange(0.,L*10.**(-3.),ystep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
#s=np.zeros([Nx,Ny,Nt])

A=10.**(-23)#1.*10.**(-16.)/yeartosec
rho=910.
g=9.81
n=3.
a=0.3/yeartosec
Z8=(a*5.*L**4)/(2.*A*(rho*g)**3.)# is approx. 3,898.0 m.... 
Z=Z8**(1./float(8.))


#s = H + b
for l in range (0,Nt-1):
    print(l)
#    for f in range (1,Nx-1):
#        for g in range (1,Ny-1):
#            s[f,g,l] = b[f,g,l] + H[f,g,l]
    for y in range (1,Nx-1):
        for k in range (1,Ny-1):
#            s[y,k,l] = b[y,k,l] + H[y,k,l]
#            s[y+1,k,l] = b[y+1,k,l] + H[y+1,k,l]
#            s[y-1,k,l] = b[y-1,k,l] + H[y-1,k,l]
#            s[y,k-1,l] = b[y,k-1,l] + H[y,k-1,l]
#            s[y+1,k-1,l] = b[y+1,k-1,l] + H[y+1,k-1,l]
#            s[y-1,k-1,l] = b[y-1,k-1,l] + H[y-1,k-1,l]
#            s[y,k+1,l] = b[y,k+1,l] + H[y,k+1,l]
#            s[y+1,k+1,l] = b[y+1,k+1,l] + H[y+1,k+1,l]
#            s[y-1,k+1,l] = b[y-1,k+1,l] + H[y-1,k+1,l]
            xgradplusx=(1./4.)*((H[y+1,k,l]-H[y,k,l])/xstep+(H[y+1,k+1,l]-H[y,k+1,l])/xstep)**2
            ygradplusx=(1./4.)*((H[y,k+1,l]-H[y,k,l])/ystep+(H[y+1,k+1,l]-H[y+1,k,l])/ystep)**2
            avsplusx=((H[y,k,l]+H[y+1,k,l]+H[y,k+1,l]+H[y+1,k+1,l])/4.)**5
            dplusplus=(xgradplusx+ygradplusx)*avsplusx
            xgradminx=(1./4.)*((H[y+1,k-1,l]-H[y,k-1,l])/xstep+(H[y+1,k,l]-H[y,k,l])/xstep)**2
            ygradminx=(1./4.)*((H[y,k,l]-H[y,k-1,l])/ystep+(H[y+1,k,l]-H[y+1,k-1,l])/ystep)**2
            avsminx=((H[y,k-1,l]+H[y+1,k-1,l]+H[y,k,l]+H[y+1,k,l])/4.)**5
            dplusmin=(xgradminx+ygradminx)*avsminx
            

            xgradplusy=(1./4.)*((H[y,k-1,l]-H[y-1,k-1,l])/xstep+(H[y,k,l]-H[y-1,k,l])/xstep)**2
            ygradplusy=(1./4.)*((H[y-1,k,l]-H[y-1,k-1,l])/ystep+(H[y,k,l]-H[y,k-1,l])/ystep)**2
            avsplusy=((H[y-1,k-1,l]+H[y,k-1,l]+H[y-1,k,l]+H[y,k,l])/4.)**5
            dminmin=(xgradplusy+ygradplusy)*avsplusy
            xgradminy=(1./4.)*((H[y,k,l]-H[y-1,k,l])/xstep+(H[y,k+1,l]-H[y-1,k+1,l])/xstep)**2
            ygradminy=(1./4.)*((H[y-1,k+1,l]-H[y-1,k,l])/ystep+(H[y,k+1,l]-H[y,k,l])/ystep)**2
            avsminy=((H[y-1,k,l]+H[y,k,l]+H[y-1,k+1,l]+H[y,k+1,l])/4.)**5
            dminplus=(xgradminy+ygradminy)*avsminy
            part1x=(tstep/(xstep**2.))*(1./2.)*(dplusplus+dplusmin)*(H[y+1,k,l]-H[y,k,l])
            part2x=(tstep/(xstep**2.))*(1./2.)*(dminplus+dminmin)*(H[y,k,l]-H[y-1,k,l])
            #gradxdt=part1x-part2x
            part1y=(tstep/(ystep**2.))*(1./2.)*(dplusplus+dminplus)*(H[y,k+1,l]-H[y,k,l])
            part2y=(tstep/(ystep**2.))*(1./2.)*(dplusmin+dminmin)*(H[y,k,l]-H[y,k-1,l])

            H[y,k,l+1]=H[y,k,l]+a*tstep+((part1x-part2x)+(part1y-part2y))*2.*((rho*g)**3.)*A/5.
            if H[y,k,l+1]<0.:
                H[y,k,l+1]=0.
            #s[y,k,l+1] = b[y,k,l+1] + H[y,k,l+1]

#s[:,:,:]=b[:,:,:]+H[:,:,:]
                

#%%      
X, Y = np.meshgrid(xplot, yplot) 
Z=H[:,:,-1]           
plt1=plt.figure()
plt.contour(X,Y,Z)#, [levels], **kwargs)
plt.colorbar()
#plt.plot(s[:,:,int(0.7*Nt)])

locx=0.7
locxt=int(locx*Nx)
locy=0.7
locyt=int(locy*Ny)
plt.figure()
plt.plot (timeplot,H[locxt,locyt,:]) 
plt.title('2D ice sheet model timeseries, location:'+str(locx)+'Nx'+','+str(locy)+'Ny')
plt.grid()
plt.ylabel('s')
plt.xlabel('Time (yr)')  

print ('Z is:', Z)
   