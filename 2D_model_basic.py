# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 14:50:54 2018

@author: angelo en inge
"""

import matplotlib.pyplot as plt
import numpy as np
plt.close("all")

yeartosec=365.*24.*3600.
L=1.5*10.**6./4. #length domain in meter
xstep=(1.5/50.)*10.**6. #x step in meter, de 100. was eerst 50.
ystep=xstep
Tend=20000.*yeartosec #final time in seconds, years if yeartosec=1.
tstep=2.*yeartosec/10. #t step in seconds
Nx=int(L/xstep) #number of x steps
Ny=Nx
Nt=int(Tend/tstep) #number of t steps
b=np.zeros([Nx,Ny]) #bottom topography empty list
H=np.zeros([Nx,Ny,Nt]) #bottom topography empty list
d=np.zeros([(Nx-1),(Ny-1)]) #bottom topography empty list
#h=np.empty([Nx,Ny,Nt]) #bottom topography empty list
#for k in range (0,Nt):
#    b[:,:,k]=[i for i in np.arange(0.,L*10.**(-3.),xstep*10.**(-3.))]
timeplot=[i for i in np.arange(0.,(Tend/yeartosec),(tstep/yeartosec))]#time list in years
#, can be used for plotting, no calculations (or with conversion)!
xplot=[i for i in np.arange(0.,L*10.**(-3.),xstep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
yplot=[i for i in np.arange(0.,L*10.**(-3.),ystep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
s=np.empty([Nx,Ny,Nt])

A=1.*10.**(-16.)/yeartosec
rho=910.
g=9.81
n=3.
a=0.3/yeartosec
#Z8=(a*5.*L**4)/(2.*A*(rho*g)**3.)# is approx. 3,898.0 m.... 
#Z=Z8**(1./float(8.))
dNy=int(Ny-1)
dNx=int(Nx-1)


h=b+H[:,:,0]
for l in range (0,Nt-1):
    print(l)    

#    h=b+H
    d1=(1./4.)*(((h[1:,0:-1]-h[0:-1,0:-1])/xstep)+((h[1:,1:]-h[0:-1,1:])/xstep))**2.\
    +(1./4.)*(((h[0:-1,1:]-h[0:-1,0:-1])/ystep)+((h[1:,1:]-h[1:,0:-1])/ystep))**2.
    H5=(H[0:-1,0:-1,l]**5.+H[1:,0:-1,l]**5+H[0:-1,1:,l]**5.+H[1:,1:,l]**5)/4.
    d=d1*H5
    ddds=(1/xstep**2)*((1./2.)*(d[1:,1:]+d[1:,0:-1])*(h[2:,1:-1]-h[1:-1,1:-1])-(1./2.)*(d[0:-1,1:]+d[0:-1,0:-1])*(h[1:-1,1:-1]-h[0:-2,1:-1]))\
    +(1./ystep**2)*((1./2.)*(d[1:,1:]+d[0:1,1:])*(h[1:-1,2:]-h[1:-1,1:-1])-(1./2.)*(d[1:,0:-1]+d[0:-1,0:-1])*(h[1:-1,1:-1]-h[1:-1,0:-2]))
    H[1:-1,1:-1,(l+1)]=H[1:-1,1:-1,l]+a*tstep+(2.*(rho*g)**3.*A/5.)*ddds*tstep #let op!!! d is een element korter!!
    
#    H[0,:,:]=0.
#    H[int(Nx-1),:,:]=0.
#    H[:,0,:]=0.
#    H[:,(Ny-1),:]=0.
    
    
    
    H[:,:,l+1] = np.where(H[:,:,l+1]<0,0,H[:,:,l+1])
    #H[:,:,l+1] = np.where(H[:,:,l+1]>0,H[:,:,l],0)
    
#    for y in range (1,Nx-1):
#        for k in range (1,Ny-1):
#            if H[y,k,l]<0.:
#                H[y,k,l]=0.
    h[:,:]=b+H[:,:,l+1]
#        for k in range (1,Ny-1):
#            xgradplusx=(1./4.)*((s[y+1,k,l]-s[y,k,l])/xstep+(s[y+1,k+1,l]-s[y,k+1,l])/xstep)**2
#            ygradplusx=(1./4.)*((s[y,k+1,l]-s[y,k,l])/ystep+(s[y+1,k+1,l]-s[y+1,k,l])/ystep)**2
#            avsplusx=((s[y,k,l]+s[y+1,k,l]+s[y,k+1,l]+s[y+1,k+1,l])/4.)**5
#            dplusplus=(xgradplusx+ygradplusx)*avsplusx
#            xgradminx=(1./4.)*((s[y+1,k-1,l]-s[y,k-1,l])/xstep+(s[y+1,k,l]-s[y,k,l])/xstep)**2
#            ygradminx=(1./4.)*((s[y,k,l]-s[y,k-1,l])/ystep+(s[y+1,k,l]-s[y+1,k-1,l])/ystep)**2
#            avsminx=((s[y,k-1,l]+s[y+1,k-1,l]+s[y,k,l]+s[y+1,k,l])/4.)**5
#            dplusmin=(xgradminx+ygradminx)*avsminx
#            
#
#
#            
#            xgradplusy=(1./4.)*((s[y,k-1,l]-s[y-1,k-1,l])/xstep+(s[y,k,l]-s[y-1,k,l])/xstep)**2
#            ygradplusy=(1./4.)*((s[y-1,k,l]-s[y-1,k-1,l])/ystep+(s[y,k,l]-s[y,k-1,l])/ystep)**2
#            avsplusy=((s[y-1,k-1,l]+s[y,k-1,l]+s[y-1,k,l]+s[y,k,l])/4.)**5
#            dminmin=(xgradplusy+ygradplusy)*avsplusy
#            xgradminy=(1./4.)*((s[y,k,l]-s[y-1,k,l])/xstep+(s[y,k+1,l]-s[y-1,k+1,l])/xstep)**2
#            ygradminy=(1./4.)*((s[y-1,k+1,l]-s[y-1,k,l])/ystep+(s[y,k+1,l]-s[y,k,l])/ystep)**2
#            avsminy=((s[y-1,k,l]+s[y,k,l]+s[y-1,k+1,l]+s[y,k+1,l])/4.)**5
#            dminplus=(xgradminy+ygradminy)*avsminy
#            part1x=(tstep/(xstep**2.))*(1./2.)*(dplusplus+dplusmin)*(s[y+1,k,l]-s[y,k,l])
#            part2x=(tstep/(xstep**2.))*(1./2.)*(dminplus+dminmin)*(s[y,k,l]-s[y-1,k,l])
#            gradxdt=part1x-part2x
#            part1y=(tstep/(ystep**2.))*(1./2.)*(dplusplus+dminplus)*(s[y,k+1,l]-s[y,k,l])
#            part2y=(tstep/(ystep**2.))*(1./2.)*(dplusmin+dminmin)*(s[y,k,l]-s[y,k-1,l])
#
#            s[y,k,l+1]=s[y,k,l]+tstep+((part1x-part2x)+(part1y-part2y))

#%%      
X, Y = np.meshgrid(xplot[:-1], yplot[:-1]) 
Z=H[:,:,-1]#int(0.7*Nt)]           
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


   