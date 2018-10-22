# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 14:50:54 2018

@author: angelo en inge
"""

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
plt.close("all")

#%%
r=6.37e6
ndeg=360.
pi=3.14159265359
correction_factor=2*np.pi*r/ndeg
print(23.*correction_factor)
print(40.*correction_factor*np.cos(pi*72./180.))


#%%

yeartosec=365.*24.*3600.
Lx=700.*10.**3.#/4. #length domain in meter
Ly=1500.*10.**3.#/4.
xstep=(1.5/50.)*10.**6. #x step in meter, de 100. was eerst 50.
ystep=xstep
Tend=20000.*yeartosec #final time in seconds, years if yeartosec=1.
tstep=2.*yeartosec#/10. #t step in seconds
Nx=int(Lx/xstep) #number of x steps
Ny=int(Ly/ystep) #number of y steps
Nt=int(Tend/tstep) #number of t steps
b=np.ones([Nx,Ny]) #bottom topography empty list
b[-7:-3,19:36]=1500.
b[2:11,3:13]=500.
b[4:13,15:46]=-200.
b[7:9,18:44]=-500.
b[-12:,:14]=-1000.#oceaan
a=np.zeros([Nx,Ny]) #bottom topography empty list
Hinit=2000.
H=np.ones([Nx,Ny,Nt])*Hinit#bottom topography empty list
#H[:,:]=2000.
d=np.zeros([(Nx-1),(Ny-1)]) #bottom topography empty list
timeplot=[i for i in np.arange(0.,(Tend/yeartosec),(tstep/yeartosec))]#time list in years
#, can be used for plotting, no calculations (or with conversion)!
xplot=[i for i in np.arange(0.,Lx*10.**(-3.),xstep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
yplot=[i for i in np.arange(0.,Ly*10.**(-3.),ystep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
s=np.empty([Nx,Ny,Nt])

A=1.*10.**(-16.)/yeartosec
rho=910.
g=9.81
n=3.
#a=0.3/yeartosec
dadh=0.003/yeartosec
ELA=1500.
amax=1./yeartosec #amax=3.17097919838e-08
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
    a[:,:]=(h[:,:]-ELA)*dadh
    a[:,:]=np.where(a[:,:]>amax,amax,a[:,:])
#    print(a[14,14])
    H[1:-1,1:-1,(l+1)]=H[1:-1,1:-1,l]+a[1,-1]*tstep+(2.*(rho*g)**3.*A/5.)*ddds*tstep #let op!!! d is een element korter!!
#    H[1:-1,1:-1,(l+1)]=H[1:-1,1:-1,l]+a*tstep+(2.*(rho*g)**3.*A/5.)*ddds*tstep #let op!!! d is een element korter!!
    
    H[:,:,l+1] = np.where(H[:,:,l+1]<0,0,H[:,:,l+1])
    H[-12:,:14]=0.
#    Hfl=H[:,:,l+1].flatten
#    print(np.sum(Hfl))
    
    h[:,:]=b+H[:,:,l+1]

#%%
X, Y = np.meshgrid(xplot[:-1], yplot[:]) 
Z=np.transpose(b[:,:])         
plt1=plt.figure()
plt.contourf(X,Y,Z)#, [levels], **kwargs)
plt.colorbar()
#plt.savefig("base.png")
      
X, Y = np.meshgrid(xplot[:-1], yplot[:]) 
Z=np.transpose(H[:,:,-1])         
plt1=plt.figure()
plt.contourf(X,Y,Z)#, [levels], **kwargs)
plt.colorbar()
plt.savefig("icesheetcontourHinit"+str(Hinit)+"ELA"+str(ELA)+"tstest.png")

X, Y = np.meshgrid(xplot[:-1], yplot[:])#-1 
Z=np.transpose(H[:,:,-1])#int(0.7*Nt)] 
plt2 = plt.figure()
ax = plt2.gca(projection = '3d')
surf = ax.plot_surface(X, Y, Z, cmap=cm.viridis, linewidth=0, antialiased=False)
#Axes3D.plot_surface(X,Y,Z)
plt.savefig("icesheet3DHinit"+str(Hinit)+"ELA"+str(ELA)+"tstest.png")

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
plt.savefig("icesheettimeseriesHinit"+str(Hinit)+"ELA"+str(ELA)+"tstest.png")


   