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
resfac = 2. # resolution factor
res = int(resfac)


#%%

yeartosec=365.*24.*3600.
Lx=700.*10.**3.#/4. #length domain in meter
Ly=1500.*10.**3.#/4.
xstep=(1.5/(50.*res))*10.**6. #x step in meter, de 100. was eerst 50.
ystep=xstep
Tend=20000.*yeartosec #final time in seconds, years if yeartosec=1.
tstep=2.*yeartosec/10. #t step in seconds
Nx=int(Lx/(1*xstep)) #number of x steps
Ny=int(Ly/(1*ystep)) #number of y steps
Nt=int(Tend/tstep) #number of t steps
b=np.ones([Nx,Ny]) #bottom topography empty list


b[-10*res:-2*res,18*res:37*res] = 250.
b[-9*res:-3*res,19*res:36*res] = 500.
b[-7*res:-2*res,20*res:37*res] = 500.
b[-8*res:-3*res,19*res:36*res] = 1000.
b[-7*res:-3*res,19*res:36*res]=1500.

b[1*res:11*res,1*res:15*res] = 250.
b[2*res:11*res,2*res:14*res] = 500.
b[3*res:6*res,3*res:11*res]=1000.
b[6*res:10*res,3*res:13*res] = 1000.

b[2*res:14*res,13*res:49*res] = 50.
b[3*res:13*res,14*res:48*res] = 100.
b[4*res:12*res,15*res:47*res] = 150.
b[5*res:11*res,16*res:46*res] = 200.
b[6*res:10*res,17*res:45*res] = 250.
b[7*res:9*res,18*res:44*res] = 300.

b[12*res:24*res,0*res:18*res] = -500. #sea water

a=np.zeros([Nx,Ny]) #mass balance empty list
#Hinit=2000.
H=np.zeros([Nx,Ny,Nt])#*Hinit#bottom topography empty list
#H[:,:]=2000.
d=np.zeros([(Nx-1),(Ny-1)]) #bottom topography empty list
timeplot=[i for i in np.arange(0.,(Tend/yeartosec),(tstep/yeartosec))]#time list in years
#, can be used for plotting, no calculations (or with conversion)!
xplot=[i for i in np.arange(0.,Lx*10.**(-3.),xstep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
yplot=[i for i in np.arange(0.,Ly*10.**(-3.),ystep*10.**(-3.))]#x list in km
#, can be used for plotting, no calculations (or with conversion)!
#s=np.empty([Nx,Ny,Nt])

A=1.*10.**(-16.)/yeartosec
rho=910.
g=9.81
n=3.
#a=0.3/yeartosec
dadh=0.003/yeartosec
ELA=0.
amax=1./yeartosec #amax=3.17097919838e-08


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
    H[1:-1,1:-1,(l+1)]=H[1:-1,1:-1,l]+a[1:-1,1:-1]*tstep+(2.*(rho*g)**3.*A/5.)*ddds*tstep #let op!!! d is een element korter!!
#    H[1:-1,1:-1,(l+1)]=H[1:-1,1:-1,l]+a*tstep+(2.*(rho*g)**3.*A/5.)*ddds*tstep #let op!!! d is een element korter!!
    
    H[:,:,l+1] = np.where(H[:,:,l+1]<0,0,H[:,:,l+1])
    H[:,:,l+1] = np.where(b[:,:]<0.,0.,H[:,:,l+1])
    
#    H[-12:,:14]=0.
#    Hfl=H[:,:,l+1].flatten
#    print(np.sum(Hfl))
    
    h[:,:]=b+H[:,:,l+1]

#%%
X, Y = np.meshgrid(xplot[:-1], yplot[:]) 
Z=np.transpose(b[:,:])         
plt1=plt.figure(figsize = (6,8))
plt.contourf(X,Y,Z)#, [levels], **kwargs)
plt.colorbar()
#plt.savefig("base.png")
      
X, Y = np.meshgrid(xplot[:-1], yplot[:]) 
Z=np.transpose(H[:,:,-1])         
plt1=plt.figure(figsize = (6,8))
plt.contourf(X,Y,Z)#, [levels], **kwargs)
plt.colorbar()
plt.savefig("icesheetcontourres"+str(res)+".png")

X, Y = np.meshgrid(xplot[:-1], yplot[:])#-1 
B=np.transpose(b[:,:])#int(0.7*Nt)] 
plt2 = plt.figure(figsize = (7,8))
ax = plt2.gca(projection = '3d')
surf = ax.plot_surface(X, Y, B, cmap=cm.viridis, linewidth=0, antialiased=False)
plt.savefig("icebedrockres" + str(res) + ".png")

X, Y = np.meshgrid(xplot[:-1], yplot[:])#-1 
Z=np.transpose(H[:,:,-1])#int(0.7*Nt)] 
plt2 = plt.figure(figsize = (7,8))
ax = plt2.gca(projection = '3d')
surf = ax.plot_surface(X, Y, Z, cmap=cm.viridis, linewidth=0, antialiased=False)
for i in range(0,360,15):
    print(ax.azim)
    ax.view_init(azim=i)
    plt.savefig("icesheet3DHres"+str(res)+"rot" + str(i)+".png")

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
plt.savefig("icesheettimeseriesres"+str(res)+".png")


   