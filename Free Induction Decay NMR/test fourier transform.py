#test fourier transform
import numpy as np
import matplotlib.pyplot as plt
n=1000
m=1000
omega=2*np.pi
t=np.linspace(0,10,n)
dt=t[1]
Mxy=(np.cos(-omega*t+0.5*np.pi)+1j*np.sin(-omega*t+0.5*np.pi))*np.exp(-t)
v=np.linspace(0,10,m)
matrix=np.zeros((n,m),dtype='complex_')
V,T=np.meshgrid(v,t)
matrix=dt*np.exp(np.pi*1j*(2*V*T-0.5))
spectrum=np.dot(Mxy,matrix)
plt.plot(v,spectrum.imag)
#plt.plot(t,Mxy)
plt.show()


