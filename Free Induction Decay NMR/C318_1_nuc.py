#318 C318_1_nuc
import numpy as np
import scipy
import matplotlib.pyplot as plt
n=5001
m=2000
t=np.linspace(0,10,n)
dt=t[1]#step interval
Mt=np.zeros((3,len(t)))
Minitial=np.array([0,10,0])
Bz=2#T
gamma=8#HzT-1
C=np.zeros((3,3))
C[0][1]=gamma*Bz
C[1][0]=gamma*-Bz
C=C*2*np.pi#in radians per second
Mxy=np.zeros(n,dtype = 'complex_')
for i in range(len(t)):
    Mt[:,i]=np.dot(scipy.linalg.expm(t[i]*C),Minitial)
    Mxy[i]=Mt[0][i]+1j*Mt[1][i]
#fourier transform of Mxy
v=np.linspace(0,20,m)
matrix=np.zeros((n,m),dtype='complex_')
V,T=np.meshgrid(v,t)
matrix=dt*np.exp(np.pi*1j*(2*V*T-0.5))
spectrum=np.dot(Mxy,matrix)
'''
#FID
plt.plot(t,Mt[0],label='x-component')
plt.plot(t,Mt[1],label='y-component')
plt.plot(t,Mt[2],label='z-component')
plt.title('Magnetisation components vs time')
plt.xlabel('time /s')
plt.ylabel('component')
plt.legend()
plt.show()
'''
'''
#path of the magnetisation
x=Mxy.real
y=Mxy.imag
z=t
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, color='red')
ax.set_xlabel('X magnetisation')
ax.set_ylabel('Y magnetisation')
ax.set_zlabel('time')
plt.show()
'''
#Fourier transformed spectrum
plt.plot(v,spectrum.real,label='real')
plt.plot(v,spectrum.imag,label='imaginary')
plt.xlabel('frequency /Hz')
plt.ylabel('absorbance')
plt.title(f"spectrum with gyromagnetic ratio of {gamma} Hz/T")
plt.legend()
plt.show()
