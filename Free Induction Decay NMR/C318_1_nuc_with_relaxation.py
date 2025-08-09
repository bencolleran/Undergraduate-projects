#318 C318_1_nuc_with_relaxation
import numpy as np
import scipy
import matplotlib.pyplot as plt
n=5001
m=5001
t=np.linspace(0,10,n)
dt=t[1]#step interval
Mt=np.zeros((4,len(t)))
Mt[3,:]=1
Meq=10
Tone=20
Ttwo=10
Minitial=np.array([0,Meq,0,1])
Bz=2#T
gamma=0.5#HzT-1
C=np.zeros((4,4))
C[0][1]=gamma*Bz
C[1][0]=gamma*-Bz
C[0][0]=-1/Ttwo
C[1][1]=-1/Ttwo
C[2][2]=-1/Tone
C[2][3]=Meq/Tone
#print(C)
C=C*2*np.pi#in radians per second
Mxy=np.zeros(n,dtype = 'complex_')
for i in range(len(t)):
    Mt[:,i]=np.dot(scipy.linalg.expm(t[i]*C),Minitial)
    Mxy[i]=Mt[0][i]+1j*Mt[1][i]
#Optional apodisation, causes signal to decay faster
param=5
for i in range(len(t)):
    Mxy[i]=Mxy[i]*np.exp(-param*t[i])

#Fourier transform of Mxy
v=np.linspace(0,10,m)
matrix=np.zeros((n,m),dtype='complex_')
V,T=np.meshgrid(v,t)
matrix=dt*np.exp(np.pi*1j*(2*V*T-0.5))
spectrum=np.dot(Mxy,matrix)
'''
print(np.trapz(spectrum.real,v,v[1]))
'''
'''
#path of the magnetisation
x=Mt[0]
y=Mt[1]
z=Mt[2]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, color='red')
ax.set_xlabel('X magnetisation')
ax.set_ylabel('Y magnetisation')
ax.set_zlabel('Z magnetisation')
plt.show()
'''
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

#Fourier transformed spectrum
plt.plot(v,spectrum.real,label='real')
plt.plot(v,spectrum.imag,label='imaginary')
plt.xlabel('frequency /Hz')
plt.ylabel('absorbance')
plt.title(f"spectrum with apodisation coefficient of {param}")
plt.legend()
plt.show()