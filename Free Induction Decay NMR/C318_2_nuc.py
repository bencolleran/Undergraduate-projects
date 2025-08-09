#318 C318_2_nuc
import numpy as np
import scipy
import matplotlib.pyplot as plt
n=5001
m=5001
t=np.linspace(0,10,n)
dt=t[1]#step interval
Mt=np.zeros((7,len(t)))
Mt[6,:]=1
MeqA=10
MeqB=10
ToneA=20
ToneB=20
TtwoA=10
TtwoB=20
Minitial=np.array([0,MeqA,0,0,MeqB,0,1])
Bz=2#T
gammaA=2#HzT-1
gammaB=8#HzT-1
C=np.zeros((7,7))
C[0][1]=gammaA*Bz
C[1][0]=gammaA*-Bz
C[0][0]=-1/TtwoA
C[1][1]=-1/TtwoA
C[2][2]=-1/ToneA
C[2][6]=MeqA/ToneA
C[3][4]=gammaB*Bz
C[4][3]=gammaB*-Bz
C[3][3]=-1/TtwoB
C[4][4]=-1/TtwoB
C[5][5]=-1/ToneB
C[5][6]=MeqB/ToneB
#print(C)
C=C*2*np.pi#in radians per second
MxyTot=np.zeros(n,dtype = 'complex_')
for i in range(len(t)):
    Mt[:,i]=np.dot(scipy.linalg.expm(t[i]*C),Minitial)
    MxyTot[i]=Mt[0][i]+Mt[3][i]+1j*(Mt[1][i]+Mt[4][i])

v=np.linspace(0,20,m)
matrix=np.zeros((n,m),dtype='complex_')
V,T=np.meshgrid(v,t)
matrix=dt*np.exp(np.pi*1j*(2*V*T-0.5))
spectrum=np.dot(MxyTot,matrix)
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
'''
#projection on xy plane
plt.plot(MxyTot.real,MxyTot.imag)
plt.xlabel('X-projection')
plt.ylabel('Y-projection')
plt.title(f"spectrum with gyromagnetic ratio of {gammaA} and {gammaB} Hz/T")
plt.show()
'''

plt.plot(v,spectrum.real,label='real')
#plt.plot(v,spectrum.imag,label='imaginary')
plt.xlabel('frequency /Hz')
plt.ylabel('absorbance')
plt.title(f"spectrum with T2B {TtwoB} s")
#plt.legend()
plt.show()