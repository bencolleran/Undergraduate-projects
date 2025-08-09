#318 C318_2_nuc_exchange
import numpy as np
import scipy
import matplotlib.pyplot as plt
R=8.314462618153
temperature=425
k=(3.41*10**12)*np.exp(-85*1000/(R*temperature))#s-1
kf=k
kr=k
n=5001
m=5001
t=np.linspace(0,10,n)
dt=t[1]#step interval
Mt=np.zeros((7,len(t)))
Mt[6,:]=1
MeqA=10
MeqB=10
ToneA=100
ToneB=100
TtwoA=50
TtwoB=50
Minitial=np.array([0,MeqA,0,0,MeqB,0,1])
Bz=2#T
gammaA=2#HzT-1
gammaB=8#HzT-1
C=np.zeros((7,7))
C[0][1]=gammaA*Bz
C[1][0]=gammaA*-Bz
C[0][0]=-1/TtwoA-kf
C[1][1]=-1/TtwoA-kf
C[2][2]=-1/ToneA-kf
C[2][6]=MeqA/ToneA
C[3][4]=gammaB*Bz
C[4][3]=gammaB*-Bz
C[3][3]=-1/TtwoB-kr
C[4][4]=-1/TtwoB-kr
C[5][5]=-1/ToneB-kr
C[5][6]=MeqB/ToneB
C[0][3]=kr
C[1][4]=kr
C[2][5]=kr
C[3][0]=kf
C[4][1]=kf
C[5][2]=kf
#print(C)
C=C*2*np.pi#in radians per second
MxyA=np.zeros(n,dtype = 'complex_')
MxyB=np.zeros(n,dtype = 'complex_')
MxyTot=np.zeros(n,dtype = 'complex_')
for i in range(len(t)):
    Mt[:,i]=np.dot(scipy.linalg.expm(t[i]*C),Minitial)
    MxyA[i]=Mt[0][i]+1j*Mt[1][i]
    MxyB[i]=Mt[3][i]+1j*Mt[4][i]
    MxyTot[i]=Mt[0][i]+Mt[3][i]+1j*(Mt[1][i]+Mt[4][i])

v=np.linspace(0,20,m)
matrix=np.zeros((n,m),dtype='complex_')
V,T=np.meshgrid(v,t)
matrix=dt*np.exp(np.pi*1j*(2*V*T-0.5))
spectrumA=np.dot(MxyA,matrix)
spectrumB=np.dot(MxyB,matrix)
spectrumTot=np.dot(MxyTot,matrix)
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

#FID
plt.plot(t,Mt[0],label='x-component')
plt.plot(t,Mt[1],label='y-component')
#plt.plot(t,Mt[2],label='z-component')
#plt.title('Magnetisation components vs time')
plt.title(f"FID recorded at T={temperature} K")
plt.xlabel('time /s')
plt.ylabel('component')
plt.legend()
plt.show()

'''
#plt.plot(v,spectrumA.real,label='A')
#plt.plot(v,spectrumB.real,label='B')
plt.plot(v,spectrumTot.real,label='Total')
#plt.plot(v,spectrum.imag,label='imaginary')
plt.xlabel('frequency /Hz')
plt.ylabel('absorbance')
plt.title(f"spectrum recorded at T={temperature} K")
plt.legend()
plt.show()
'''
