import h5py
import os
import numpy as np
import math
import scipy
from scipy import signal
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline as spline
from scipy.fftpack import fft, fftshift ,ifft,rfft,fftfreq,rfftfreq
c=2.9979e10
G=6.67408e-8
Msun=1.989e33
Length = G*Msun/c**2
Time = Length/c
Frequency=1/Time


#define functions for the analysis
def fre_do(x,y,mass):
    fd=fft(y)
    N=len(y)
    if (N % 2) == 1:
        N=N+1
    T=x[1]-x[0]
    xf = np.linspace(0.0, 1.0/(2.0*T), int(N/2))/mass
    fq=fftfreq(len(y))
    mask=fq>=0
    fd=2.0*(fd/N)
    fd=fd[mask]
    fd=abs(fd)
    return xf,fd

def analyze(rhM,time,mass):



    peaks,prop=scipy.signal.find_peaks(abs(rhM))
    ampls=rhM[peaks]
    merg=np.amax(abs(ampls))
    merg=np.where(abs(ampls)==merg)
    merg=int(merg[0])
    t0=peaks[merg]

    ampl=rhM[t0:]
    tim=time[t0:]

    #ampl=rhM
    #tim=time

    tuk=signal.tukey(len(ampl))
    dat=ampl*tuk

    fq,fd=fre_do(tim,dat,mass)

    mx=np.where(fd==np.amax(fd))[0][0]
    freq=fq[mx]
    amp=fd[mx]

    return fq,fd,tim,dat

EOS=['15H','125H','H','HB','B','SFHo']
MASS=['135_135','125_146','125_125','121_151','118_155','117_156','116_158','112_140','107_146']


if os.path.exists('results'):
    pass
else:
    os.mkdir('results')


#do the analysis and save the plots
for eos in EOS:
    if eos!='SFHo':
        for mas in MASS:
            f=open('data/'+eos+'_'+mas,'r')
            lines=f.readlines()
            result1=[]
            result2=[]
            for x in lines:
                for i in range(len(x.split(' '))):
                    if x.split(' ')[i]!='':
                        result1.append(x.split(' ')[i])
                        for j in range(i+1,len(x.split(' '))):
                            if x.split(' ')[j]!='':
                                result2.append(x.split(' ')[j])
                                break
                        break

            time=np.zeros(len(result1))
            strain=np.zeros(len(result1))
            for i in range(len(result1)):
                time[i]=float(result1[i])
                strain[i]=float(result2[i])

            mastot=float(mas.split('_')[0])/100+float(mas.split('_')[1])/100
            freq2,amp2,tim,post=analyze(strain,time,mastot)
            fig=plt.figure()
            plt.subplot(212)
            plt.plot((freq2*Frequency),amp2)
            plt.xlim(0,3000)
            plt.xlabel('Frequency (Hz)')
            plt.legend(['Postmerger only'])
            plt.subplot(222)
            plt.plot(tim,post)
            plt.title('Postmerger')
            plt.subplot(221)
            plt.plot(time,strain)
            plt.title('Time Domain')
            plt.savefig('results/'+eos+'_'+mas+'.jpg')
            
