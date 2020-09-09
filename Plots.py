import h5py
import os
import numpy as np
import math
import scipy
import pywt
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


EOS=['15H','125H','H','HB','B']
MASS=['135_135','125_146','125_125','121_151','118_155','117_156','116_158','112_140','107_146']
eq=['135_135','125_125']

#Define functions to calculate expected frequencies

#for q=1
def f20(M,R6):
    return 8.943+4.059*M-1.332*R6-.358*(M**2)-.182*R6*M+.048*(R6**2)

def fspir(M,R8):
    return 6.264+1.929*M-.645*R8+.881*(M**2)-.311*R8*M+.03*(R8**2)

def fpeak(M,R6):
    return 13.822-0.576*M-1.375*R6+.479*(M**2)-.073*R6*M+.044*(R6**2)

#for all cases
def f20_a(M,R6):
    return 9.586+4.09*M-1.427*R6+.048*(M**2)-.261*R6*M+.055*(R6**2)

def fspir_a(M,R8):
    return 5.846+1.75*M-.555*R8+1.002*(M**2)-.316*R8*M+.026*(R8**2)

def fpeak_a(M,R8):
    return 10.942-.369*M-.987*R8+1.095*(M**2)-.201*R8*M+.036*(R8**2)


#define functions for the analysis
def fre_do(x,y,mass):
    fd=fft(y)
    N=len(y)
    if (N % 2) == 1:
        N=N+1
    T=x[1]-x[0]
    xf = np.linspace(0.0, 1.0/(2.0*T), int(N/2))#/mass
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

    tuk=signal.tukey(len(ampl),0.03)
    dat=ampl*tuk

    fq,fd=fre_do(tim,dat,mass)

    mx=np.where(fd==np.amax(fd))[0][0]
    freq=fq[mx]
    amp=fd[mx]

    return fq,fd,tim,dat

#calculate R1.6 and R1.8 for all the EOS
m_r1=np.load('tid_def/15H.npy')
m_r2=np.load('tid_def/125H.npy')
m_r3=np.load('tid_def/H.npy')
m_r4=np.load('tid_def/HB.npy')
m_r5=np.load('tid_def/B.npy')

mx=np.amax(m_r1[0])
idx=np.where(m_r1[0]==mx)
idx=idx[0][0]
#cs1=spline(m_r1[0][1:idx],k_l1[0][1:idx])
cs11=spline(m_r1[0][1:idx],m_r1[1][1:idx])

mx=np.amax(m_r2[0])
idx=np.where(m_r2[0]==mx)
idx=idx[0][0]
#cs2=spline(m_r2[0][1:idx],k_l2[0][1:idx])
cs21=spline(m_r2[0][1:idx],m_r2[1][1:idx])

mx=np.amax(m_r3[0])
idx=np.where(m_r3[0]==mx)
idx=idx[0][0]
#cs3=spline(m_r3[0][1:idx],k_l3[0][1:idx])
cs31=spline(m_r3[0][1:idx],m_r3[1][1:idx])

mx=np.amax(m_r4[0])
idx=np.where(m_r4[0]==mx)
idx=idx[0][0]
#cs4=spline(m_r4[0][1:idx],k_l4[0][1:idx])
cs41=spline(m_r4[0][1:idx],m_r4[1][1:idx])

mx=np.amax(m_r5[0])
idx=np.where(m_r5[0]==mx)
idx=idx[0][0]
#cs5=spline(m_r5[0][1:idx],k_l5[0][1:idx])
cs51=spline(m_r5[0][1:idx],m_r5[1][1:idx])

r68=np.zeros((len(EOS),2))
i=0
for eos in EOS:
    if eos=='15H':
        r68[i,0]=cs11(1.6)*Length/1.0e5
        r68[i,1]=cs11(1.8)*Length/1.0e5

    elif eos=='125H':
        r68[i,0]=cs21(1.6)*Length/1.0e5
        r68[i,1]=cs21(1.8)*Length/1.0e5
    elif eos=='H':
        r68[i,0]=cs31(1.6)*Length/1.0e5
        r68[i,1]=cs31(1.8)*Length/1.0e5

    elif eos=='HB':
        r68[i,0]=cs41(1.6)*Length/1.0e5
        r68[i,1]=cs41(1.8)*Length/1.0e5

    elif eos=='B':
        r68[i,0]=cs51(1.6)*Length/1.0e5
        r68[i,1]=cs51(1.8)*Length/1.0e5

    i=i+1




#do the analysis and save the plots
#make the folders
if os.path.exists('results'):
    pass
else:
    os.mkdir('results')

if os.path.exists('results/q1'):
    pass
else:
    os.mkdir('results/q1')

if os.path.exists('results/q1/linear'):
    pass
else:
    os.mkdir('results/q1/linear')

if os.path.exists('results/q1/log'):
    pass
else:
    os.mkdir('results/q1/log')

if os.path.exists('results/all_q'):
    pass
else:
    os.mkdir('results/all_q')

if os.path.exists('results/all_q/linear'):
    pass
else:
    os.mkdir('results/all_q/linear')

if os.path.exists('results/all_q/log'):
    pass
else:
    os.mkdir('results/all_q/log')

if os.path.exists('results/3fig'):
    pass
else:
    os.mkdir('results/3fig')

if os.path.exists('results/3fig/linear'):
    pass
else:
    os.mkdir('results/3fig/linear')

if os.path.exists('results/3fig/log'):
    pass
else:
    os.mkdir('results/3fig/log')

if os.path.exists('results/spec'):
    pass
else:
    os.mkdir('results/spec')


#q=1


for eos in EOS:
        if eos=='15H':
            nmb=1
        elif eos=='125H':
            nmb=2
        elif eos=='H':
            nmb=3
        elif eos=='HB':
            nmb=4
        elif eos=='B':
            nmb=5

        for mas in eq:
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

            mas1=float(mas.split('_')[0])/100
            mas2=float(mas.split('_')[1])/100
            mastot=float(mas.split('_')[0])/100+float(mas.split('_')[1])/100
            q=mas1/mas2
            Mc=pow(q/pow(1+q,2),3/5)*mastot


            freq2,amp2,tim,post=analyze(strain,time,mastot)
            f_2=f20(Mc,r68[nmb-1,0])
            f_s=fspir(Mc,r68[nmb-1,1])
            f_p=fpeak(Mc,r68[nmb-1,0])
            f_0=2*f_p-f_2
            #print(f_2_a,f_s_a,f_p_a)



            fig=plt.figure()
            plt.plot((freq2*Frequency),amp2)
            ax=plt.subplot()
            ax.axvline(x=(f_p*Mc)*1000,color='r',label='peak')
            ax.axvspan((f_p*Mc)*1000-196, (f_p*Mc)*1000+196, alpha=0.3, color='grey')
            ax.axvline(x=(f_2*Mc)*1000,color='g',label='2-0')
            ax.axvspan((f_2*Mc)*1000-229, (f_2*Mc)*1000+229, alpha=0.3, color='yellow')
            ax.axvline((f_s*Mc)*1000,color='orange',label='spiral')
            ax.axvspan((f_s*Mc)*1000-286, (f_s*Mc)*1000+286, alpha=0.3, color='cyan')
            ax.axvline((f_0*Mc)*1000,linestyle="--",color='grey',label='2+0')
            plt.xlim(0,5000)
            plt.xlabel('Frequency (Hz)')
            plt.legend()
            plt.savefig('results/q1/linear/'+eos+'_'+mas+'.jpg')
            plt.close

            fig=plt.figure()
            plt.plot((freq2*Frequency),amp2)
            ax=plt.subplot()
            ax.axvline(x=(f_p*Mc)*1000,color='r',label='peak')
            ax.axvspan((f_p*Mc)*1000-196, (f_p*Mc)*1000+196, alpha=0.3, color='grey')
            ax.axvline(x=(f_2*Mc)*1000,color='g',label='2-0')
            ax.axvspan((f_2*Mc)*1000-229, (f_2*Mc)*1000+229, alpha=0.3, color='yellow')
            ax.axvline((f_s*Mc)*1000,color='orange',label='spiral')
            ax.axvspan((f_s*Mc)*1000-286, (f_s*Mc)*1000+286, alpha=0.3, color='cyan')
            ax.axvline((f_0*Mc)*1000,linestyle="--",color='grey',label='2+0')
            plt.xlim(0,5000)
            plt.xlabel('Frequency (Hz)')
            plt.yscale('log')
            plt.ylim(1e-3,1)
            plt.legend()
            plt.savefig('results/q1/log/'+eos+'_'+mas+'.jpg')
            plt.close()




#all cases


for eos in EOS:
        if eos=='15H':
            nmb=1
        elif eos=='125H':
            nmb=2
        elif eos=='H':
            nmb=3
        elif eos=='HB':
            nmb=4
        elif eos=='B':
            nmb=5

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

            mas1=float(mas.split('_')[0])/100
            mas2=float(mas.split('_')[1])/100
            mastot=float(mas.split('_')[0])/100+float(mas.split('_')[1])/100
            q=mas1/mas2
            Mc=pow(q/pow(1+q,2),3/5)*mastot


            freq2,amp2,tim,post=analyze(strain,time,mastot)
            f_2_a=f20_a(Mc,r68[nmb-1,0])
            f_s_a=fspir_a(Mc,r68[nmb-1,1])
            f_p_a=fpeak_a(Mc,r68[nmb-1,1])
            f_0_a=2*f_p_a-f_2_a
            #print(f_2_a,f_s_a,f_p_a)



            fig=plt.figure()
            plt.plot((freq2*Frequency),amp2)
            ax=plt.subplot()
            ax.axvline(x=(f_p_a*Mc)*1000,color='r',label='peak')
            ax.axvspan((f_p_a*Mc)*1000-196, (f_p_a*Mc)*1000+196, alpha=0.3, color='grey')
            ax.axvline(x=(f_2_a*Mc)*1000,color='g',label='2-0')
            ax.axvspan((f_2_a*Mc)*1000-229, (f_2_a*Mc)*1000+229, alpha=0.3, color='yellow')
            ax.axvline((f_s_a*Mc)*1000,color='orange',label='spiral')
            ax.axvspan((f_s_a*Mc)*1000-286, (f_s_a*Mc)*1000+286, alpha=0.3, color='cyan')
            ax.axvline((f_0_a*Mc)*1000,linestyle="--",color='grey',label='2+0')
            plt.xlim(0,5000)
            plt.xlabel('Frequency (Hz)')
            plt.legend()
            plt.savefig('results/all_q/linear/'+eos+'_'+mas+'.jpg')
            plt.close()

            fig=plt.figure()
            plt.subplot(212)
            plt.plot((freq2*Frequency),amp2)
            plt.xlim(0,5000)
            plt.xlabel('Frequency (Hz)')
            plt.legend(['Postmerger only'])
            plt.subplot(222)
            plt.plot(tim,post)
            plt.title('Postmerger')
            plt.subplot(221)
            plt.plot(time,strain)
            plt.title('Time Domain')
            plt.savefig('results/3fig/linear/'+eos+'_'+mas+'.jpg')
            plt.close()

            fig=plt.figure()
            plt.plot((freq2*Frequency),amp2)
            ax=plt.subplot()
            ax.axvline(x=(f_p_a*Mc)*1000,color='r',label='peak')
            ax.axvspan((f_p_a*Mc)*1000-196, (f_p_a*Mc)*1000+196, alpha=0.3, color='grey')
            ax.axvline(x=(f_2_a*Mc)*1000,color='g',label='2-0')
            ax.axvspan((f_2_a*Mc)*1000-229, (f_2_a*Mc)*1000+229, alpha=0.3, color='yellow')
            ax.axvline((f_s_a*Mc)*1000,color='orange',label='spiral')
            ax.axvspan((f_s_a*Mc)*1000-286, (f_s_a*Mc)*1000+286, alpha=0.3, color='cyan')
            ax.axvline((f_0_a*Mc)*1000,linestyle="--",color='grey',label='2+0')
            plt.xlim(0,5000)
            plt.xlabel('Frequency (Hz)')
            plt.yscale('log')
            plt.ylim(1e-3,1)
            plt.legend()
            plt.savefig('results/all_q/log/'+eos+'_'+mas+'.jpg')
            plt.close()

            fig=plt.figure()
            plt.subplot(212)
            plt.plot((freq2*Frequency),amp2)
            plt.yscale('log')
            plt.xlim(0,5000)
            plt.ylim(10**(-6),1)
            plt.xlabel('Frequency (Hz)')
            plt.legend(['Postmerger only'])
            plt.subplot(222)
            plt.plot(tim,post)
            plt.title('Postmerger')
            plt.subplot(221)
            plt.plot(time,strain)
            plt.title('Time Domain')
            plt.savefig('results/3fig/log/'+eos+'_'+mas+'.jpg')
            plt.close()

            fc =f_p_a
            dt=(tim[1]-tim[0])*Time*1000
            band = 2.5
            wavelet = 'cmor'+str(band)+'-'+str(fc)
            widths = fc/np.linspace(fc-1.0, fc+1.0, 400)/dt
            cwtmatr, freqs = pywt.cwt(post, widths, wavelet, dt)
            power = abs(cwtmatr)
            fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
            ax.pcolormesh(tim, freqs, power,cmap='jet')
            plt.savefig('results/spec/'+eos+'_'+mas+'.jpg')
            plt.close()
