# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:40:48 2022

@author: Robert
"""
import numpy as np
import scipy.io.wavfile as wavfile
import matplotlib.pyplot as plt
from math import ceil

fs,audio_dat = wavfile.read('improved.wav')

normalised = np.zeros(len(audio_dat))  #creates an array of for all the audio data samples after the samples have been normalised in the for loop
N = len(audio_dat) #number of samples

for i in range(N): 
    normalised[i] = audio_dat[i]/2**15 #normalises the audio data between +1 and -1 by dividing by the greatest possible value
    #2**15 as the bit-format is signed 16-bit integer  

    
#calculates fft. divided by the number of samples as  fft misses this step
audio_fftdat = np.fft.fft(normalised)/N

def highpass_filter(inputarr, pass_freq, fsampling):
    f_pass = pass_freq
    f_end = fsampling/2
    fs = fsampling
    n1 = ceil(N*f_pass/fs)
    n2 = ceil(N*f_end/fs)
    filtered = inputarr
    filtered[n1:n2+1] = 0
    filtered[len(filtered)-n2:len(filtered)-n1+1] = 0
    return filtered
    
def distort():
 np.tanh

#def mixer():
    

filter_dat = highpass_filter(audio_fftdat, 2500, fs)

#distortion needs to be done in the time domain
dat_for_distortion = np.fft.ifft(filter_dat)*N*2**15

faxis = np.linspace(0,fs,N)
taxis = np.linspace(0,N/fs, N)

###############################PLOTTING#######################################

fig1 = plt.figure(1)
#plotting time domain w/ log scales
plt.plot(taxis, normalised)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Signal in Time Domain: improved.wav' )
plt.legend() # Shows the plot legend
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
#plt.savefig('TIME%s.pdf' %letter, format = 'pdf')
plt.show() 
#fig1.savefig('TIME%s.pdf' %letter, format = 'pdf')

fig2 = plt.figure(2)
#plotting frequency domain w/ linear scales
plt.plot(faxis, abs(filter_dat))
plt.xlim(0,fs)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('FFT spectrum Linear Scale: improved.wav')
plt.legend() # Shows the plot legend
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
#plt.savefig('FFTLIN%s.pdf' %letter, format = 'pdf')
plt.show() 
#fig2.savefig('FFTLIN%s.pdf' %letter, format = 'pdf')

dB = 20*np.log10((audio_fftdat))
dB_Max = 20*np.log10(max(audio_fftdat))

fig3 = plt.figure(3)
#plotting frequency domain w/ log scales
plt.plot(faxis, dB - dB_Max) #2**15
plt.xlim(1,fs)
plt.xscale('log')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dB)')
plt.title('FFT spectrum Log Scale: improved.wav')
plt.legend() # Shows the plot legend
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
#plt.savefig('FFTLOG%s.pdf' %letter, format = 'pdf')
plt.show() 
#fig3.savefig('FFTLOG%s.pdf' %letter, format = 'pdf')

fig4 = plt.figure(4)
#plotting time domain w/ log scales
plt.plot(taxis, dat_for_distortion)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Signal in Time Domain: improved.wav' )
plt.legend() # Shows the plot legend
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
#plt.savefig('TIME%s.pdf' %letter, format = 'pdf')
plt.show() 
#fig1.savefig('TIME%s.pdf' %letter, format = 'pdf')


