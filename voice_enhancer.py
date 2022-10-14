# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 16:44:46 2022

@author: Robert
"""


import numpy as np
import scipy.io.wavfile as wavfile
import matplotlib.pyplot as plt
from math import ceil

def harm_amp(inputarr, freq1, freq2, fsampling): 
    
    f1 = freq1
    f2 = freq2
    fs = fsampling
    n1 = ceil(N*f1/fs)
    n2 = ceil(N*f2/fs)
    yenhanced = inputarr
    yenhanced[n1:n2+1] = 5*yenhanced[n1:n2+1]
    yenhanced[len(yenhanced)-n2:len(yenhanced)-n1+1] = 5*yenhanced[len(yenhanced)-n2:len(yenhanced)-n1+1]
    return yenhanced
    

letter = 'original'
fs,audio_dat = wavfile.read('%s.wav' %letter) #creates a tuple which contains the sampling frequency (fs) and the audio data

N = len(audio_dat) #number of samples
normalised = audio_dat/2**15
#creates an array of for all the audio data samples after the samples have been normalised in the for loop
#normalises the audio data between +1 and -1 by dividing by the greatest possible value
#2**15 as the bit-format is signed 16-bit integer 
    
#calculates fft. divided by the number of samples as  fft misses this step
audio_fftdat = np.fft.fft(normalised)/N

faxis = np.linspace(0,fs,N)
taxis = np.linspace(0,N/fs, N)

#TRY TO IDENTIFY PEAKS TO AMPLIFY IN THE LOG FREQUENCY SCALE
consonant_range = harm_amp(audio_fftdat,2000,8000,fs)    



better_audio = np.real(np.fft.ifft(consonant_range)*N*2**15)   #this produced a replica of original file, the file was named newbie
wavfile.write("improved234.wav",fs, better_audio.astype(np.int16)) #.astype(np.int16)




###############################PLOTTING#######################################

fig1 = plt.figure(1)
#plotting time domain w/ log scales
plt.plot(taxis, normalised)
plt.xlabel('Time (s)')
plt.ylabel('Normalised Amplitude')
plt.title('Signal in Time Domain: %s.wav' %letter)
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.show() 
fig1.savefig('TIME%s.pdf' %letter, format = 'pdf')

fig2 = plt.figure(2)
#plotting frequency domain w/ linear scales
plt.plot(faxis, abs(audio_fftdat))
plt.xlim(0,fs/2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('FFT spectrum Linear Scale: %s.wav' %letter)
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.show() 
#fig2.savefig('FFTLIN%s.pdf' %letter, format = 'pdf')

dB = 20*np.log10((audio_fftdat))
dB_Max = 20*np.log10(max(audio_fftdat))

fig3 = plt.figure(3)
#plotting frequency domain w/ log scales
plt.plot(faxis, dB - dB_Max) #2**15
plt.xlim(80,fs/2)
plt.xscale('log')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dB)')
plt.title('FFT spectrum Log Scale: %s.wav' %letter)
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.show() 
fig3.savefig('FFTLOG%s.pdf' %letter, format = 'pdf')


