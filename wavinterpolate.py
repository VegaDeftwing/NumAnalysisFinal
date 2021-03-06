#!/bin/python
import sympy
from scipy.io import wavfile
import numpy as np
from rich import print
import pretty_errors
import random
from matplotlib import pyplot as plt 
import math
import soundfile as sf


#####################################################################################
# Important Variables
#####################################################################################

# number of samples that will be decimated and reconsturcted
samples_to_injest = 50000
downsample_level = 2

assert (samples_to_injest%2==0),"Samples to injest must be an even number!"
assert (downsample_level%2==0),"Downsample level must be an even number!"

#####################################################################################
# DEFINE THE INTERPOLATION FUNCTIONS
#####################################################################################
#
# Each of these will take in the .wav segment, as well as where the 0'd out range
# starts and ends
# This -might- need to know if the file is 8, 16, or 24 bit as well, to compensate
# for the extra byte that numpy adds on 24 bit wavs.
# Keep in mind samples_to_injest is the number of samples both before and after, so
# it needs divided by two to look forward and ahead.
# the range that's interpolated will be
# wav[(zstart - samples_to_injest/2:zstart),:] and wav[zend:zend + samples_to_injest/2,:]
# as we don't want to 'learn' on the range we've just 0'd out.
# BUT we do need to keep in mind the x/time value jump, so that the interpolation
# doesn't think these two ranges are contiuous

def LinearInterpolate(samples_to_injest, zstart, zend):
    print("Running Linear Spline Interpolation")

    div2interp = np.zeros((int(samples_to_injest/downsample_level)))

    i = zstart
    j = 0
    while i < zstart+samples_to_injest:
        div2interp[j] = wav[i] 
        i += downsample_level
        j += 1

    xp = np.arange(zstart,zend,downsample_level)
    yp = div2interp
    xn = np.arange(zstart,zend-downsample_level,1)
    linearWav = np.copy(wav)

    x = sympy.symbols('x')

    y = []

    for i in range(1,len(xp)):
        y.append (((xp[i] - x) / (xp[i] - xp[i-1]))*yp[i-1] + ((x - xp[i-1])/(xp[i] - xp[i-1]))*yp[i])

    for i in range(zstart,zend-downsample_level):
        linearWav[i] = y[((i-zstart)//downsample_level)].subs(x,(i))

    return linearWav

def QuadInterpolate(samples_to_injest, zstart, zend):
    print("Running Quadratic Spline Interpolation")

    div2interp = np.zeros((int(samples_to_injest/downsample_level)))

    i = zstart
    j = 0
    while i < zstart+samples_to_injest:
        div2interp[j] = wav[i] 
        i += downsample_level
        j += 1

    np.set_printoptions(formatter={'int':str})

    xp = np.arange(zstart,zend,downsample_level)
    yp = div2interp
    xn = np.arange(zstart,zend-downsample_level,1)
    quadWav = np.copy(wav)

    x = sympy.symbols('x')

    y = []
    z = []
    z.append(0)
    
    for i in range(0,len(xp)-1):
        z.append ((-1)*z[i] + 2*((yp[i+1]-yp[i])/(xp[i+1]-xp[i])))
        
    for i in range(0, len(xp)-1):
        y.append ((((z[i+1]-z[i])/(2*(xp[i+1]-xp[i])))*(x-xp[i])**2)+z[i]*(x-xp[i])+yp[i])

    for i in range(zstart,zend-downsample_level):
        quadWav[i] = y[((i-zstart)//downsample_level)].subs(x,(i))

    return quadWav

def RCubeInterpolate(samples_to_injest, zstart, zend):
    print("Running Cubic Spline Interpolation")

    div2interp = np.zeros((int(samples_to_injest/downsample_level)))

    i = zstart
    j = 0
    while i < zstart+samples_to_injest:
        div2interp[j] = wav[i] 
        i += downsample_level
        j += 1

    np.set_printoptions(formatter={'int':str})

    xp = np.arange(zstart,zend,downsample_level)
    yp = div2interp
    rCubeWav = np.copy(wav)
    
    x = sympy.symbols('x')
    
    y = []
    b = []
    c = np.zeros(len(xp))
    d = []
    e = []
    alp = []
    r = 2+math.sqrt(3)
    h = downsample_level
    
    e.append(3*r/(2*(h**2))*(yp[1]-yp[0]))
    for i in range(1, len(xp)-1):
        e.append((3/(h**2))*(yp[i-1]-2*yp[i]+yp[i+1]))
    e.append(0)
    
    alp.append(e[0]/r)
    for i in range(1, len(xp)-1):
        alp.append((e[i]-alp[i-1])/r)
    alp.append(0)
    
    for i in reversed(range(0,len(xp)-1)):
        c[i] = alp[i]-(c[i+1]/r)
        
    for i in range(0,len(xp)-1):
        b.append ((yp[i+1]-yp[i])/h-((2*c[i]+c[i+1])*h)/3)
        
    for i in range(0,len(xp)-1):
        d.append((1/(3*h))*(c[i+1]-c[i]))

    for i in range(0,len(xp)-1):
        y.append (yp[i]+b[i]*x+c[i]*(x**2)+d[i]*(x**3))

    for i in range(zstart,zend-downsample_level):
        rCubeWav[i] = y[((i-zstart)//downsample_level)].subs(x,(i % downsample_level))
        
    return rCubeWav

#####################################################################################
# MAIN
#####################################################################################

def PlotWavs(length, start, end, mainWav, linearWav, quadWav, rCubeWav):
    #TODO save the image
    extra_space = 100 #how many samples to show before and after the 0'd out samples
    fig, axs = plt.subplots(4,2)
    fig.suptitle("Waveform Interpolation")
    x = np.arange(0,length+extra_space*2,1)
    #base waveform
    axs[0,0].set_title("Input Waveform")
    axs[0,0].plot(x, mainWav[start-extra_space:end+extra_space])
    axs[0,0].axvspan(extra_space, length+extra_space, color='red', alpha=.1)
    #interpolated waveforms
    axs[0,1].set_title("Linear Spline Interpolation")
    axs[0,1].plot(x, linearWav[start-extra_space:end+extra_space], 'tab:orange')
    axs[1,0].set_title("Quadratic Spline Interpolation")
    axs[1,0].plot(x, quadWav[start-extra_space:end+extra_space], 'tab:green')
    axs[1,1].set_title("R-Cubic Spline Interpolation")
    axs[1,1].plot(x, rCubeWav[start-extra_space:end+extra_space], 'tab:red')

    # Multi Graph Comparison
    axs[2,0].set_title("Compare all splines")
    axs[2,0].plot(x, mainWav[start-extra_space:end+extra_space]-quadWav[start-extra_space:end+extra_space], 'tab:green')
    axs[2,0].plot(x, mainWav[start-extra_space:end+extra_space]-linearWav[start-extra_space:end+extra_space], 'tab:orange')
    axs[2,0].plot(x, mainWav[start-extra_space:end+extra_space]-rCubeWav[start-extra_space:end+extra_space], 'tab:red')

    #resulting difference
    axs[2,1].set_title("Linear Spline Interpolation Difference")
    axs[2,1].plot(x, mainWav[start-extra_space:end+extra_space]-linearWav[start-extra_space:end+extra_space], 'tab:orange')
    axs[3,0].set_title("Quadratic Spline Interpolation Difference")
    axs[3,0].plot(x, mainWav[start-extra_space:end+extra_space]-quadWav[start-extra_space:end+extra_space], 'tab:green')
    axs[3,1].set_title("R-Cubic Spline Interpolation Difference")
    axs[3,1].plot(x, mainWav[start-extra_space:end+extra_space]-rCubeWav[start-extra_space:end+extra_space], 'tab:red')

    for ax in axs.flat:
        ax.set(xlabel='Time', ylabel='Amplitude')
    for ax in axs.flat:
        ax.label_outer()

    plt.show()



def SaveWavs(linearWav,quadWav,rCubeWav):
    print("You can now go listen to the file to determine the quality of interpolation")


# input24.wav is 24bit signed pcm, input16 is signed 16 bit, input8 is unsigned 8bit
input_wav = 'NATEST24.wav'

# Get the wave file data
samplerate, wav = wavfile.read(input_wav)
print(samplerate)
tempwav = np.zeros((wav.shape[0]))
tempwav = wav[:,0]
wav = tempwav
print(wav)
np.set_printoptions(formatter={'int':hex})
print(f"sample rate  = {samplerate}")
print(f"raw data     = {wav[:]}")

print(f"""\nlooking at a single sample and the way we're reading in the
data, threre might be extra 0's depending on the sample bit depth.
.wav files are commonly 8, 16, or 24 bit ints or 32bit float.
We'll avoid floats, so of the int types both 24's are stored
as int32's in numpy. 
For this file, samples are {type(wav[0])} internally\n""")

zstart = random.randrange(samples_to_injest/2,(wav.shape[0]-samples_to_injest-samples_to_injest),1)
#zstart = 3801749
zend = zstart + samples_to_injest
print(f"samples {zstart} to {zend}] will be downsampled and interpolated")

for i in range(wav.shape[0]-1, wav.shape[0]-1, -1):
    wav[i,0] = wav[i,0] | or_val
    if i==wav.shape[0]-(1+num_bits):
        wav[i,0] = wav[i,0] & ~or_val

sf.write("InputWaveSegment.wav", wav[zstart:zend], samplerate, 'PCM_24')

linearWav = LinearInterpolate(samples_to_injest, zstart, zend)
sf.write("outputLinear.wav", linearWav[zstart:zend], samplerate, 'PCM_24')
print("Linear Output .wav Written!")
quadWav = QuadInterpolate(samples_to_injest, zstart, zend)
sf.write("outputQuad.wav", quadWav[zstart:zend], samplerate, 'PCM_24')
print("Quadratic Output .wav Written!")
rCubeWav = RCubeInterpolate(samples_to_injest, zstart, zend)
sf.write("outputCube.wav", rCubeWav[zstart:zend], samplerate, 'PCM_24')
print("Cubic Output .wav Written!")

print("finished!")

PlotWavs(samples_to_injest, zstart, zend, wav, linearWav, quadWav, rCubeWav)
