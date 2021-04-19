#!/bin/python
from scipy.io import wavfile
import scipy.io
from scipy.interpolate import CubicSpline
import numpy as np
from rich import print
import pretty_errors
import random
from matplotlib import pyplot as plt 


#####################################################################################
# Important Variables
#####################################################################################

# number of samples that will be 0'd out and replaced with the interpolated 'guess'
length_to_interpolate = 100
# number of samples used to 'learn' the interpolation, half forwards, half backwards
samples_to_injest = 1000

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

    linearWav = np.copy(wav)
    linearWav[zstart:zend,:] = 0 
    if type(wav[0,0]) is np.int32:
        print("32bit internally, LSB's are off by a byte, Interpolation will only be over the 24 valid bits")
        #TODO interpolation function with byte compensation
    else:
        print(f"{type(wav[0,0])} internally, Interpolation can procede as normal")
        #TODO interpolation function, running normally

    return linearWav

def QuadInterpolate(samples_to_injest, zstart, zend):
    print("Running Quadratic Spline Interpolation")

    quadWav = np.copy(wav)
    quadWav[zstart:zend,:] = 0 
    if type(wav[0,0]) is np.int32:
        print("32bit internally, LSB's are off by a byte, Interpolation will only be over the 24 valid bits")
        #TODO interpolation function with byte compensation
    else:
        print(f"{type(wav[0,0])} internally, Interpolation can procede as normal")
        #TODO interpolation function, running normally

    return quadWav

def RCubeInterpolate(samples_to_injest, zstart, zend):
    print("Running R Cubic Spline Interpolation")

    rCubeWav = np.copy(wav)
    rCubeWav[zstart:zend,:] = 0
    if type(wav[0,0]) is np.int32:
        print("32bit internally, LSB's are off by a byte, Interpolation will only be over the 24 valid bits")
        #TODO interpolation function with byte compensation
    else:
        print(f"{type(wav[0,0])} internally, Interpolation can procede as normal")
        #TODO interpolation function, running normally

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
    axs[0,0].plot(x, mainWav[start-extra_space:end+extra_space,0])
    axs[0,0].axvspan(extra_space, length+extra_space, color='red', alpha=.1)
    #interpolated waveforms
    axs[0,1].set_title("Linear Spline Interpolation")
    axs[0,1].plot(x, linearWav[start-extra_space:end+extra_space,0], 'tab:orange')
    axs[1,0].set_title("Quadratic Spline Interpolation")
    axs[1,0].plot(x, quadWav[start-extra_space:end+extra_space,0], 'tab:green')
    axs[1,1].set_title("R-Cubic Spline Interpolation")
    axs[1,1].plot(x, rCubeWav[start-extra_space:end+extra_space,0], 'tab:red')

    # Establish a baseline using built in fuction
    # This shows that most interpolations will give at least a phase shift.
    xp1 = np.arange(0,samples_to_injest/2,1)
    xp2 = np.arange(length+samples_to_injest,length+samples_to_injest+samples_to_injest/2)
    xp  = np.concatenate((xp1,xp2))
    yp1 = mainWav[int(start - samples_to_injest/2):start,0]
    yp2 = mainWav[end:int(end + samples_to_injest/2),0]
    yp  = np.concatenate((yp1,yp2))
    xn = np.arange(length/2,length/2+length,1)
    interp = CubicSpline(xp,yp)
    lazyWav = np.copy(wav)
    lazyWav[start:end,0] = interp(xn)
    axs[2,0].set_title("Interpolation using SciPy Cubic Spline, baseline")
    axs[2,0].plot(x, lazyWav[start-extra_space:end+extra_space,0], 'tab:blue')
    axs[2,0].plot(x, mainWav[start-extra_space:end+extra_space,0]-lazyWav[start-extra_space:end+extra_space,0], 'tab:orange') 

    #resulting difference
    axs[2,1].set_title("Linear Spline Interpolation Difference")
    axs[2,1].plot(x, mainWav[start-extra_space:end+extra_space,0]-linearWav[start-extra_space:end+extra_space,0], 'tab:orange')
    axs[3,0].set_title("Quadratic Spline Interpolation Difference")
    axs[3,0].plot(x, mainWav[start-extra_space:end+extra_space,0]-quadWav[start-extra_space:end+extra_space,0], 'tab:green')
    axs[3,1].set_title("R-Cubic Spline Interpolation Difference")
    axs[3,1].plot(x, mainWav[start-extra_space:end+extra_space,0]-rCubeWav[start-extra_space:end+extra_space,0], 'tab:red')

    for ax in axs.flat:
        ax.set(xlabel='Time', ylabel='Amplitude')
    for ax in axs.flat:
        ax.label_outer()

    plt.show()



def SaveWavs(linearWav,quadWav,rCubeWav):
    #TODO make something to save the modified .wav files
    print("You can now go listen to the file to determine the quality of interpolation")


# input24.wav is 24bit signed pcm, input16 is signed 16 bit, input8 is unsigned 8bit
input_wav = 'NATEST24.wav'

# Get the wave file data
samplerate, wav = wavfile.read(input_wav)
np.set_printoptions(formatter={'int':hex})
print(f"sample rate      = {samplerate}")
print(f"raw data (left)  = {wav[:, 0]}")
print(f"raw data (right) = {wav[:, 1]}")

print(f"""\nlooking at a single sample and the way we're reading in the
data, threre might be extra 0's depending on the sample bit depth.
.wav files are commonly 8, 16, or 24 bit ints or 32bit float.
We'll avoid floats, so of the int types both 24's are stored
as int32's in numpy. 
For this file, samples are {type(wav[0,0])} internally\n""")

# We need to pick a sample range at a random starting point to 
# intentiaonally zero-out, the interpolation will act on the data before 
# and after this range to give us a function that we can try to repair.
#
#  1 -|         ,-'''-.
#     |      ,-'       `-.           
#     |    ,'             `.
#     |  ,'                 `.
#     | /                     \
#     |/                       \
# ----+-------------------------\--------------------------
#     |          __           __ \          __           /  __
#     |          ||/2         ||  \        3||/2        /  2||
#     |                            `.                 ,'
#     |                              `.             ,'
#     |                                `-.       ,-'
# -1 -|                                   `-,,,-'
#
# So, here we could se the values in this sine to something like this, where
# the values around the top of the sine wave get dropped to 0.
#
#  1 -|           
#     |      ,-'        -.           
#     |    ,'             `.
#     |  ,'                 `.
#     | /                     \
#     |/                       \
# ----+--------||||||||---------\--------------------------
#     |          __           __ \          __           /  __
#     |          ||/2         ||  \        3||/2        /  2||
#     |                            `.                 ,'
#     |                              `.             ,'
#     |                                `-.       ,-'
# -1 -|                                   `-,,,-'
#
# The goal would be to 'fix' the file to restore the sine wave.
# -------------------------------------------------------------
# get a random segment over which we want to interpolate
# leaving enough room on both sides to 'learn' from
zstart = random.randrange(samples_to_injest/2,(wav.shape[0]-length_to_interpolate-samples_to_injest),1)
zend = zstart + length_to_interpolate
print(f"samples {zstart} to {zend}] will be replaced with the interploated values")

linearWav = LinearInterpolate(samples_to_injest, zstart, zend)
quadWav = QuadInterpolate(samples_to_injest, zstart, zend)
rCubeWav = RCubeInterpolate(samples_to_injest, zstart, zend)

PlotWavs(length_to_interpolate, zstart, zend, wav, linearWav, quadWav, rCubeWav)

#TODO we need some sort of evaluation metric, maybe a mix of looking at the difference between the waves,
# the integral, and something to account for phase shift?