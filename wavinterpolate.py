#!/bin/python
from scipy.io import wavfile
import scipy.io
import numpy as np
from rich import print
import pretty_errors

#####################################################################################
# DEFINE THE INTERPOLATION FUNCTIONS
#####################################################################################
#
# Each of these will take in the .wav segment, as well as where the 0'd out range
# starts and ends
# This -might- need to know if the file is 8, 16, or 24 bit as well, to compensate
# for the extra byte that numpy adds on 24 bit wavs.
# [TODO] are we allowed to use ones from pervious homework?

def linearInterpolate(wav_segment, zstart, zend):
    print("Running Linear Interpolation")
    if type(wav[0,0]) is np.int32:
        print("32bit internally, LSB's are off by a byte, Interpolation will only be over the 24 valid bits")
        # [TODO] interpolation function with byte compensation
    else:
        print(f"{type(wav[0,0])} internally, Interpolation can procede as normal")
        # [TODO] interpolation function, running normally

def NDDInterpolate(wav_segment, zstart, zend):
    print("Running Newtons Divided Differences Interpolation")
    if type(wav[0,0]) is np.int32:
        print("32bit internally, LSB's are off by a byte, Interpolation will only be over the 24 valid bits")
        # [TODO] interpolation function with byte compensation
    else:
        print(f"{type(wav[0,0])} internally, Interpolation can procede as normal")
        # [TODO] interpolation function, running normally

def LagrangeInterpolate(wav_segment, zstart, zend):
    print("Running Lagrange Polynomial Interpolation")
    if type(wav[0,0]) is np.int32:
        print("32bit internally, LSB's are off by a byte, Interpolation will only be over the 24 valid bits")
        # [TODO] interpolation function with byte compensation
    else:
        print(f"{type(wav[0,0])} internally, Interpolation can procede as normal")
        # [TODO] interpolation function, running normally

#####################################################################################
# 
#####################################################################################

# input24.wav is 24bit signed pcm, input16 is signed 16 bit, input8 is unsigned 8bit
input_wav = 'input24.wav'

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

# [TODO] We need to pick a sample range at a random starting point to 
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

if type(wav[0,0]) is np.int32:
    print("32bit internally, LSB's are off by a byte, Interpolation will only be over the 24 valid bits")
    # [TODO] do the 3 interpolation methods
else:
    print(f"{type(wav[0,0])} internally, Interpolation can procede as normal")
    # [TODO] do the 3 interpolation methods


# [TODO] write .wav file out

print("You can now go listen to the file to determine the quality of interpolation")

# [TODO] show the resulting waveforms after interpolation with Matplotlib against the input wavform