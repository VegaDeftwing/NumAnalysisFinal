# Numerical Analysis Final
Final project for CSCE-440 Numerical Analysis at University of Nebraska - Lincoln

to run the code, you must first install some python libraries you will need
`sympy`, `scipy.io`, `numpy`, `rich`, `prett_errors`, `random`, `matplotlib`, and `soundfile`

if you are having issues, please ensure you are running the latest release of `scipy.io`

`wavinterpolate.py` hard codes some important values. In order to change the input wav, please redefine the name or path to the input .wav file on line 211, currently this line is set to `input_wav = 'NATEST24.wav'`, and this .wav file is inculded in this GitHub repository (or .zip if that's how you have gotten this code). Three input files are provide, `NATEST24.wav`,`NATEST16.wav`, and `NATEST8.wav` which all have their respective bit-depths in the name.

To change the number of samples that will be used for interpolation, change `downsample_level = 2` on line 19 to another value, such as 4 or 8. This value reflects how many samples will be dropped and replaced with interpolated samples, so at 8, 7 out of every 8 samples will de dropped, with only the 8th samples left to 'learn' the interpolation function.

Then just run the code with `python wavinterpolate.py` with a working directory that has the necessary files. Output files named `outputLinear.wav`, `outputQuad.wav`, and `outputCube.wav` will be generated, and a window should open showing the graphical views of the waveforms.


Code released under the MIT License:

---

Copyright 2021 Vega Carlson, Tyler Paul, and Sawyer Smith

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

---

Audio (NATEST*.wav) by Vega Carlson and released under CC BY-SA 4.0 