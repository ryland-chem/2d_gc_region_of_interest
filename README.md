# GC×GC-MS Region of Interest Selection

## Ryland T. Giebelhaus, A. Paulina de la Mata, and James J. Harynuk*

### University of Alberta, Edmonton, Canada
#### james.harynuk@ualberta.ca

## 1.0 About
This is a Matlab algorithm developed to locate regions of interest (ROIs) within GC×GC-MS datasets. In brief, this algorithm takes a tensor, representing a GC×GC-MS chromatogram, and extracts the ROIs as smaller tensors for further chemometric analysis. This is especially handy for when a low number of components is desired. This algorithm is an extension of the 1D ROI algorithm described in Giebelhaus _et al._ [1] with some minor modifications. We found GC×GC-MS chromatograms to be significantly more sparse than 1D GC-MS chromatograms (added benefit of increased peak capacity), therefore the _pseudo_ F-ratio was having a challenging time dealing with sparse chromatograms (< 50 analytes). To address this, we adjusted the algorithm to first crudely select ROIs and isolate the noise (non-ROIs). We then calculated the standard deviation of the noise and autoscaled to this standard deviation for all moving windows. This is applied over multiple itterations and is user specified.

## 2.0 Use
Download entire repository and unzip. Use the gcxgcROIMain.m function to perform ROI selection.

Users input data as an *I* by *J* by *K* tensor into the function (chromTensor), where *I* is acquisitions along first dimension, *J* is acquisitons along second dimension, and *K* is the _m/z_ of each ion.

Currently in v0.8 for testing.

### 2.1 Inputs

* **chromTensor**: *I* by *J* by *K* tensor represnting the chromatogram.
* **wndw**: Moving window size, ideal window size is approximatley the average peak width along the 2nd dimension. See [1] for more details.
* **cutOff**: Probability cutoff, on a range of 0-1. Typically 0.7 is a good place to start.
* **mode**: Depending on the data needs and sparcity of the tensor. Mode 1 just performs the standard ROI selection from [1] but in two dimensions. Mode 2 performs the itterations approach, which is helpful for sparce data.
* **plot**: Signals if the chromatogram should be output at the end, 1 to output and 0 to not print figure.
* **itters**: Number of itterations to perform if mode is 2. In our experience, 3-4 itterations is appropriate. Additional itterations do not appear to improve the data. Fill with any integer if mode is 1.

### 2.2 Outputs

dataOut is the only output, which is a structure containing a number of fields described below.

* **noiseDropped**: Original tensor with all noise regions removed and set to zero.
* **arrayPvals**: Array of p-values at every *I* and *J* in the chromatogram.
* **pValCutOff**: Array of p-values at every *I* and *J* in the chromatogram above the p-value cutoff.
* **ticDataReshaped**: *I* by *J* matrix representing the total ion count (TIC) of the input chromatogram.
* **originalTensor**: Original input tensor.
* **boolArray**: Boolean *I* by *J* matrix, where 1 means an acquisition is in an ROI and 0 means it is not.
* **labMatrix**: *I* by *J* matrix of all ROIs, with each ROI numbered.
* **numROIs**: Number of ROIs in data.
* **metaData**: Some metadata describing the user input parameters.
* **indivROIs**: A structure with each ROI extracted as an individual tensor. This allows users to access each tensor to perform chemometrics on it.

## 3.0 References

[1] Ryland T. Giebelhaus, Michael D. Sorochan Armstrong, A. Paulina de la Mata, and James J. Harynuk. Untargeted region of interest selection for gas chromatography - mass spectrometry data using a pseudo F-ratio moving window, _Journal of Chromatography A_, **2022**, 1682, 463499; 10.1016/j.chroma.2022.463499