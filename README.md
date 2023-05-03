# GC×GC-MS Region of Interest Selection

## Ryland T. Giebelhaus, Michael D. Sorochan Armstrong, A. Paulina de la Mata, and James J. Harynuk*

### University of Alberta, Edmonton, Canada
#### james.harynuk@ualberta.ca

## 1.0 About
This is a Matlab algorithm developed to locate regions of interest (ROIs) within GC×GC-MS datasets. In brief, this algorithm takes a tensor, representing a GC×GC-MS chromatogram, and extracts the ROIs as smaller tensors for further chemometric analysis. This is especially handy for when a low number of components is desired. This algorithm is an extension of the 1D ROI algorithm described in Giebelhaus _et al._ [1] with some minor modifications. We found GC×GC-MS chromatograms to be significantly more sparse than 1D GC-MS chromatograms (added benefit of increased peak capacity), therefore the _pseudo_ F-ratio was having a challenging time dealing with sparse chromatograms (< 50 analytes). To address this, we adjusted the algorithm to first crudely select ROIs and isolate the noise (non-ROIs). We then calculated the standard deviation of the noise and autoscaled to this standard deviation for all moving windows. This is applied over multiple itterations and is user specified.

 
