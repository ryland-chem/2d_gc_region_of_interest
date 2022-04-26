%Convert the image to Black and White (pValCutOff)
BWpVal = pValCutOff > 0;

%calculate connected components
conComp = bwconncomp(BWpVal);

%Creating label matrix for each ROI
labMatrix = labelmatrix(conComp);

%number of ROIs
numROIs = max(labMatrix(:));

%Display the image of unique ROIs
imshow(label2rgb(labMatrix,'jet','k','shuffle'));