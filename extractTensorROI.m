function [extractTensor] = extractTensorROI(dataOut, roiNum)
    
    %get the tensor from the structure
    tensorData = dataOut.noiseDropped;

    %get the size of the tensor
    sz = size(tensorData);
    
    %make an array to store the tensor in
    layerTensor = zeros(sz(1), sz(2)); %#ok
    extractTensor = zeros(sz(1), sz(2), sz(3));

    %make the mask
    maskRoi = dataOut.labMatrix == roiNum;
    maskNotRoi = ~maskRoi;


    for i = 1:sz(3)

        layerTensor = tensorData(:,:,i);

        layerTensor(maskNotRoi) = 0;

        extractTensor(:,:,i) = layerTensor;

    end

end