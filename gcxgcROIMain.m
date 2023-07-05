%%Region of interest selection for 2D GCxGC-MS data using pseudo-fisher
%%ratios with connected components segmentation
%
%(c) 2022 Ryland T. Giebelhaus, A. Paulina de la
%Mata, James J. Harynuk
%
%Utilisation of a moving window function to calculate the f ratios for a
%particular region of interest. Smaller windows are more sensitive to
%smaller features, but may split larger regions of interest. Larger windows
%are less sensitive to small regions and create larger regions of interest.
%
%Code does the following
%-takes tensor IxJxK where I is 2nd dimension scans, J is 1st dimension
%scans, and K is ion intensities for m/z's in scan, moving window size
%(typically 10-45) and probability threshold (0.7 is a good place to
%start)
%-Returns array of probability values (0-1) across entire chromatographic
%plane and probabilities above probability threshold
%-TIC of input data
%-matrix of where each ROI is located
%-number of ROIs found
%-tensor of uploaded chromatogram with the noise (non ROIs) dropped.

%%%inputs
%%chromTensor: IxJxK where I is 2nd dimension scans, J is 1st dimension
%scans, and K is ion intensities for m/z's in scan
%%wndw: Moving window size, ideally the average peak width, start with 10
%and work up from there
%%cutOff: probability threshold, start with 0.7 and change depending on
%desired confidence level
%to do the gcxgc roi algorithm with one itteration, autoscaling with the
%standard deviation of the data in the moving window, set mode = 1. Set
%itters to any value

%to autoscale to the stdev of the noise of the chromatogram, set mode = 2,
%select how many times to itterate over the function (>= 2).

%to plot data set plot = 1, not plot = 0

%%%outputs
%%arrayPvals: Returns array of probability values (0-1) across entire
%chromatographic plane
%%pValCutOff: probability values above the probability threshold across
%chromatographic plane
%%ticDataReshaped: TIC formatted in shape of chromatographic plane
%%labMatrix: IxJ matrix showing whether a scan is in an ROI (>0) and which
% ROI a given scan is located in.
%%numROIs: The number of distinct ROIs found using the user inputs
%%noiseDropped: IxJxK tensor of uploaded data except the nonROI scans are
%dropped (set to zero)

function [dataOut] = gcxgcROIMain(chromTensor, wndw, cutOff, mode, plot, itters)

    %itterate over the gcxgc roi function multiple times, starting with
    %doing the normal roi selection by autoscaling to the stdev of the
    %wndw, then move to itterating over the gcxgcroi function multiple
    %times, using the stdev of the noise from the previous itteration to
    %autoscale to.
    
    %calculate the stdev of the entire tensor first
    sizeTensor = size(chromTensor);
    specdata = reshape(chromTensor, [sizeTensor(1)*sizeTensor(2), sizeTensor(3)]);
    stDevIn1 = std(specdata);
    
    %two modes, 1 and 2
    if mode == 1

        [arrayPvals, pValCutOff, ticDataReshaped, noiseDropped, boolArray] = gcxgcfroii(specdata, wndw, cutOff, 1, [], sizeTensor);

    elseif mode == 2

        for i = 1:itters
        
            if i == 1
        
                [arrayPvals, pValCutOff, ticDataReshaped, ~, boolArray] = gcxgcfroii(specdata, wndw*2, cutOff, 2, stDevIn1, sizeTensor);
        
            else
                
                %unfold arrayPvals
                vectPvals = reshape(arrayPvals, 1, []); %#ok

                sizeData = sizeTensor;
            
                %how many scans per mod
                scansPerMod = sizeData(1);

                %ions per scan
                ionsPerScan = sizeData(3);

                %preallocate noiseRegions
                noiseRegions = zeros(scansPerMod*sizeData(2), ionsPerScan);
    
                %unfold arrayPvals
                vectPvals = reshape(arrayPvals, 1, []);
                
                %how many modulctions
                numbMods = sizeData(2)/scansPerMod;
                
                %should be fast anyways but loop over the vector of pvals and find the
                %places that are noise and 
                for kk = 1:numbMods*scansPerMod
                    
                    %if less than the cutoff (noise)
                    if vectPvals(kk) < cutOff
            
                        noiseRegions(kk,:) = specdata(kk,:);
            
                    end
            
                end
                
                %remove all rows with 0s
                noiseRegions = noiseRegions(any(noiseRegions,2),:);
                
                %calculate the stDevNoise for each mass channel
                stDevNoise = std(noiseRegions);

                %[stDevNoise] = calcStdNoise(arrayPvals, chromTensor, cutOff);
                [arrayPvals, pValCutOff, ticDataReshaped, noiseDropped, boolArray] = gcxgcfroii(specdata, wndw, cutOff, 2, stDevNoise, sizeTensor);
        
            end
        
        end

    end

    %package up all the data into a structure
    dataOut.noiseDropped = noiseDropped;
    dataOut.arrayPvals = arrayPvals;
    dataOut.pValCutOff = pValCutOff;
    dataOut.ticDataReshaped = ticDataReshaped;
    dataOut.originalTensor = chromTensor;
    dataOut.boolArray = boolArray;

    %if user wants to plot data then plot data
    %display data
    [labMatrix, numROIs] = dispROIsgcxgc(dataOut, plot);

    dataOut.labMatrix = labMatrix;
    dataOut.numROIs = numROIs;

    %some info about what was done
    metaData = table(wndw, cutOff, mode, itters);
    dataOut.metaData = metaData;

    %going to put all the ROIs into their own tensors
    if dataOut.numROIs > 0

        for i = 1:dataOut.numROIs
            
            %change the field name in the struct each time
            fieldName = sprintf('ROI_%d', i);
            
            %run this function to get a tensor
            [extractTensor] = extractTensorROI(dataOut, i);
            
            %run this function to make the ROI into a smaller tensor
            [subtensor] = extract_blob(extractTensor);
            
            %add to the structure subTensors
            subTensors.(fieldName) = subtensor;

            %add to dataOut
            dataOut.indivROIs = subTensors;
    
        end
        
    end

end