  %%Region of interest selection for 2D GCxGC-MS data using pseudo-fisher
%%ratios with connected components segmentation
%
%(c) 2022 Ryland T. Giebelhaus, Michael D.S. Armstrong, A. Paulina de la
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


function [arrayPvals, pValCutOff, ticDataReshaped, noiseDropped, arrayBools] = gcxgcfroii(specdata, wndw, cutOff, StdevIMP, stdIn, sizeTensor)

%bool to print graph
%at the start so user can input then let run
%prompt = 'Output graph (y/n)';
%choicePrint = input(prompt, 's');

% choicePrint = 'y';

%need to transform the tensor into a linear array
%sizeTensor = size(chromTensor);

%specdata = reshape(chromTensor, [sizeTensor(1)*sizeTensor(2), sizeTensor(3)]);

%take standard deviation of everything
stDevAll = std(specdata);

%how many scans per mod
scansPerMod = sizeTensor(1);

%index 1 (where to start measuring from)
indxCounter = [];
indxCounter(1) = 1;

%where to stop counting
indxCounter(2) = sizeTensor(1);

%calculate size of the specData
sizeData = size(specdata);

%number of scans
numbScans = sizeData(1);

%ions per scan
ionsPerScan = sizeData(2);

%how many modulctions
numbMods = numbScans/scansPerMod;

% %prepopulating array for speed
arrayPvals = [];
arrayBools = [];

    apvIndx = 1;

while indxCounter(2) <= numbScans

    
    data = specdata(indxCounter(1):indxCounter(2),:);

    %Initialisation
    sz = size(data);
    
    %Start at 1:wndw
    indx(1) = 1;
    indx(2) = wndw;
    
    mat = [];
    
    iter = 1;

    %inputs for while loop:
    %indx
    %sz
    %StdevIMP
    %data
    %stDevAll
    %stdIn
    %f
    
    while indx(2) <= sz(1)
        
        %C is the autoscaled matrix

        if StdevIMP == 1

            C = (data(indx(1):indx(2),:) - mean(data(indx(1):indx(2),:)))./std(data(indx(1):indx(2),:));

        elseif StdevIMP == 0

            C = (data(indx(1):indx(2),:) - mean(data(indx(1):indx(2),:)))./stDevAll;

        elseif StdevIMP == 2

            C = (data(indx(1):indx(2),:) - mean(data(indx(1):indx(2),:)))./stdIn;

        end
        
        %C
        C(isnan(C)) = 0;
        
        %s are the singular values
        s = svds(C,2);
        
        %f is the fisher ratio
        f(iter) = s(1)^2/s(2)^2; %#ok
        
        %matrix of probabilities
        mat(indx(1):indx(2),iter) = fcdf(f(iter),wndw,wndw); %#ok
        
        %increase the iteration number, change the region where the window is
        %active
        iter = iter + 1;
        indx(1) = indx(1) + 1;
        indx(2) = indx(2) + 1;
        
    end

    x2mat = zeros(sz(1), wndw*2);
    
    for ii = 1:sz(1)
        
        %ii is the observation number we are looping through
        %the degrees of freedom is the number of observations, since there is no
        %class information
        dof(ii) = sum(mat(ii,:) > 0); %#ok
        
        %retrieve the nonzero elements. mat is an off-diagonal matrix.
        nzx2 = nonzeros(mat(ii,:));
        nzx2 = nzx2';

        [x2mat] = chi2invROI_mex(nzx2, ii, x2mat);
        
%         for qq = 1:max(size(nzx2))
%             %every probability is translated to a chi2 value with a dof of 1.
%             x2mat(ii,qq) = chi2inv(nzx2(qq),1);
%         end
        
        %sum the chi2 values
        x2vec(ii) = sum(x2mat(ii,:),2); %#ok
        
        %translate the chi2 values back to
        pv(ii) = chi2cdf(x2vec(ii),dof(ii)); %#ok
        
    end
    
    indxCounter = indxCounter + scansPerMod;
    
    %Change the orientation to a column
    pvals = pv';

    %need to do some boolean work
    boolVect = pvals;
    boolVect(boolVect < cutOff) = 0;
    boolVect(boolVect > 0) = 1;

    [~, ~, boolVect] = spacesBetweenCleanUp(boolVect);

    arrayBools = [arrayBools boolVect']; %#ok
    
    arrayPvals = [arrayPvals pvals]; %#ok
    
    %converts values < cutoff to 0
    pValCutOff = arrayPvals;
    pValCutOff(pValCutOff < cutOff) = 0;

    apvIndx = 1 + apvIndx;
  
end %end of the while loop.
    
%need to drop data points that are below the pValCutOff
%itterate over the entire array
%snake through array, cause it is faster

%first calculate the TIC
%empty array for TIC data
ticData = zeros(1, size(specdata, 1));

for i = 1:numbScans
    
    ticData(i) = sum(specdata(i,:));

end

%transpose to make it go the right direction
ticData = ticData';

%need to reshape so it is array that is printable with image()
ticDataReshaped = reshape(ticData, scansPerMod, []);

%connected components segmentation here

%Convert the image to Black and White (pValCutOff)
%BWpVal = pValCutOff > 0;
BWpVal = arrayBools > 0;

%calculate connected components
conComp = bwconncomp(BWpVal);

%Creating label matrix for each ROI
labMatrix = labelmatrix(conComp);
% 
% %number of ROIs
% numROIs = max(labMatrix(:));

% %display a graph if the user requests one. Shows the ROIS overlaid on the
% %TIC of the input chromatogram.
% if choicePrint == 'y'
%     
%     Lrgb = label2rgb(labMatrix,'jet','w','shuffle');
%     
%     clims = [10 5E5];
%     colormap jet;
%     imagesc((ticDataReshaped),clims);
%     set(gca,'YDir','normal');
%     ylabel("2nd Dimension Acquisitions"); xlabel("1st Dimension Acquisitions");
%     hold on
%     himage = imshow(Lrgb);
%     himage.AlphaData = 0.8;
% 
% else
%     
% end

%itterate over the entire chromatogram and drop noise regions giving output
%of just ROIs
%populate tensor with zeros for speed
noiseDropped = zeros(scansPerMod*numbMods, ionsPerScan);
noiseRegions = zeros(scansPerMod*numbMods, ionsPerScan);

for ii = 1:numbMods*scansPerMod
    
    if labMatrix(ii) > 0
       
        noiseDropped(ii,:) = specdata(ii,:);
        
    else
        
        noiseDropped(ii,:) = zeros(1, ionsPerScan);
        noiseRegions(ii,:) = specdata(ii,:);
        
    end
    
end

%refold the noiseDropped into a tensor
noiseDropped = reshape(noiseDropped, [scansPerMod, numbMods, ionsPerScan]);

end
