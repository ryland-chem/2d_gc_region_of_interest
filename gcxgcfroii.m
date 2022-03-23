%%Region of interest selection for 2D GCxGC-MS data using pseudo-fisher ratios
%
%(c) 2021 Michael Sorochan Armstrong & Ryland T. Giebelhaus
%
%Utilisation of a moving window function to calculate the f ratios for a
%particular region of interest. Smaller windows are more sensitive to
%smaller features, but may split larger regions of interest. Larger windows
%are less sensitive to small regions and create larger regions of interest.
%Currently this returns a vector of probabilities for one-dimensional data,
%and it is up to the user to determine the significance level.
%
%Currently not optimised for speed.
%
%v1.01

%working on 1D for a bit more then coming back to 2D
%we need a way to seperate features within a ROI

function [arrayPvals, pValCutOff, ticData, ticDataReshaped] = gcxgcfroii(chromTensor, wndw, cutOff)

%need to transform the tensor into a linear array
sizeTensor = size(chromTensor);

specdata = reshape(chromTensor, [sizeTensor(1)*sizeTensor(2), sizeTensor(3)]);

%how many scans per mod
scansPerMod = sizeTensor(1);

%index 1 (where to start measuring from)
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

while indxCounter(2) <= numbScans
    
    data = specdata(indxCounter(1):indxCounter(2),:);

    %Initialisation
    sz = size(data);
    
    %Start at 1:wndw
    indx(1) = 1;
    indx(2) = wndw;
    
    mat = [];
    
    iter = 1;
    
    %Using a while loop instead of doing arithmetic; change this later
    while indx(2) <= sz(1)
        
        %C is the autoscaled matrix
        C = (data(indx(1):indx(2),:) - mean(data(indx(1):indx(2),:)))./std(data(indx(1):indx(2),:));
        
        %C
        C(isnan(C)) = 0;
        
        %s are the singular values
        s = svds(C,2);
        
        %f is the pseudo fisher ratio
        f(iter) = s(1)^2/s(2)^2; %#ok
        
        %matrix of probabilities
        mat(indx(1):indx(2),iter) = fcdf(f(iter),wndw,wndw); %#ok
        
        %increase the iteration number, change the region where the window is
        %active
        iter = iter + 1;
        indx(1) = indx(1) + 1;
        indx(2) = indx(2) + 1;
        
    end
    
    for ii = 1:sz(1)
        
        %ii is the observation number we are looping through
        %the degrees of freedom is the number of observations, since there is no
        %class information
        dof(ii) = sum(mat(ii,:) > 0); %#ok
        
        %retrieve the nonzero elements. mat is an off-diagonal matrix.
        nzx2 = nonzeros(mat(ii,:));
        
        for qq = 1:max(size(nzx2))
            %every probability is translated to a chi2 value with a dof of 1.
            x2mat(ii,qq) = chi2inv(nzx2(qq),1); %#ok
        end
        
        %sum the chi2 values
        x2vec(ii) = sum(x2mat(ii,:),2); %#ok
        
        %translate the chi2 values back to
        pv(ii) = chi2cdf(x2vec(ii),dof(ii)); %#ok
        
    end
    
    indxCounter = indxCounter + scansPerMod;
    
    %Change the orientation to a column
    pvals = pv';
    
    arrayPvals = [arrayPvals pvals];
    
    %converts values < cutoff to 0
    pValCutOff = arrayPvals;
    pValCutOff(pValCutOff < cutOff) = 0;
  
end
    
%need to drop data points that are below the pValCutOff
%itterate over the entire array
%snake through array, cause it is faster

%first calculate the TIC
%empty array for TIC data
ticData = [];

for i = 1:numbScans
    
    ticData(i) = sum(specdata(i,:));

end

%transpose to make it go the right direction
ticData = ticData';

%need to reshape so it is array that is printable with image()
ticDataReshaped = reshape(ticData, scansPerMod, []);

%need to use loop and logic to drop regions that are < pValVCutoff
%in the array of P vals

%total number of datapoints
totalData = numbScans * ionsPerScan;

%have to flip all arrays
% arrayPvals = flip(arrayPvals);
% pValCutOff = flip(pValCutOff);
% ticDataReshaped = flip(ticDataReshaped);

% % %%%%%%watershed segmentation here%%%%%
% % %calculate the gradient, find regions (edges) with extreme changes
% gmag = imgradient(pValCutOff);
% 
% %need to mark objects in the foreground
% se = strel('disk', 4000);
% Io = imopen(pValCutOff, se);
% 
% %opening by reconstruction
% Ie = imerode(pValCutOff, se);
% Iobr = imreconstruct(Ie, pValCutOff);
% 
% %followed by closing, helps to fill in regions that look like dougnuts
% Ioc = imclose(Io, se);
% 
% Iobrd = imdilate(Iobr, se);
% Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
% Iobrcbr = imcomplement(Iobrcbr);
% 
% %regional maxima calc
% %this affects the watershed here
% fgm = imregionalmax(Iobrcbr);
% 
% se2 = strel(ones(5,5));
% fgm2 = imclose(fgm,se2);
% fgm3 = imerode(fgm2,se2);
% 
% %again affects the watershed
% fgm4 = bwareaopen(fgm3, 5);
% 
% %compute background markers
% bw = imbinarize(Iobrcbr);
% 
% %watershed
% D = bwdist(bw);
% %DL has our different regions of interest in it
% DL = watershed(D);
% 
% %overlay probability cutoff with the watershed
% Lrgb = label2rgb(DL, 'jet','w','shuffle');
% image(pValCutOff, 'CDataMapping', 'scaled');
% hold on
% himage = imshow(Lrgb);
% himage.AlphaData = 0.3;

end
