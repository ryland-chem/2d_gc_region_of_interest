function [labMatrix, numROIs] = dispROIsgcxgc(dataOut, plot)

    ticDataReshaped = dataOut.ticDataReshaped;
    %boolArray = dataOut.boolArray;
    %pValCutOff = dataOut.pValCutOff;
    
    boolArray = dataOut.boolArray;

    figure;

    %Convert the image to Black and White (pValCutOff)
    BWpVal = boolArray > 0;
    
    %calculate connected components
    conComp = bwconncomp(BWpVal);
    
    %Creating label matrix for each ROI
    labMatrix = labelmatrix(conComp);
    
    %number of ROIs
    numROIs = max(labMatrix(:));

    if plot == 1
        
        Lrgb = label2rgb(labMatrix,'jet','w','shuffle');
        
        clims = [1000 5E5];
        colormap jet;
        imagesc((ticDataReshaped),clims);
        set(gca,'YDir','normal');
        ylabel("2nd Dimension Acquisitions"); xlabel("1st Dimension Acquisitions");
        hold on
        himage = imshow(Lrgb);
        himage.AlphaData = 0.8;

    end

end