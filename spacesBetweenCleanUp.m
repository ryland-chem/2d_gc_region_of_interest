function [StartStop, newSS2, newBool] = spacesBetweenCleanUp(boolCutOff)

    %find where ROIs start and end
    startROI = find(diff(boolCutOff) == 1);
    endROI = find(diff(boolCutOff) == -1);

    %case for when the roi starts immediatley (in the first acquisiton)
    if boolCutOff(1) == 1

        startROI = [1; startROI];

    end

    %case for when the roi ends at the end
    if boolCutOff(size(boolCutOff, 1)) == 1

        endROI = [endROI; size(boolCutOff, 1)];

    end

    %calculate the start and stop numbers
    StartAndStop = vertcat(startROI, endROI);
    StartAndStop = sort(StartAndStop);
    StartStop = StartAndStop;
    newSS = StartStop;

    beginROI = 3;
    endROI = 2;

    while endROI < size(StartStop, 1)
    
        if StartStop(beginROI) - StartStop(endROI) <= 3
    
            newSS(endROI:beginROI) = 0;
    
        end
    
        beginROI = beginROI + 2;
        endROI = endROI + 2;
    
    end
    
    newSS(newSS == 0) = [];
    
    roiStart = 1;
    roiEnd = 2;
    
    newSS2 = newSS;
    
    while roiEnd <= size(newSS, 1)
    
        if newSS(roiEnd) - newSS(roiStart) <= 8
    
            newSS2(roiStart:roiEnd) = 0;
    
        end
    
        roiStart = roiStart + 2;
        roiEnd = roiEnd + 2;
    
    end
    
    newSS2(newSS2 == 0) = [];
    
    %make a zero vector
    zeroVect = zeros(1, size(boolCutOff, 1));

    j = 1;
    k = 2;

    while j < size(newSS2, 1)

        zeroVect(newSS2(j):newSS2(k)) = 1;

        j = j + 2;
        k = k + 2;

    end

    newBool = zeroVect;

end