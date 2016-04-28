% timeAvgVolumeFractions.m

clear all; close all; clc;
tic
boxX = 2.0;
boxY = 1.5;
numVboxes = round(boxY*20);
numHboxes = round(boxX*20);
numTotBoxes = numHboxes*numVboxes;
yInc = boxY/numVboxes;
xInc = boxX/numHboxes;

StartLine = 1;

% Get the data
filename = 'data/fric0.02_P575.txt';
dat = importdata(filename);
% Parse the data, each line is a time slice
H = zeros(length(dat(:,1))-StartLine+1,numVboxes,numHboxes);

disp(['Number of iterations: ',num2str(length(dat(:,1)))]);

itr = 1;
cntr = 1;
for lineNum = StartLine:1:length(dat(:,1))
    box = [0 0];    
    for i=1:1:(length(dat(lineNum,:))-1)/2
        box = [0 0];
        inc = 1;
        x = dat(lineNum,2*i);
        y = dat(lineNum,2*i+1);
        for j=1:1:numHboxes
            if isnan(x)
                break;
            end
            box(2) = j;
            finder = x-j*xInc;
            if finder < 0
                break;
            end
            if j == numHboxes
                inc = 0;
            end
        end
        for j=1:1:numVboxes
            if isnan(y)
                break;
            end
            box(1) = j;
            finder = y-j*yInc;
            if finder < 0
                break;
            end
            if j == numVboxes
                inc = 0;
            end
        end
        if(box(1)>0 && box(2)>0)
            H(itr,box(1),box(2)) = H(itr,box(1),box(2)) + inc;
        end
    end
    itr = itr + 1;
    cntr = cntr+1;
    if cntr > 200
        disp(itr);
        cntr = 1;
        toc
    end
end

Hnorm = mean(H);
Hnorm = Hnorm/max(max(Hnorm));

% Center of ij box is ((j-1)*xInc+xInc/2),((i-1)*yInc+yInc/2)
itr = 1;
%figure

figure
for j=1:1:numHboxes
    for i=1:1:numVboxes
        vFrac = Hnorm(1,i,j);            
        p = patch([(j-1)*xInc j*xInc j*xInc (j-1)*xInc],[(i-1)*yInc (i-1)*yInc i*yInc i*yInc],[vFrac 0 1-vFrac]);
        %hold on
    end
end
%hold off
axis([0 boxX 0 boxY]);
pbaspect([1 boxY/boxX 1]);
title(['Time Averaged Volume Fraction,  Num = ' num2str(sum(sum(H(100,:,:))))]);


totalTime = toc