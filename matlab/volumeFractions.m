% volumeFractions.m

clear all; close all; clc;
tic
boxX = 4.0;
boxY = 3.0;
numVboxes = round(boxY*10);
numHboxes = round(boxX*10);
numTotBoxes = numHboxes*numVboxes;
yInc = boxY/numVboxes;
xInc = boxX/numHboxes;
%H = zeros(numVboxes,numHboxes);

StartLine = 600;
framerate = 60;
movieFile = 'bin/testVol3.avi';

% Get the data
filename = 'data/auto_small_ts0.01_P600.txt';
dat = importdata(filename);
% Parse the data, each line is a time slice
H = zeros(length(dat)-StartLine,numVboxes,numHboxes);

itr = 1;
for lineNum = StartLine:1:length(dat)
    time(itr) = dat(lineNum,1);
    box = [0 0];
    % H(itr,:,:) = zeros(numVboxes,numHboxes);
    for i=1:1:(length(dat(lineNum,:))-1)/2
        box = [0 0];
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
        end
        if(box(1)>0 && box(2)>0)
            H(itr,box(1),box(2)) = H(itr,box(1),box(2)) + 1;
        end
    end

    %Hnorm = H/max(max(H));

    % Center of ij box is ((j-1)*xInc+xInc/2),((i-1)*yInc+yInc/2)
%     for j=1:1:numHboxes
%         for i=1:1:numVboxes
%             vFrac = Hnorm(i,j);
%             patch([(j-1)*xInc j*xInc j*xInc (j-1)*xInc],[(i-1)*yInc (i-1)*yInc i*yInc i*yInc],[vFrac 0 1-vFrac]);
%         end
%     end
%     axis([0 boxX 0 boxY]);
%     pbaspect([1 boxY/boxX 1]);
%     title(['Time = ', num2str(time,'%4.2f'), ' seconds:  Num = ' num2str(sum(sum(H)))]);
%     
%     Mov(itr) = getframe(gcf);         % get the frame and compile it into the movie file
    itr = itr + 1;
end

Hnorm = H/max(max(max(H)));

% Center of ij box is ((j-1)*xInc+xInc/2),((i-1)*yInc+yInc/2)
itr = 1;
%figure
for k = StartLine:1:length(dat)
    figure
    for j=1:1:numHboxes
        for i=1:1:numVboxes
            vFrac = Hnorm(k,i,j);            
            p = patch([(j-1)*xInc j*xInc j*xInc (j-1)*xInc],[(i-1)*yInc (i-1)*yInc i*yInc i*yInc],[vFrac 0 1-vFrac]);
            %hold on
        end
    end
    %hold off
    axis([0 boxX 0 boxY]);
    pbaspect([1 boxY/boxX 1]);
    title(['Time = ', num2str(time(itr),'%4.2f'), ' seconds:  Num = ' num2str(sum(sum(H(itr,:,:))))]);
    
    Mov(itr) = getframe(gcf);         % get the frame and compile it into the movie file
    itr = itr + 1;
    close(gcf)
end

writerObj = VideoWriter(movieFile); % write the movie to a file
writerObj.FrameRate = framerate; writerObj.Quality = 100; % optional
open(writerObj); writeVideo(writerObj,Mov); close(writerObj);
totalTime = toc