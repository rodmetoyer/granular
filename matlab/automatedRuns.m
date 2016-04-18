% Automated runs

close all; clear all; clc;

numParticles = 537;
numSims = 3;
increment = 500;
runTimes = NaN(1,numSims); % RunTimes do not include the postprocess step (i.e. the movie making time is not included)
simTime = 30.0;
tStep = 0.01;
movie = true;
data = true;
numPlot = NaN(1,numSims);
boxX = 4.0;
boxY = 3.0;

for i=1:1:numSims
    close all;
    runTimes(i) = particle2dFunc(numParticles,simTime,tStep,movie,data,boxX,boxY)
    numPlot(i) = numParticles;
    numParticles = numParticles + increment;
end

runTimes_hrs = runTimes/3600

figure
plot(numPlot,runTimes,'o');