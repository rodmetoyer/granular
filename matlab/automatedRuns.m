% Automated runs

close all; clear all; clc;

numParticles = 110;
numSims = 4;
increment = 53;
runTimes = NaN(1,numSims);
simTime = 30.0;
tStep = 0.05;
movie = true;
data = true;
numPlot = NaN(1,numSims);

for i=1:1:numSims
    close all;
    runTimes(i) = particle2dFunc(numParticles,simTime,tStep,movie,data)
    numPlot(i) = numParticles;
    numParticles = numParticles + increment;
end

figure
plot(numPlot,runTimes,'o');