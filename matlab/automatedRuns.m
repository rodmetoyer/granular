% Automated runs

close all; clear all; clc;

numParticles = 575;
numSims = 1;
increment = 500;
runTimes = NaN(1,numSims); % RunTimes do not include the postprocess step (i.e. the movie making time is not included)
simTime = 30.0;
tStep = 0.02;
movie = true;
data = true;
numPlot = NaN(1,numSims);
boxX = 2.5;
boxY = 1.5;
initXvel = 0.5;
initDisp = 0.8;
AngVel = 0.750;

for i=1:1:numSims
    close all;
    runTimes(i) = particle2dFuncFriction(numParticles,simTime,tStep,movie,data,boxX,boxY,initXvel,AngVel,initDisp)
    %runTimes(i) = particle2dFuncEven(numParticles,simTime,tStep,movie,data,boxX,boxY)
    numPlot(i) = numParticles;
    numParticles = numParticles + increment;
end

runTimes_hrs = runTimes/3600

figure
plot(numPlot,runTimes,'o');