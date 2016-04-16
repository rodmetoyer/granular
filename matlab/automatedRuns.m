% Automated runs

close all; clear all; clc;

id = 1;
numParticles = 9;
numSims = 10;
increment = 21;
runTime = NaN(1,numSims);

for i=1:1:numSims
    close all;
    file = strcat('test',num2str(id),'num',num2str(numParticles));
    runTime(i) = particle2d(numParticles,file);
    numParticles = numParticles + increment;
    id = id + 1;
end