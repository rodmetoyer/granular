% Automated runs

close all; clear all; clc;

id = 1;
numParticles = 100;
numSims = 10;
increment = 100;
runTime = NaN(1,numSims);

for i=1:1:numSims
    close all;
    file = strcat('autotestID',num2str(id),'numPs');
    runTime(i) = particle2d(numParticles,file)
    numParticles = numParticles + increment;
    id = id + 1;
end