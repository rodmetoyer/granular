% particle2dSim.m
% Main simulation file for the 2d particle simulation

%% Set-up the simulation
clear all; close all; clc;

simTime = 10;
timeStep = 0.1;
numSteps = round(simTime/timeStep)+1;
doYouWantMovie = false;
movieFile = 'test.avi';
frameRate = 20;         % frame rate of the movie file
speedReduction = 1.0;   % reduce the frame rate by a constant value

% The box: always starts at (0,0)
boxx = 10.0; % x-coord of box upper border
boxy = 10.0; % y-coord of box upper border

time = 0.0;

% Particles
% Same for all particles
numParticles = 100;
mass    = 1.0;
radius  = 0.1;
spring  = 1.0;
damp    = 1.0;
fricCo1 = 1.0;
fricCo2 = 1.0;
% Random distribution in the first 1/3 of the box
%rng(77); % Set the seed. Built-in rng good enough
xCoord = randperm(100*numParticles,numParticles);
xCoord = 2*radius + xCoord/(100*numParticles)*(boxx/3-2*radius);
%rng(81); % Set the seed again to get different numbers
yCoord = randperm(100*numParticles,numParticles);
yCoord = 2*radius + yCoord/(100*numParticles)*(boxy-2*radius);
% No need to scatter the angular position
thetas = zeros(1,numParticles);

xVecs = [xCoord;yCoord;thetas];
% Velocities: FOr now all in the x direction. TODO make them stochastic
vx0 = 2.0;
vVecs = vx0*[ones(1,numParticles);zeros(2,numParticles)];

% Make the particles
for i=1:1:numParticles
    ParticleArray(i) = Sphere2d(mass,radius,spring,damp,fricCo1,fricCo2,xVecs(:,i),vVecs(:,i));
end
    
% Show the box and the particles
figure
plot(ParticleArray(1).position(1),ParticleArray(1).position(2),'o','MarkerEdgeColor','none','MarkerFaceColor','k');
hold on
for i=2:1:numParticles
    plot(ParticleArray(i).position(1),ParticleArray(i).position(2),'o','MarkerEdgeColor','none','MarkerFaceColor','k');
end
axis([0 boxx 0 boxy])
hold off
%% Run the simulation

% compute the forces/acceleration on all of the particles and advance time
for j=2:1:numSteps
    
    for i=1:1:numParticles
        ParticleArray(i).force = 0;
        
        ParticleArray(i) = ParticleArray(i).computeAcceleration;
        ParticleArray(i) = ParticleArray(i).advanceStatesEuler(timeStep);
        xCoord(j,i) = ParticleArray(i).position(1);
        yCoord(j,i) = ParticleArray(i).position(2);
    end
    time(j) = time(j-1) + timeStep;
end

%% Postprocess the data
% Make a movie of the motion
if doYouWantMovie
    figure;
    movegui(gcf);
    clear Mov; % Just in case
    n = 0;
    size = 10;
    for i=1:1:numSteps
        for j = 1:1:numParticles
            plot(xCoord(i,j),yCoord(i,j),'o','MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor','k');
            hold on
        end
        hold off                        % turn off the plot
        axis([0 boxx 0 boxy]);       % set axis 
        grid off;                        % turn on the grid
        title(['Time = ', num2str(time(i)), ' seconds']); % put current time in the title
        Mov(i) = getframe(gcf);         % get the frame and compile it into the movie file
    end
    writerObj = VideoWriter(movieFile); % write the movie to a file
    writerObj.FrameRate = frameRate/speedReduction; writerObj.Quality = 100; % optional
    open(writerObj); writeVideo(writerObj,Mov); close(writerObj);
    
    % Move the movie to the bin
    % TODO(Rodney) gitignore the bin
    if exist('bin','dir') == 0
        mkdir('bin');
    end
    movefile(movieFile,['bin/' movieFile]);
end