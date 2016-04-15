% particle2dSim.m
% Main simulation file for the 2d particle simulation

%% Set-up the simulation
clear all; close all; clc;

simTime = 30;
timeStep = 0.1;
numSteps = round(simTime/timeStep)+1;
doYouWantMovie = false;
test = 'test200';
movieFile = strcat(test,'.avi');
dataFile = strcat(test,'.txt');
frameRate = 20;         % frame rate of the movie file
speedReduction = 1.0;   % reduce the frame rate by a constant value

% The box: always starts at (0,0)
boxx = 10.0; % x-coord of box upper border
boxy = 10.0; % y-coord of box upper border

time = 0.0;

% Particles
% Same for all particles
totalNumParticles = 200; % This is roughly the total number you want
numParticles = 1; % This changes, this is not the final total
num = numParticles;
spawnPeriod = 5; % Spawning period
spawnRate = round(totalNumParticles/(round(spawnPeriod/timeStep))); % Number to spawn per time step
mass    = 1.0;
radius  = 0.1;
spring  = 1.0;
damp    = 1.0;
fricCo1 = 1.0;
fricCo2 = 1.0;
% Random distribution in the first 1/5 of the box
rng(77); % Set the seed. Built-in rng good enough
xCoord = randperm(100*numParticles,numParticles);
xCoord = 2*radius + xCoord/(100*numParticles)*(boxx/5-2*radius);
rng(81); % Set the seed again to get different numbers
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
itr = 1;
for j=1:1:numSteps
    % Spawn particles for the first few seconds to get a steady density
    if(j<(round(spawnPeriod/timeStep)+1))
        xxtra = randperm(100*numParticles,spawnRate);
        % Respawn somewhere in the first 1/5 of the box
        xxtra = 2*radius + xxtra/(100*numParticles)*(boxx/5-2*radius);
        yxtra = randperm(100*numParticles,spawnRate);
        yxtra = 2*radius + yxtra/(100*numParticles)*(boxy-2*radius);
        xVecs = [xxtra;yxtra;zeros(1,spawnRate)];
        vVecs = vx0*[ones(1,spawnRate);zeros(2,spawnRate)];
        for i=1:1:spawnRate
            ParticleArray(numParticles + 1) = Sphere2d(mass,radius,spring,damp,fricCo1,fricCo2,xVecs(:,i),vVecs(:,i));
            numParticles = numParticles + 1;
            num(itr) = numParticles;
        end
    end
    num(itr) = numParticles;
    itr = itr+1;
    for i=1:1:numParticles
        ParticleArray(i).force = 0;
        
        ParticleArray(i) = ParticleArray(i).computeAcceleration;
        ParticleArray(i) = ParticleArray(i).advanceStatesEuler(timeStep);
        xCoord(j,i) = ParticleArray(i).position(1);
        yCoord(j,i) = ParticleArray(i).position(2);
        % If it is outside the box, respawn
        if(xCoord(j,i)>boxx || yCoord(j,i)>boxy || xCoord(j,i)*yCoord(j,i)<0.0)
            xCoord(j,i) = randperm(100*numParticles,1);
            % Respawn somewhere in the first 1/5 of the box
            xCoord(j,i) = 2*radius + xCoord(j,i)/(100*numParticles)*(boxx/5-2*radius);
            %rng(81); % Set the seed again to get different numbers
            yCoord(j,i) = randperm(100*numParticles,1);
            yCoord(j,i) = 2*radius + yCoord(j,i)/(100*numParticles)*(boxy-2*radius);
            xVecs = [xCoord(j,i);yCoord(j,i);0];
            vVecs = vx0*[1;0;0];
            ParticleArray(i) = Sphere2d(mass,radius,spring,damp,fricCo1,fricCo2,xVecs,vVecs);
        end
    end
    time(j) = j*timeStep;
    
end

%% Postprocess the data
% Make a file of the positional data
% Time x1 y1 x2 y2 x3 y3 ...
for i=1:1:numSteps
    for j=1:1:numParticles
        data(i,2*j-1) = xCoord(i,j);
        data(i,2*j) = yCoord(i,j);
    end
end
data = [time',data];
csvwrite(dataFile,data);

% Make a movie of the motion
if doYouWantMovie
    figure;
    movegui(gcf);
    clear Mov; % Just in case
    for i=1:1:numSteps
        for j = 1:1:numParticles
            plot(xCoord(i,j),yCoord(i,j),'o','MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor','k');
            hold on
        end
        hold off                        % turn off the plot
        axis([0 boxx 0 boxy]);       % set axis 
        grid off;                        % turn on the grid
        title(['Time = ', num2str(time(i)), ' seconds:  Num = ' num2str(num(i))]); % put current time in the title
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