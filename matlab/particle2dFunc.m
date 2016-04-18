function runTime = particle2dFunc(totalNumParticles,simTime,timeStep,doYouWantMovie,doYouWantDataFile,boxx,boxy)

tic;
% Simulation properties
%totalNumParticles = 100; % This is the total number of particles you want
%simTime = 10;
%timeStep = 0.01;
numSteps = round(simTime/timeStep)+1;
%doYouWantMovie = false;
%doYouWantDataFile = false;
test = strcat('auto_small_ts',num2str(timeStep),'_P',num2str(totalNumParticles));
movieFile = strcat(test,'.avi');
dataFile = strcat(test,'.txt');
frameRate = 20;         % frame rate of the movie file
speedReduction = 1.2;   % reduce the frame rate by a constant value
plotRadiusScaler = 85;

% World properties: the box always starts at (0,0)
%boxx = 10.0; % x-coord of box upper border
%boxy = 6.0; % y-coord of box upper border

% Particle properties
% Same for all particles
xVelocity = 5.0;
mass    = 0.1;
radius  = 0.025;
spring  = 25.0;
damp    = 1.0;
fricCo1 = 1.0;
fricCo2 = 1.0;

% Add obstruction
obRadius = 0.5;
obstruction = sphereObstruction([boxx-2.0*obRadius,boxy/2],obRadius);
   
% Show the box and the obstruction(s) without particles
figure
plot(obstruction.midpoint(1),obstruction.midpoint(2),'o','MarkerSize',obstruction.radius*plotRadiusScaler','MarkerEdgeColor','none','MarkerFaceColor','r');
axis([0 boxx 0 boxy]);
pbaspect([1 boxy/boxx 1]);

% Initialize the simulation
numParticles = 0; % This changes, this is not the final total
numArray = NaN(1,numSteps);
spawnPeriod = simTime/10; % Rough spawning period
spawnRate_perStep = max(1,round(totalNumParticles/((spawnPeriod/timeStep)))); % Number to spawn per time step
shortage = totalNumParticles - spawnRate_perStep*timeStep*spawnPeriod;
spawnPeriod = spawnPeriod+shortage/spawnRate_perStep;
check = round(spawnPeriod/timeStep);
itr = 1;
time = NaN(1,numSteps);
force = [0;0;0]; %Spawn force
ParticleArray = Sphere2d(totalNumParticles);
xCoord = NaN(1,totalNumParticles);
yCoord = NaN(1,totalNumParticles);

%% Run the simulation
for j=1:1:numSteps
    % Spawn particles for the first few seconds to get a steady density
    if(j<check+1)
        if(numParticles < totalNumParticles)
            %Particles spawn at the x=0 line.
            xxtra = zeros(1,spawnRate_perStep) + rand(1,spawnRate_perStep)*xVelocity*timeStep;
            yxtra = randperm(100*totalNumParticles,spawnRate_perStep);
            yxtra = 2*radius + yxtra/(100*totalNumParticles)*(boxy-2*radius);
            xVecs = [xxtra;yxtra;zeros(1,spawnRate_perStep)];
            vVecs = xVelocity*[ones(1,spawnRate_perStep);zeros(2,spawnRate_perStep)];
            for i=1:1:spawnRate_perStep
                ParticleArray(numParticles + 1) = Sphere2d(mass,radius,spring,damp,fricCo1,fricCo2,xVecs(:,i),vVecs(:,i),force);
                numParticles = numParticles + 1;
            end
        end
    end
    
    for i=1:numParticles
        % Compute forces
        ParticleArray(i).force = [0;0;0];
        if (ParticleArray(i).position(1)+ParticleArray(i).radius) > (obstruction.midpoint(1)-obstruction.radius)
            if ParticleArray(i).position(2)+ParticleArray(i).radius > obstruction.midpoint(2)-obstruction.radius && ParticleArray(i).position(2)+ParticleArray(i).radius < obstruction.midpoint(2)+obstruction.radius
                % force
                xComp = (obstruction.midpoint(1)-ParticleArray(i).position(1));
                yComp = (obstruction.midpoint(2)-ParticleArray(i).position(2));
                distance = sqrt(xComp^2+yComp^2);
                maxDist = obstruction.radius+ParticleArray(i).radius;
                if distance < maxDist
                    % force is spring and damper
                    xCompNorm = xComp/distance;
                    yCompNorm = yComp/distance;
                    % Speed of the particle is the relative speed for the
                    % fixed obstruction
                    distdot = xComp/distance*ParticleArray(i).velocity(1)+yComp/distance*ParticleArray(i).velocity(2);
                    ParticleArray(i).force = -[xCompNorm*(ParticleArray(i).spring*(maxDist-distance)+ParticleArray(i).damp*distdot);yCompNorm*(ParticleArray(i).spring*(maxDist-distance)+ParticleArray(i).damp*distdot);0];
                else
                    ParticleArray(i).force = [0;0;0];
                end
            end
        end
        % Loop through all of the other particles - bad way to do this!!!
        for k=1:numParticles
            if i==k
                continue; % I don't play with myself
            end
            % Is he close to me?
            xCompPart = (ParticleArray(k).position(1)-ParticleArray(i).position(1));
            yCompPart = (ParticleArray(k).position(2)-ParticleArray(i).position(2));
            distPart = sqrt(xCompPart^2+yCompPart^2);
            maxDistPart = ParticleArray(i).radius + ParticleArray(k).radius;
            if distPart < maxDistPart
                % Yes - then he's going to add to my force maybe
                xCompPartNorm = xCompPart/distPart;
                yCompPartNorm = yCompPart/distPart;
                % Now relative matters
                relVx = ParticleArray(k).velocity(1) - ParticleArray(i).velocity(1);
                rexVy = ParticleArray(k).velocity(2) - ParticleArray(i).velocity(2);
                distPartDot = xCompPartNorm*relVx+yCompPartNorm*rexVy;
                ParticleArray(i).force = ParticleArray(i).force...
                    -[xCompPartNorm*(ParticleArray(i).spring*(maxDistPart-distPart)-ParticleArray(i).damp*distPartDot);...
                    yCompPartNorm*(ParticleArray(i).spring*(maxDistPart-distPart)-ParticleArray(i).damp*distPartDot);...
                    0];
            end
        end        
        
        xCoord(j,i) = ParticleArray(i).position(1);
        yCoord(j,i) = ParticleArray(i).position(2);
        
        % Advance the particle states
        %ParticleArray(i) = ParticleArray(i).computeAcceleration;
        ParticleArray(i).acceleration = ParticleArray(i).force/ParticleArray(i).mass;
        %ParticleArray(i) = ParticleArray(i).advanceStatesEuler(timeStep);
        ParticleArray(i).position = ParticleArray(i).position + ParticleArray(i).velocity*timeStep;
        ParticleArray(i).velocity = ParticleArray(i).velocity + ParticleArray(i).acceleration*timeStep;
        
        % If a particle is outside the box, respawn on the zero line
        if(xCoord(j,i)>boxx || yCoord(j,i)>boxy || xCoord(j,i)*yCoord(j,i)<0.0)
            xCoord(j,i) =  rand*xVelocity*timeStep;
            %rng(81); % Set the seed again to get different numbers
            yCoord(j,i) = randperm(100*totalNumParticles,1);
            yCoord(j,i) = 2*radius + yCoord(j,i)/(100*totalNumParticles)*(boxy-2*radius);
            xVecs = [xCoord(j,i);yCoord(j,i);0];
            vVecs = xVelocity*[1;0;0];
            ParticleArray(i) = Sphere2d(mass,radius,spring,damp,fricCo1,fricCo2,xVecs,vVecs,force);
        end
        
    end
    % Advance time
    time(j) = j*timeStep-timeStep; 
    % Note how many particles are currently in the system
    numArray(itr) = numParticles;
    itr = itr+1;
end
runTime = toc;
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

if doYouWantDataFile
    if exist('data','dir') == 0
        mkdir('data');
    end
    csvwrite(strcat('data/',dataFile),data);
end

% Make a movie of the motion
plotInterval = 1/(frameRate*timeStep);
plotInterval = floor(plotInterval/speedReduction);
% Note that the movie may be slightly slower than real time depending on
% the time step used. This is great for really small time steps.
itr = 1;
if doYouWantMovie
    figure;
    movegui(gcf);
    clear Mov; % Just in case
    for i=1:plotInterval:numSteps
        for j = 1:1:numParticles
            plot(xCoord(i,j),yCoord(i,j),'o','MarkerSize',radius*plotRadiusScaler,'MarkerEdgeColor','none','MarkerFaceColor',[0.1 j/numParticles 1-j/numParticles]);
            hold on
        end
        plot(obstruction.midpoint(1),obstruction.midpoint(2),'o','MarkerSize',obstruction.radius*plotRadiusScaler','MarkerEdgeColor','none','MarkerFaceColor','r');
        hold off                        % turn off the plot
        axis([0 boxx 0 boxy]);       % set axis 
        pbaspect([1 boxy/boxx 1]);
        grid off;                        % turn on the grid
        title(['Time = ', num2str(time(i),'%4.2f'), ' seconds:  Num = ' num2str(numArray(i))]); % put current time in the title
        Mov(itr) = getframe(gcf);         % get the frame and compile it into the movie file
        itr = itr + 1;
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

