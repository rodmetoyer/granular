classdef particle
    %Particle Particle is the highest level particle class
    %   The Particle class is the highest level particle class. Objects of
    %   the Particle class are point objects. Subordinate classes may have
    %   geometry. Particles must have mass.
    
    % Public properties - anyone can change these
    properties (GetAccess=public,SetAccess=public)
        force = [0.0; 0.0; 0.0]
    end
    
    properties (GetAccess=public,SetAccess=protected)
        % Scalar
        mass = 0.1
        % Vector
        position = [0.0; 0.0; 0.0]
        velocity = [0.0; 0.0; 0.0]
    end
    
    % Protected properties - only this class
    properties (Dependent)
        acceleration
    end
    
    methods
        % ctor
        function thisParticle = particle(mass,position,velocity)
            if nargin > 0
                thisParticle.mass = mass;
                thisParticle.position = position;
                thisParticle.velocity = velocity;
            end
        end
        
        function acceleration=get.acceleration(obj)
            acceleration = obj.force/obj.mass;
        end
        
        function obj = advanceStates(obj,dt)
            obj.position = obj.position + obj.velocity*dt;
            obj.velocity = obj.velocity + obj.acceleration*dt;
        end
        
    end    
end

%% Helper funcitons
