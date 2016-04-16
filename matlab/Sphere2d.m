classdef Sphere2d
    %Sphere2d is a two dimensional spherical particle class
    
    % Public properties - anyone can change these
    properties
        mass
        radius
        spring
        damp
        frictionCo1
        frictionCo2
        position
        velocity
        acceleration
        force
    end
    
%     properties (Dependent)
%         acceleration
%     end
    
    methods
        % ctor
        function obj = Sphere2d(mass,radius,spring,damp,fricCo1,fricCo2,xVec,vVec,force)            
            if(nargin > 0)
                if(nargin == 1)
                obj = repmat(obj,1,mass);
                else
                obj.mass        = mass;                
                obj.radius      = radius;
                obj.spring      = spring;
                obj.damp        = damp;
                obj.frictionCo1 = fricCo1;
                obj.frictionCo2 = fricCo2;                
                obj.position    = xVec;
                obj.velocity    = vVec;
                obj.force       = force;
                end
            end
        end
        
        function obj = advanceStatesEuler(obj,timeStep)
            obj.position = obj.position + obj.velocity*timeStep;
            obj.velocity = obj.velocity + obj.acceleration*timeStep;
        end
        
        function obj = computeForce(obj,others)
            % force on obj due to collision with others
            
        end
        
        function obj = computeAcceleration(obj)
            obj.acceleration = obj.force/obj.mass;
        end
        
    end    
end
