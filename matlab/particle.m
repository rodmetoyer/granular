classdef particle
    %Particle Particle is the highest level particle class
    %   The Particle class is the highest level particle class. Objects of
    %   the Particle class are point objects. Subordinate classes may have
    %   geometry.
    
    % Public properties - anyone can change these
    properties
        position = [0.0; 0.0; 0.0]
        velocity = [0.0; 0.0; 0.0]
    end
    
    % Protected properties - only this class
    properties (Access = protected)
        
    end
    
    methods
        % ctor
        function this = particle(position,velocity)
            if nargin > 0
                this.position = position;
                this.velocity = velocity;
            end
        end
        
    end
    
end

