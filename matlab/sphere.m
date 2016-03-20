classdef sphere < particle
    %sphereParticle This class is derived from the particle class
    %   The sphereParticle class is derived from the particle class. It has
    %   some additional properties because it has geometry.
    
    properties
        q = [0.0;0.0;0.0;0.0]
        qd = [0.0;0.0;0.0;0.0]
    end
    
    properties (Access = protected)
       radius
       k
       c
    end
    
    methods
        %ctor
        function this = sphereParticle(position,velocity,q,qd)
            if nargin > 0
                this.position = position;
                this.velocty = velocity;
                this.q = q;
                this.qd = qd;
            end
        end
    end
    
end

