classdef sphere < particle
    %sphereParticle This class is derived from the particle class
    %   The sphereParticle class is derived from the particle class. It has
    %   some additional properties because it has geometry.
        
    properties (GetAccess=public,SetAccess=protected)
       radius = 1.0
       q = [0.0;0.0;0.0;0.0]
       qd = [0.0;0.0;0.0;0.0]
       k = 1.0
       c = 1.0
    end
    
    methods
        %ctor
        function sphereObj = sphere(mass,radius,position,velocity,q,qd)
         if nargin == 0
            mass = 0.1;
            position = [0.0; 0.0; 0.0];
            velocity = [0.0; 0.0; 0.0];
         end
         sphereObj@particle(mass,position,velocity);
         if nargin > 0 % Use value if provided
            sphereObj.radius = radius;
            sphereObj.q = q;
            sphereObj.qd = qd;
         end
        end
    
    end
    
end

