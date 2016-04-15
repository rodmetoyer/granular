classdef sphereObstruction
    %polygon is a point-defined obstruction obstructing the flow
    
    properties
        midpoint
        radius
    end
    
    methods
        %ctor
        function obj = sphereObstruction(midpoint,radius)
            if(nargin > 0)
                obj.midpoint = midpoint;
                obj.radius = radius;
            end
        end
        
        
    end
end