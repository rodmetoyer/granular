classdef polygon
    %polygon is a point-defined obstruction obstructing the flow
    
    properties
        points
    end
    
    methods
        %ctor
        function obj = polygon(points)
            if(nargin > 0)
                obj.points = points;                
            end
        end
        
        
    end
end