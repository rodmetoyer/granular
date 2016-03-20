

clear all; close all; clc; clear class;

for i=1:1:100
    pointParticle(i) = particle();
end

for i=1:1:1
    sphereParticle(i) = sphere([0.0;0.0;0.0],[0.0;0.0;0.0],[0.0;0.0;0.0;0.0],[0.0;0.0;0.0;0.0],100);
end