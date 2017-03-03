function [r,angles] = rectilinearity3D(mesh, ~)
% [r,angles] = rectilinearity3D(mesh, unoptimised(optional))
% Computed the retilinearity measure of a 3D mesh.
% Variables:
% r - output rectilinearity.
% angles - best angles that optimise rectilinearity measure.
% mesh - input mesh.
% unoptimised (optional) - if specified then angles are not optimised.
%
% David Pickup

if nargin > 1
    angles = [0 0 0];
    minP = projectedSurfaceAreaMex(angles,mesh.TRIV,mesh.X,mesh.Y,mesh.Z);
else
    % Initialise minimum projected perimeter.
    minP = Inf;
    finalangles = [0 0 0];

    % Perform angle optimisation five times to avoid a local minimum.
    for i = 1:30
        angles = fminsearch(@(x) projectedSurfaceAreaMex(x,mesh.TRIV,mesh.X,mesh.Y,mesh.Z), [rand()*2*pi, rand()*2*pi, rand()*2*pi]);
        P = projectedSurfaceAreaMex(angles,mesh.TRIV,mesh.X,mesh.Y,mesh.Z);
        if P < minP
            minP = P;
            finalangles = angles;
        end
    end
    angles = finalangles;
end

% Initialise the total surface area of the mesh.
S = 0;

% Iterate through all triangles.
for i = 1:size(mesh.TRIV,1)
    % Get current triangle.
    tri = mesh.TRIV(i,:);

    % Get verteces for current triangle.
    v1 = [mesh.X(tri(1)) mesh.Y(tri(1)) mesh.Z(tri(1))];
    v2 = [mesh.X(tri(2)) mesh.Y(tri(2)) mesh.Z(tri(2))];
    v3 = [mesh.X(tri(3)) mesh.Y(tri(3)) mesh.Z(tri(3))];

    % Get edge lengths of current triangle.
    a = pdist([v1; v2]);
    b = pdist([v1; v3]);
    c = pdist([v2; v3]);

    % Calculate area of triangle using Heron's formula and add to total surface area.
    S = S + 0.25*sqrt((a^2 + b^2 + c^2)^2 - 2*(a^4 + b^4 + c^4));
end

% Compute final rectilinearity measure.
r = 3*((S/minP) - (2/3));

return;
