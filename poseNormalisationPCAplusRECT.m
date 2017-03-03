function mesh = poseNormalisationPCAplusRECT(mesh)
% mesh = poseNormalisationPCAplusRECT(mesh)
% Normalises the pose of a mesh using PCA and rectilinearity.
% Variables:
% mesh - mesh to be normalised.
%
% David Pickup 2013

%% Translate mesh so its centroid is at the origin.
mesh.X = mesh.X - mean(mesh.X);
mesh.Y = mesh.Y - mean(mesh.Y);
mesh.Z = mesh.Z - mean(mesh.Z);

%% Scale mesh so the largest radial distance is equal to 1.
distance = sqrt(mesh.X.^2 + mesh.Y.^2 + mesh.Z.^2);
scale = 1/max(distance);
mesh.X = mesh.X .* scale;
mesh.Y = mesh.Y .* scale;
mesh.Z = mesh.Z .* scale;
clear distance;

%% Use PCA to rotate mesh so its principle componants are aligned with
% the X, Y and Z axes.
try
    meshPCA = mesh;
    [coeff,score] = princomp([mesh.X,mesh.Y,mesh.Z]);
    clear coeff;
    meshPCA.X = score(:,1);
    meshPCA.Y = score(:,2);
    meshPCA.Z = score(:,3);
catch
    error('Error computing PCA embedding.');
end

%% Use rectilinearity to rotate mesh so its orientation give the highest
% rectilinearity score.
meshRect = mesh;
[r angles] = rectilinearity3D(meshRect);
% Get angles.
alpha = angles(1);
beta = angles(2);
gamma = angles(3);
% Create X rotation matrix.
Rx = [	1 0 0;...
	0 cos(alpha) -sin(alpha);...
	0 sin(alpha) cos(alpha)];
% Create Y rotation matrix.
Ry = [	cos(beta) 0 sin(beta);...
	0 1 0;...
	-sin(beta) 0 cos(beta)];
% Create Z rotation matrix.
Rz = [	cos(gamma) -sin(gamma) 0;...
	sin(gamma) cos(gamma) 0;...
	0 0 1];
% Get mesh vertices.
v = zeros(3, numel(meshRect.X));
v(1,:) = meshRect.X(:);
v(2,:) = meshRect.Y(:);
v(3,:) = meshRect.Z(:);
% Rotate mesh vertices.
v = Rx*Ry*Rz*v;
% Put vertices back into mesh.
meshRect.X(:) = v(1,:);
meshRect.Y(:) = v(2,:);
meshRect.Z(:) = v(3,:);

%% Render silhoettes of the two results and choose the result with the
% fewest black pixels.
vs = [  1 0 0;...
        0 1 0;...
        0 0 1   ];
imagesPCA = depthImages(meshPCA, vs, [256 256]);
imagesRect = depthImages(meshRect, vs, [256 256]);
% imagesPCA = renderDepthImages(meshPCA.TRIV(:,1), meshPCA.TRIV(:,2), ...
%     meshPCA.TRIV(:,3), meshPCA.X, meshPCA.Y, meshPCA.Z, vs(:,1), ...
%     vs(:,2), vs(:,3), 256);
% imagesRect = renderDepthImages(meshRect.TRIV(:,1), meshRect.TRIV(:,2), ...
%     meshRect.TRIV(:,3), meshRect.X, meshRect.Y, meshRect.Z, vs(:,1), ...
%     vs(:,2), vs(:,3), 256);

pcaScore = 0;
rectScore = 0;
for i = 1:3
    pcaScore = pcaScore + numel(find(imagesPCA{i}<1));
    rectScore = rectScore + numel(find(imagesRect{i}<1));
end

if pcaScore < rectScore
    mesh = meshPCA;
    'PCA'
else
    mesh = meshRect;
    'Rectilinearity'
end

return;