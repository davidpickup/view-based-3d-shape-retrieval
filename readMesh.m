function mesh = readMesh(fname)
% mesh = readMesh(fname)
% Read in a mesh from a file.
% Variables:
% mesh - output mesh.
% fname - input file name.
%
% David Pickup 2013

[vertex,faces] = read_mesh(fname);

S = size(vertex);
if S(1) > S(2)
    mesh.X = vertex(:,1);
    mesh.Y = vertex(:,2);
    mesh.Z = vertex(:,3);

    mesh.TRIV = faces;
else
    mesh.X = vertex(1,:)';
    mesh.Y = vertex(2,:)';
    mesh.Z = vertex(3,:)';

    mesh.TRIV = faces';
end

return;