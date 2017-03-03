function [edgeTables, vertexTables, vs, tris] = geodesicSphere(Nd)
% [edgeTables, vertexTables, vs, tris] = geodesicSphere(Nd)
% Computed a geodesic sphere by subdividing a octahedron.
% Variables:
% edgeTables - edge table for each subdivision layer.
% vertexTables - vertex table for each subdivision layer.
% vs - vertices for final geodesic sphere.
% tris - triangles for final geodesic sphere.
% Nd - number of times to subdivide the octahedron.
%
% David Pickup

% Initialise vertex tables.
vertexTables = {};

% Initialise octahedron vertices.
vs = [1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1];
size(vs)

% Initialise octahedron edges.
edgeTables = {[ 1,3;...
                5,3;...
                1,5;...
                6,3;...
                1,6;...
                1,4;...
                5,4;...
                6,4;...
                2,3;...
                2,5;...
                2,6;...
                2,4]};

% Initialise octahedron triangles.
tris = [    1,2,3;
            6,7,3;
            9,2,10;
            12,7,10;
            5,4,1;
            11,4,9;
            6,8,5;
            12,8,11];

% Initialise current vertex number.
vNum = 7;

% Iterate through the number of subdivisions.
for i = 1:Nd
    % Initialise current vertex table.
    vertexTables = {vertexTables{:}, zeros(size(vs,1), size(vs,1))};
    vertexTables{i}(:) = -1;
    
    % Initialise a new edgeTable and triangles.
    newTris = zeros(4*size(tris,1),3);
    es = edgeTables{i};
    esNew = zeros((size(newTris,1)*3)./2,2);
    
    % Remove duplicate edges just for vertex creation.
    esTmp = unique(es,'rows');
    
    % Iterate through all edges.
    for j = 1:size(esTmp,1)
        % Create new vertices at the centre of the edge
        vs(vNum,:) = (vs(esTmp(j,1),:) + vs(esTmp(j,2),:))./2;
        vs(vNum,:) = vs(vNum,:) ./ norm(vs(vNum,:));
            
        % Add new vertices to vertex table.
        vertexTables{i}(esTmp(j,1),esTmp(j,2)) = vNum;
        vertexTables{i}(esTmp(j,2),esTmp(j,1)) = vNum;
        
        % Increment number of vertices.
        vNum = vNum + 1;
    end
    clear esTmp;
    
    % Initialise triangle and edge counters.
    tNum = 1;
    eNum = 1;
    
    % Iterate through all triangles.
    for j = 1:size(tris,1)
        % Get current triangle and vertex table.
        tri = tris(j,:);
        vTab = vertexTables{i};

        % Get the three vertices of the triangle.
        v1 = vTab(es(tri(1),1),es(tri(1),2));
        v2 = vTab(es(tri(2),1),es(tri(2),2));
        v3 = vTab(es(tri(3),1),es(tri(3),2));
        
        % Add new edges to edge table.
        esNew(eNum,:) = [v1, es(tri(1),1)];
        esNew(eNum+1,:) = [v1, es(tri(1),2)];
        esNew(eNum+2,:) = [v2, es(tri(2),1)];
        esNew(eNum+3,:) = [v2, es(tri(2),2)];
        esNew(eNum+4,:) = [v3, es(tri(3),1)];
        esNew(eNum+5,:) = [v3, es(tri(3),2)];
        esNew(eNum+6,:) = [v1, v2];
        esNew(eNum+7,:) = [v2, v3];
        esNew(eNum+8,:) = [v3, v1];
        
        % Add new triangles to triangle table.
        if es(tri(1),1) == es(tri(2),1)
            newTris(tNum,:) = [eNum,eNum+6,eNum+2];
        elseif es(tri(1),1) == es(tri(2),2)
            newTris(tNum,:) = [eNum,eNum+6,eNum+3];
        elseif es(tri(1),1) == es(tri(3),1)
            newTris(tNum,:) = [eNum,eNum+8,eNum+4];
        elseif es(tri(1),1) == es(tri(3),2)
            newTris(tNum,:) = [eNum,eNum+8,eNum+5];
        end
        
        if es(tri(1),2) == es(tri(2),1)
            newTris(tNum+1,:) = [eNum+1,eNum+6,eNum+2];
        elseif es(tri(1),2) == es(tri(2),2)
            newTris(tNum+1,:) = [eNum+1,eNum+6,eNum+3];
        elseif es(tri(1),2) == es(tri(3),1)
            newTris(tNum+1,:) = [eNum+1,eNum+8,eNum+4];
        elseif es(tri(1),2) == es(tri(3),2)
            newTris(tNum+1,:) = [eNum+1,eNum+8,eNum+5];
        end
        
        if es(tri(2),1) == es(tri(3),1)
            newTris(tNum+2,:) = [eNum+2,eNum+4,eNum+7];
        elseif es(tri(2),1) == es(tri(3),2)
            newTris(tNum+2,:) = [eNum+2,eNum+5,eNum+7];
        elseif es(tri(2),2) == es(tri(3),1)
            newTris(tNum+2,:) = [eNum+3,eNum+4,eNum+7];
        elseif es(tri(2),2) == es(tri(3),2)
            newTris(tNum+2,:) = [eNum+3,eNum+5,eNum+7];
        end
        
        newTris(tNum+3,:) = [eNum+6,eNum+7,eNum+8];
        
        tNum = tNum + 4;
        eNum = eNum + 9;
    end
    
    % Add new set of edges to list of edgeTables.
    edgeTables = {edgeTables{:},esNew};
    
    % Update triangles.
    tris = newTris;
end

% Remove duplicate edges.
for i = 1:size(edgeTables,2)
    edgeTables{i} = unique(edgeTables{i},'rows');
end

return;