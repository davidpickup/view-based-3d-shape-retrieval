function images = depthImages(mesh, views, resolution)
% images = depthImages(mesh, views, resolution)
% Renders a depth image of a mesh from each of a set of views.
% Variables:
% images - set of rendered depth images.
% mesh - mesh structure.
% views - set of camera viewpoints.
% resolution - resolution of images.
%
% David Pickup 2013

% Initialise set of images.
images = cell(size(views,1),1);

% Get mesh data.
verts = [mesh.X,mesh.Y,mesh.Z]';
nTris = size(mesh.TRIV,1);
mesh.TRIV = mesh.TRIV';
nViews = size(views,1);

% Is the script running in OpenGL Psychtoolbox? Abort, if not.
AssertOpenGL;

% Find the screen to use for display:
screenid=max(Screen('Screens'));

% Initialise OpenGL.
InitializeMatlabOpenGL;

% Open a window.
win = Screen('OpenWindow',screenid,0,[100,100,100+resolution(1),100+resolution(2)]);

% Iterate through the number of views.
for i = 1:nViews
    % Setup OpenGL.
    Screen('BeginOpenGL',win);
    
    glClearColor (1.0, 1.0, 1.0, 0.0);
	glClearDepth(1.0);
    glEnable(GL.DEPTH_TEST);

    % Setup projection matrix.
    glMatrixMode(GL.PROJECTION);
    glLoadIdentity;
    glOrtho(-1, 1, -1, 1, 0.0, 2.0);

    % Position light.
    glMatrixMode(GL.MODELVIEW);
    glLoadIdentity;
    
    if (((i == 3 || i == 4) && nViews > 3) || (nViews == 3 && i == 2))
            % Setup the current camera viewpoint.
            gluLookAt(views(i,1), views(i,2), views(i,3), 0, 0, 0, 0, 0, 1);
    else
            % Setup the current camera viewpoint.
            gluLookAt(views(i,1), views(i,2), views(i,3), 0, 0, 0, 0, 1, 0);
    end

    %glPolygonMode(GL.FRONT_AND_BACK, GL.FILL);


    % Clear OpenGL buffers.
    glClear;

    glEnableClientState(GL.VERTEX_ARRAY);
    glVertexPointer(3, GL.DOUBLE, 0, verts);

    glDrawElements(GL.TRIANGLES, nTris*3, GL.UNSIGNED_INT, uint32(mesh.TRIV-1));

    Screen('EndOpenGL', win);
    Screen('Flip', win);
    
    images{i} = glReadPixels(0,0,resolution(1),resolution(2),GL.DEPTH_COMPONENT,GL.FLOAT);
end

% End OpenGL rendering and close window.
Screen('CloseAll');

return;