function D = computeDistancesForDataset(location, fname, nModels, startStage)
% D = computeDistancesForDataset(location, fname, nModels, startStage)
% Computes a distance matrix for a given shape dataset.
% Variables:
% D - output distance matrix.
% location - directory location of dataset.
% fname - filename of models (must be .obj files).
% nModels - number of models in the dataset.
% startStage -  1 = compute everything.
%               2 = codebook already calculated.
%               3 = features already calculated.
%               4 = no clock matching (single measure per model).
%
% David Pickup 2013

%% Render depth images for each model.
timer = tic;

% Compute viewpoints for rendering.
[et vertexTable vs tris] = geodesicSphere(2);
clear et tris;

if startStage < 3

    % Initialise grouped set of all images.
    groupedImages = cell(nModels,1);

    % Initialised ungrouped set of all images.
    images = {};
    
    w = waitbar(0,['Computing depth images. 0/' int2str(nModels) '. ETA: inf:inf:inf']);
    timer = tic;
    % Iterate through all models.
    for i = 1:nModels
        mesh = readMesh(fullfile([location fname int2str(i-1) '.obj']));
        mesh = normaliseMesh(mesh);
        mesh = poseNormalisationPCAplusRECT(mesh);
        currImages = depthImages(mesh, vs, [256 256]);
        
        % Add its images to total sets.
        groupedImages{i} = currImages;
        images = {images{:} currImages{:}};
        seconds = toc(timer);
        seconds = (seconds/i) * (nModels-i);
        [hours,minutes,seconds] = secondsToTime(seconds);
        waitbar(i/nModels,w,['Computing depth images. ' int2str(i) '/' int2str(nModels) '. ETA: ' int2str(hours) ':' int2str(minutes) ':' int2str(seconds)]);
    end
    close(w);
    clear w
    clear mesh
    clear meshOrig
    clear conFactor
    clear vs
end

%% If needed, compute feature codebook, otherwise load it from file.
if startStage < 2
    features = extractImageFeatures(images);
    C = int32(buildCodebookRandom(features, 1500));
    clear images
    save([location 'C.mat'], 'C');
    'Created and saved codebook.'
else
    load([location 'C.mat']);
    'Loaded codebook from file.'
end
clear images

%% Compute a feature histogram for each model and save.

if startStage < 3
    % Initialise list of feature histograms.
    Hists = cell(nModels,1);

    % Iterate through all models.
    for i = 1:nModels
        % Compute the feature histogram.
        H = buildFeatureHistograms(groupedImages{i}, C);
        
        % Add the histogram to total list.
        Hists{i} = H;
        %['Created histogram for model ' int2str(i)]
    end
    save([location 'Hists.mat'], 'Hists');
    clear groupedImages
    'Created feature histograms.'
else
    load([location 'Hists.mat']);
    'Loaded histograms from file.'
end

descTime = toc(timer); % Time to compute descriptors.
clear groupedImages

%% Produce distance matrix between each pair of models.
timer = tic; % time to compare models.

% Initialise distance matrix.
D = zeros(nModels);

% Load clockmatching permutations from file.
load permutations;

% Clockmatching.
D = clockmatchModels(Hists, permutations, vertexTable, pEdgeTables);

compareTime = toc(timer); % time to compare models.
save([location 'timings.mat'], 'descTime', 'compareTime');

% Save distance matrix to file.
save([location 'D.mat'], 'D');
'Done.'

return;