function H = buildFeatureHistograms(images, C)
% H = buildFeatureHistograms(images, C)
% Builds a feature histogram for each of a set of images using SIFT
% features.
% Variables:
% H - list of feature histograms.
% images - list of input images.
% C - list of histogram words.
%
% David Pickup

% Initialise list of histograms.
H = cell(size(images));

% Iterate through all images.
for i = 1:numel(images)
    % Extract SIFT features from current image.
    [F,D] = vl_sift(single(images{i}));
    clear F;
    
    % Assign each feature to a word in the codebook.
    I = vl_ikmeanspush(D,C);
    %I = int32(kmeanspush(double(D),double(C)));
    
    % Build a sparse matrix representation for the current image histogram.
    h = spalloc(size(C,2),1,max(size(unique(I))));
    
    % Add words to the histogram.
    if numel(I) > 0
        for j = 1:numel(I)
            h(I(j)) = h(I(j)) + 1;
        end
    end
    
    % Add current histogram to list of histograms.
    H{i} = h;
end

return;