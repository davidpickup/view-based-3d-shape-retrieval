function features = extractImageFeatures(images)
% features extractImageFeatures(images)
% Extracts a set of SIFT features from a set of images.
% Variables:
% features - extracted SIFT features.
% images - list of input images.
%
% David Pickup 2013



% Initialise list of SIFT features.
features = [];

% For printing.
count = 1;

% Iterate through all images.
for i = 1:numel(images)
    % Compute SIFT featurns from current image.
    [F,D] = vl_sift(single(images{i}));
    clear F;
    
    % Add current SIFT features to list of all features.
    features = [features, D];
    
    % For printing.
    if mod(i,66) == 0
        %['SIFT for model ' int2str(count) ' computed.']
        count = count + 1;
    end
end

return;