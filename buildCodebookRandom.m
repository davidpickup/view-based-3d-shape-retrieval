function C = buildCodebookRandom(features, k)
% C = buildCodebookRandom(images, k)
% Builds a codebook from a set of features.
% Variables:
% C - output codebook.
% features - set of features.
% k - number of words in the codebook.
%
% David Pickup 2013

% Choose a set of random features for the codebook.
C = (zeros(size(features,1),k));
rs = zeros(k,1);
for i = 1:k
    r = randi(size(features,2));
    while(ismember(r,rs))
        r = randi(size(features,2));
    end
    rs(i) = r;
    C(:,i) = features(:,r);
end
'Created codebook.'
return;