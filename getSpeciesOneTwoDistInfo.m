function fracDistSpeciesOneHigher = getSpeciesOneTwoDistInfo(speciesOneDistFileName, speciesTwoDistFileName)

% This method computes takes files with the distances of PWMs of two
% species from another species and computes the fraction of the common
% motifs for which the distance from species one is greater than the
% distance from species two

% Input:
%   1.  speciesOneDistFileName: Name of the file with the distances of PWMs
%       in the first species from PWMs in the other species
%   2.  speciesTwoDistFileName: Name of the file with the distances of PWMs
%       in the second species from PWMs in the other species
% Output:
%   1.  fracDistSpeciesOneHigher: Fraction of the common motifs for which 
%       the distance from species one is greater than the distance from 
%       species two

TomtomOutBestDistsSpeciesOne = importdata(speciesOneDistFileName);
TomtomOutBestDistsSpeciesTwo = importdata(speciesTwoDistFileName);
speciesOneIntersectIndexes = [];
speciesTwoIntersectIndexes = [];
for i = 1:length(TomtomOutBestDistsSpeciesTwo.data)
    % Iterate through the motifs in species two and find those that are
    % also in species one
    TF = TomtomOutBestDistsSpeciesTwo.textdata{i, 1};
    motif = TomtomOutBestDistsSpeciesTwo.textdata{i, 2};
    speciesOneTFLocations = find(strcmp(TomtomOutBestDistsSpeciesOne.textdata(:,1), TF) == 1);
    speciesOneMotifLocations = find(strcmp(TomtomOutBestDistsSpeciesOne.textdata(:,2), motif) == 1);
    index = intersect(speciesOneTFLocations, speciesOneMotifLocations);
    if ~isempty(index)
        % The TF, motif pair is in both species, so record the index in
        % each species
        speciesOneIntersectIndexes = vertcat(speciesOneIntersectIndexes, index);
        speciesTwoIntersectIndexes = vertcat(speciesTwoIntersectIndexes, i);
    end
end
scatter(TomtomOutBestDistsSpeciesOne.data(speciesOneIntersectIndexes), TomtomOutBestDistsSpeciesTwo.data(speciesTwoIntersectIndexes), '.')
[c, p] = corr(TomtomOutBestDistsSpeciesOne.data(speciesOneIntersectIndexes), TomtomOutBestDistsSpeciesTwo.data(speciesTwoIntersectIndexes))
speciesOneHigher = TomtomOutBestDistsSpeciesOne.data(speciesOneIntersectIndexes) > TomtomOutBestDistsSpeciesTwo.data(speciesTwoIntersectIndexes);
fracDistSpeciesOneHigher = length(find(speciesOneHigher == 1))/length(speciesOneHigher)