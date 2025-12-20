function SpatialInformation = GetSpaInfo(OCCMat, SkaggsrateMat)
%Calculate spatial information [Skaggs']
%SpatialInformation = GetSpaInfo(OCCMat, SkaggsrateMat)
%
%Related function(s)
%jkABM.m
%outputs from jkABM.m are needed to calculate spatial information scores.
%
%Input format
% OCCMat	[nROW x nCOL matrix]									raw occupancy amp
% SkaggsrateMat [nROW x nCOL matrix]								Skaggs' ABM firing rate map
%
%Output format
% SpatialInformation [1 x 1]										bit/spike spatial information scores
%
%Originally from VB codes which Inah Lee has.
%Translate VB codes to matlab was done by [Jangjin Kim, July-13-2008]
%Verification was done

criteria = .0001;

SumRate = sum(sum(SkaggsrateMat, 'omitnan'), 'omitnan');
SumOcc = sum(sum(~isnan(OCCMat)));
MeanRate = SumRate / SumOcc;

H1 = sum(SkaggsrateMat(SkaggsrateMat > criteria) * -1 .* log2(SkaggsrateMat(SkaggsrateMat > criteria)) / SumOcc);
H0 = MeanRate * -1 * log2(MeanRate);

SpatialInformation = (H0 - H1) / MeanRate;