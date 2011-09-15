%%                    NuclearDensityEmbProc.m
% Joe Magliocco & Jacques Bothma                                      
% Levine Lab, UC Berkeley                        
% Functionally complete                             Last Modified: 09/13/10
%
%% Attribution:
% Feel free to use, modify and distribute this code provided that you
% attribute Joe Magliocco & Jacques Bothma for development.
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.
%
% Overview:
%   This code is meant to be a subfunction of EmbProc which takes
%   preprocessed embryo data and calculates the nuclear density.
%
% Inputs:
%   RegionPropNuc - centroids and linear indicies of nuclei
%   ConfReg - region of confidence
%   LabelNucSeg - label matrix for nuclei
%   PixelArea   - Physical Area of pixel in square microns
%   RadErrode - Number of pixels that are erroded from region of confidence
%   
% Outputs:
%   Density - base 2 logarithm of the nuclear density
%
%%

function Density = NuclearDensityEmbProc(LabelNucSeg,ConfReg,PixelArea,RadErrode)

ConfReg = imerode(ConfReg,strel('disk',RadErrode));              % Defines region to use for nuclear density
NucIndConfReg = nonzeros(unique(LabelNucSeg(ConfReg>0)));        % Index of nuclei in the regio of confidence
RegionPropNuc=regionprops(LabelNucSeg,'PixelIdxList');           % Linear indices of each nucleus

ConRegProp=[];

for i=1:length(NucIndConfReg)
    L=mean(ConfReg(RegionPropNuc(NucIndConfReg(i)).PixelIdxList));

    if L==1
    ConRegProp=[NucIndConfReg(i);ConRegProp];
    end
end

S = sum(sum(ismember(LabelNucSeg,ConRegProp))); %determines the area of the region of confidence

count = length(ConRegProp); %returns the number of centroids in the region of confidence

density = count/(S*PixelArea); %returns the density of nuclei by taking the number of nuclei in the region divided by the area of the region

Density = log(density)/log(2); %allows for nuclear densities to cluster around two peaks