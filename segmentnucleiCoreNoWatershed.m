
%%                              segmentnucleiCoreNoWatershed.m
% Jacques Bothma                                      
% Levine Lab, UC Berkeley                        
% Functionally complete                             Last Modified: 09/13/10
%
%% Attribution:
% Feel free to use, modify and distribute this code provided that you
% attribute Jacques Bothma for development.
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.
%
%
% Overview:
% This code segments the nuclei in a confocal projection of an early
% drosophila embryo to find the cores of the nuclei.
% Output is a binary image of the cores of the nuclei.
%
% Input:
%
% I -  fluorescent image of nuclei
% sizefilt -  size of the initial filter to roughly find nuclei (10- 30)
% radiusfilt -  radius of the initial filter to roughly find nuclei  (5-10)
% coursefilterthreshold - threshold for segmenting course filter (if choose
% 0 threshold will be chosen automatically using Otsu's method) 
%
% Output:
%
% NucCourseFiltT - bw image of core nuclei
%
% Comments:
% The output is used as input for other functions that do a more complete
% segmentation of the nuclei. It differs from segmentnucleiCore in that it
% doesn't use a watershed step at the end to speed it up.


function NucCourseFiltT=segmentnucleiCoreNoWatershed(I,sizefilt,radiusfilt,coursefilterthreshold)

troubleshooting=0;

%%%%%% Finding core nuclear pixels %%%%%%%%%%%

H = -fspecial('log',sizefilt,radiusfilt);                      % Defining the first filter kernel
II=single(I);                                                  % Using single format speeds things up
NucCourseFilt = imfilter(II,H,'replicate');                    % Apply first coarse filter

NucCourseFilt = (NucCourseFilt./max(NucCourseFilt(:)));        % Scale filtered result
NucCourseFilt=NucCourseFilt.*(NucCourseFilt>0);                % Remove negative pixles.

if troubleshooting
imshowbig(NucCourseFilt)
end

if coursefilterthreshold==0
NucCourseFiltT=im2bw(NucCourseFilt,graythresh(NucCourseFilt)); % Threshold the filtered nuclear image with automatic thresholding
else
NucCourseFiltT=im2bw(NucCourseFilt,coursefilterthreshold);     % Threshold the filtered nuclear image user chosen threshold
end


NucCourseFiltT=logical(NucCourseFiltT);     % Segmented core pixels

%%%%%%%%