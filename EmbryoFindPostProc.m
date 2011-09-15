%%                             EmbryoFindPostProc.m
% Jacques Bothma                                      
% Levine Lab, UC Berkeley                        
% Functionally complete                             Last Modified: 09/13/11
%
%
%% Attribution:
% Feel free to use, modify and distribute this code provided that you
% attribute Jacques Bothma for development.
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.
%
%
%
%  Important Notes:
% % This version written for Windows.
%
% Overview:
%
% This code finds the region of an image that contains an embryo by using the channel stained for nuclei. 
% It rescales the image size to speed up processing time and then opens the image to blur small features. It then filters with
% a LOG (Laplacian of Gaussian) to find edges. All edges touching the boundary are then used to define the object
% of interest, can be done with processed images. A Qhull algorithm is then used to remove concave
% irregularities in the object before it is rescaled back to original size.
%
% Input:
%
% I -  fluorescent image of nuclei
% scalingfactor -  scaling factor to reduce size of image to speed up
% processing (~0.5 for a 2048x2048 image)
%
% openingradius -  radius of opening filter (~30 for 2048)
% logsize - size of log filter(~100)
% logradius - radius of log filter(~40)
%
%
% Output:
%
% EmbryoRegion - bw image that shows region of embryo in white
%
% Comments:

function EmbryoRegion = EmbryoFindPostProc(I,scalingfactor,openingradius,logsize,logradius)

troubleshooting=0; % switch on display for trouble shooting

[ro,co]=size(I);

r = round(scalingfactor*ro);
c = round(scalingfactor*co);

IE=imresize(I,[r,c]);                                               %Reduce size of image to speed up calculations

Io = imopen(IE, strel('disk', ceil(openingradius*scalingfactor)));  %Opening image, effectively blurring it

if troubleshooting==1
imshowbig(Io)                                                       %Display
end

H=fspecial('log',ceil(logsize*scalingfactor),ceil(logradius*scalingfactor));  %Laplacian of gaussian filter kernal to find edges, large scale
NucLogFine = imfilter(single(Io),H,'replicate');                   %Filtering of image

if troubleshooting==1
figure, surf(double(NucLogFine)), shading interp
end

BW = im2bw(NucLogFine,graythresh(NucLogFine(NucLogFine>0)));       %Thresholding of filtered image to show edges 

if troubleshooting==1
imshowbig(BW)                                                      %Display
end

BW=imclose(BW,strel('disk',ceil(20*scalingfactor)));
if troubleshooting==1
imshowbig(BW)                                                      %Display
end


BWL=bwlabeln(BW);                                                  %Label black and white image


[l,w]=size(BWL);

framepixels=[1:l,w*l:-1:w*(l-1),2*l:l:l*(w-1),l+1:l:l*(w-2)+1];

LabBWLBig=ismember(BWL,(double(nonzeros(BWL(framepixels)))));      %Selecting all boundary elements on periphery

if troubleshooting==1
imshowbig(LabBWLBig)                                               %Display
end

BW=(~LabBWLBig);                                                   %Inverse of boundary to find solid objects

if troubleshooting==1
imshowbig(BW)                                                      %Display
end

LabN=bwlabeln(BW);                                                 %Labelling of solid objects

LabNBig=LabN==mode(double(nonzeros(LabN)));                        %Selection of biggest solid object

Boundaries =  cell2mat(bwboundaries(LabNBig));                     %Boundaries of solid object

ConInd = convhull(Boundaries(:,1),Boundaries(:,2));                %Pixels that correspond to convex hull
XX = Boundaries(ConInd,1);                                         %Boundaries 
YY = Boundaries(ConInd,2);


[l,w]=size(LabNBig);

LabNBig = poly2mask(YY,XX,l,w);                                   %Boundary pixels to mask


LabNBig=imresize(LabNBig,[ro,co]);                        %Resize image to original size

EmbryoRegion=LabNBig;

