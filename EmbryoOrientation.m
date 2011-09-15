%%                             EmbryoOrientation.m
% Jacques Bothma                             Last Modified: 09/13/11     
% Levine Lab, UC Berkeley                       
% Functionally complete                          
%
%
%% Attribution:
%  Feel free to use, modify and distribute this code provided that you
%  attribute Jacques Bothma for development.
%  This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
%  To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.
%%
% Overview:
%
% This code orients confocal stacks of an embryo such that the anterior part
% of the embryo (the head) is to the left and the ventral part of the
% embryo (the stomach) is to the bottom if possible. It is based partly on
% shape and partly on where the pole cells are, if it can find them.
%
% Input: I - bw image showing embryo region in white and background in
%            black.
%        Im - Image stack to be flipped.
%        NucBright -Layer of brightest nuclei to try and find pole cells.
%        ManualOrientation - logical one if user want to manually define
%        orientation. Simpy click the quadrant that should be the bottom
%        posterior.
%
% Output:
%
% EmbCorOr - Image file that has been flipped to have the correct
% orientation.
%
% Comments:
%%Inputs: Bw image
%

function EmbCorOr=EmbryoOrientation(I,Im,NucBright,ManualOrientation)

troubleshooting =0 ; 


if ManualOrientation
[yl,xl]=size(I)

imshowbig(NucBright);
hold on
[x,y] = ginput(1)
hold off
fracx=x/xl
fracy=y/yl



if fracx<0.5                          %If RR is less than LL flip image from left to right.  
   EmbCorOr=flipdim(Im,2);
else
   EmbCorOr=Im;                   % Catch to define EmbCorOr if there is no change in orientation
end
        
if fracy<0.5
   EmbCorOr=flipdim(EmbCorOr,1);  %If UU is less than DD flip image from up to down.  
end
    
    
    
else

%%%% Automatically find Pole cells %%%%%%

H = -fspecial('log',50,5); % Filter Kernel
FilterNucBright = imfilter(single(NucBright).*single(I),H,'replicate'); %Apply Filter
FilterNucBright=uint16(FilterNucBright./max(FilterNucBright(:))*(2^16-1)); %Change number format to uint16
BWNucBright = im2bw(FilterNucBright,graythresh(FilterNucBright(:))); % Threshold

imshowbig(BWNucBright)

Lab=bwlabeln(BWNucBright); % Label nucs
STATS=regionprops(Lab,'Area','Centroid');
CentArea = [cat(1,STATS.Centroid),cat(1,STATS.Area)];
SortCentArea=sortrows(CentArea,1);

LeftMost=sort(SortCentArea(1:20,3));
LeftMostSort=sort(LeftMost);

RightMost=sort(SortCentArea(end-10:end,3));
RightMostSort=sort(RightMost);

meanleft = mean(LeftMostSort(end-5:end));
meanright = mean(RightMostSort(end-5:end));

disp(['Mean nuclear size to left is: ' num2str(meanleft) ' right: ' num2str(meanright) ])



C=regionprops(uint8(I),'MajorAxisLength','MinorAxisLength'); %Region properties of black and white image

[w,l]=size(I);                                               %Sizing image

[X,Y]=meshgrid(1:l,1:w);                                     %Creating meshgrid structure to specify x and y coordinates of pixels

Ellipse=(X-l/2).^2/(C.MajorAxisLength(1)/2)^2+(Y-w/2).^2/(C.MinorAxisLength(1)/2)^2<=1;

                                                             % Defining the
                                                             % reference
                                                             % ellipse that
                                                             % has the sam
                                                             % second
                                                             % moments as
                                                             % image

DiffMat=single(I)-single(Ellipse);                           % Take the difference between bw image of
                                                             % embryo and the reference ellipse to
                                                             % look at
                                                             % deviations
                                                             % between the
                                                             % two


LL=sum(DiffMat(X<=l/2));                                     % Sum the differences between the ref ellipse and image for left half of image
RR=sum(DiffMat(X>l/2));                                      % Sum the differences between the ref ellipse and image for right half of image

UU=sum(DiffMat(Y<=w/2));                                     % Sum the differences between the ref ellipse and image for the upper half of image
DD=sum(DiffMat(Y>w/2));                                      % Sum the differences between the ref ellipse and image for the lower half of image

if meanleft>meanright                          %If RR is less than LL flip image from left to right.  
   EmbCorOr=flipdim(Im,2);
else
   EmbCorOr=Im;                   % Catch to define EmbCorOr if there is no change in orientation
end
        
if UU<DD
   EmbCorOr=flipdim(EmbCorOr,1);  %If UU is less than DD flip image from up to down.  
end

   

if troubleshooting

DISP=zeros([w,l,3]);
DISP(:,:,3)=single(I);
DISP(:,:,1)=single(Ellipse);
imshowbig(DISP)

end

end
