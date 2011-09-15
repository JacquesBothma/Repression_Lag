%%                             EveningFilter.m
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
%  Important Notes:
% % This version written for Windows.
%
% Overview:
%
% This code evens out the projected nuclear intensity. It performs a
% gaussian blur and then divides original image by the blurred image. Some
% additional steps are included to ensure that the low intensity regions of
% the image are not overemphasized.
%
% Input:
%
% Nuc -  fluorescent image of nuclei
% SizGaus -  size of gaussian filter to smooth nuclear image (~40)
% RadGaus -  radius of gaussian filter to smooth nuclear image (~20)
% varagin -  optional to define region of interest before performing the
% evening filter. Important for inhomogeneous images.
%
% Output:
%
% Nuc - Image where evening filter has been applied
%
% Comments:
%

function NucE=EveningFilter(Nuc,SizGaus,RadGaus,varagin)

troubleshooting=0;

Imageclass=class(Nuc); % Determine class of input image to be used in output 

H=fspecial('gaussian',SizGaus,RadGaus);        % Gaussian filter kernel 
FilteredImage = imfilter(single(Nuc),H,'replicate');   % Filtered image 

if troubleshooting % heat map of log filtered images
    imshowbig(label2rgb(round((FilteredImage-min(FilteredImage(:)))/(max(FilteredImage(:))-min(FilteredImage(:)))*10000),'jet',[0,0,0]))
    title('Heat map of filtered image')
end


if nargin == 3
BWfiltT=EmbryoFindPostProc(Nuc,0.5,30,100,40); % Thresholded filtered image, to suppress emphasizing low intensity regions
else
BWfiltT=varagin;
end


BWfiltTfilt = imfilter(single(BWfiltT),H,'replicate');  % Smooth edges of mask

N=single(Nuc)./(FilteredImage).*BWfiltTfilt;% Divide original image by smoothed image

Nn=N./max(N(:)); % normalize 


if strcmp(Imageclass,'uint16') % Cast image in appropriate format
    NucE=im2uint16(Nn);
elseif strcmp(Imageclass,'uint8')
    NucE=im2uint8(Nn);
else
    NucE=double(Nn);
end

if troubleshooting
    imshowbig(NucE)
end
