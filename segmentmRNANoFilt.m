%%                              segmentmRNA.m
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
%
% Overview:
% This code performs the segmentation of raw mRNA. It does a laplacian of
% gaussian filter to get rid of noise and then follows that with a watershed
% segmentation to split the fused mRNAs. 
% 
% 
% Inputs & outputs:
% The input is a raw image of the mRNA and some parameters. The output
% is the labell matrix of mRNA and a logical matrix that has ones where the intensity weighted centroids 
% are of each mRNA molecule.
%
% I - Raw mRNA image
% sizefilt -  size of the initial filter to roughly filter noise  (10-30)
% radiusfilt -  radius of the initial filter to roughly find mRNA (1-3)
% threshold - image threshold to define difference between background and
% noise, putting zero forces an automatic threshold to be determined using
% the mRNAthreshold function.
% N - Number of points in threshold range to use  (~100)
% smo- number of points to use in smoothing data (~4)
% Thresmin - Minimum threshold value to use (~0.01)
% Thresmax - Maximum threshold value to use (~1)
%
%
% Custom Functions Used:
% mRNAthreshold(I,N,smo,Thresmin,Thresmax)
%


function [LabmRNANew,mRNAICenter,threshold] = segmentmRNANoFilt(I,threshold,N,smo,ThresMin,ThresMax)

troubleshooting =1;


II=single(I);                                                  % Using single format speeds things up
II=II./max(II(:));



if troubleshooting
    imshowbig(label2rgb(round((II-min(II(:)))/(max(II(:))-min(II(:)))*10000),'jet',[0,0,0]))
    title('Heat map of image') 
end

if threshold==0
     if (mean(I(I>mean(nonzeros(I(:))) + 2*std(double(nonzeros(I(:))))))/mean(nonzeros(I(:))))<=2
     threshold=1;
     else
     threshold=mRNAthreshold(II,N,smo,ThresMin,ThresMax);
     
     figure
     imshowbig(II)
     
     end
mRNAFiltT=im2bw(II,threshold); % Threshold the filtered nuclear image with automatic thresholding
else
mRNAFiltT=im2bw(II,threshold);                      % Threshold the filtered nuclear image user chosen threshold
end

mRNAWatershed=watershed(-mRNAFiltT);                            % Watershed the filtered image
mRNAWatershedT=logical(mRNAWatershed).*logical(mRNAFiltT);     % Segmented core pixels

[LabmRNANew,NumbermRNA]=bwlabeln(mRNAWatershedT);                           % label mRNA objects



s = regionprops(LabmRNANew, {'PixelIdxList'}); %Determine the indices of the pixels that correspond to each mRNA

intensitycentroid=zeros([NumbermRNA,2]); %prepare 
 
for i=1:NumbermRNA,[r,c]=ind2sub(size(I),s(i).PixelIdxList);
     intensitycentroid(i,1)=round(sum(r.*double(I(s(i).PixelIdxList)))/sum(double(I(s(i).PixelIdxList))));
     intensitycentroid(i,2)=round(sum(c.*double(I(s(i).PixelIdxList)))/sum(double(I(s(i).PixelIdxList))));
end

mRNAICenter=zeros(size(I)); mRNAICenter(sub2ind(size(I),intensitycentroid(:,1),intensitycentroid(:,2)))=1;


