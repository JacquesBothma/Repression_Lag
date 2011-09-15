%%                             EmbryoFindNew.m
% Jacques Bothma                                    Last Modified: 09/13/11     
% Levine Lab, UC Berkeley                       
% Functionally complete                          
%
%
%% Attribution:
%  Feel free to use, modify and distribute this code provided that you
%  attribute Jacques Bothma for development.
%  This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
%  To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.

% Overview:
%
% This code finds the core region that corresponds to where embryos are in
% raw confocal projections stained for nuclei. It rescales the image size to speed up
% processing time. It uses a three pronged approach to be able to seperate structures on different length scales.
% the first step is to filter the image with a broad LOG to find edges of the embryo. A Qhull algorithm is then used to remove concave
% irregularities in the object or reveal missegmentation. If the embryo is missegmented which is defined by occupying too large or small a fraction of the 
% image area or touching the border in too many places the filter step is
% repeated with a smaller LOG filter and then dilated. If this again
% results in missegmentation an even smaller LOG filter is applied.
%
% Input:
%
% I -  fluorescent image of nuclei
% scalingfactor -  scaling factor to reduce size of image by to speed up
% processing (~0.5 for a 2048x2048 image)
%
% logsize - size of log filter(~150)
% logradius - radius of log filter(~30)
% closeradius - radius of imclose filter(~30)
% factor1 - factor by which to reduce first log filter radius if it doesn't
% work (~0.5)
% factor2 - factor by which to log filter radius if it doesn't
% work a second time (~0.4)
% dil1 - intelligent opening radius for second filter (~6)
% dil2 - intelligent opening radius for second filter (~6)
% maxperm - maximum fraction of the perimter that the embryo can occupy
% before it is considered to be missegmented (1/8)
% minarea - minimum fraction of the total area that the embryo can occupy
% before it is considered to be missegmented (0.1)
% maxarea - maximum fraction of the total area that the embryo can occupy
% before it is considered to be missegmented. (0.5)
%
% Output:
%
% EmbryoRegion - bw image that has embryo regions as ones
%
% Comments:

function EmbryoRegion = EmbryoFindNew(I,scalingfactor,logsize,logradius,closeradius,factor1,factor2,dil1,dil2,maxperm,minarea,maxarea)


warning('off','Images:initSize:adjustingMag')         % Suppress warning about magnification

troubleshooting=0;                                    % Suppreses display when not troubleshooting

%%%%% Sizing and resclaing image to speed up calculations %%%%%%%%%%%%%

[ro,co]=size(I);                                      % Size of original image
r = round(scalingfactor*ro);                          % Size to rescale image to
c = round(scalingfactor*co);
IE=imresize(I,[r,c]);                                 % Scale down size of image to speed up calculations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


framepixels=[1:r,c*r:-1:c*(r-1),2*r:c:r*(c-1),r+1:r:r*(c-2)+1]; % linear indices of pixels on boundary in rescaled image


if troubleshooting
    imshowbig(IE)                                                       %Display
end

H=fspecial('log',ceil(logsize*scalingfactor),ceil(logradius*scalingfactor));  %Laplacian of gaussian filter kernal to find edges, large scale
NucLogFine = imfilter(single(IE),H,'replicate');                              %Filtering of image


if troubleshooting % heat map of log filtered images
    imshowbig(label2rgb(round((NucLogFine-min(NucLogFine(:)))/(max(NucLogFine(:))-min(NucLogFine(:)))*10000),'jet',[0,0,0]))
end

BW = im2bw(NucLogFine,graythresh(NucLogFine(NucLogFine>0)));       %Thresholding of filtered image to show boundaries

if troubleshooting
    imshowbig(BW)                                                  % Display
end

BW=imclose(BW,strel('disk',ceil(scalingfactor*closeradius)));      % Closing image to try and join boundary elements seperated by gaps

if troubleshooting
    imshowbig(BW)                                                  % Display
end

BWL=bwlabeln(BW);                                                  % Label of boundary elements

LabBWLBig=BWL==mode(double(nonzeros(BWL)));                        % Selecting largest boundary element

LabBWLBig2=ismember(BWL,nonzeros(unique(BWL(framepixels))));       % Selecting all elements on image boundary

LabBWLBig=LabBWLBig|LabBWLBig2;                                    % Combination of largest boundary element and boundary elements on periphery

if troubleshooting
    imshowbig(LabBWLBig)                                           % Display
end

BW=(~LabBWLBig);                                                   % Inverse of boundary to find solid objects

if troubleshooting
    imshowbig(BW)                                                  % Display
end

LabN=bwlabeln(BW,4);                                               % Labelling of solid objects 4 not 8 to find contained


%%%%%%%%%%%%%%%%%% Finding biggest and brightest object %%%%%%%%%%%%

SizeAndIntensitySort=zeros(max(LabN(:)),3);    % Initialize matrix to store

for i=1:max(LabN(:));                          % Determine size and total intensity of all objects

    SizeAndIntensitySort(i,1)=i;                                    % Index
    SizeAndIntensitySort(i,2)=sum(IE(find(ismember(LabN,i))));      % Total intensity
    SizeAndIntensitySort(i,3)=length(find(ismember(LabN,i)));       % Size

end

SizeAndIntensitySort=sortrows(SizeAndIntensitySort,-2);   % Sort ascendingly according to total intensity

[dummy,Bright] = max(SizeAndIntensitySort(:,2));           % Find index that corresponds to object with highes intensity

LabNBig =LabN==SizeAndIntensitySort(Bright,1);            % bw of brightest object

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if troubleshooting
    imshowbig(LabNBig)
end

Boundaries =  cell2mat(bwboundaries(LabNBig));                    %Boundaries of solid object

ConInd = convhull(Boundaries(:,1),Boundaries(:,2));               %Pixels that correspond to convex hull
XX = Boundaries(ConInd,1);                                        %Boundaries
YY = Boundaries(ConInd,2);
LabNBig = poly2mask(YY,XX,r,c);                                   %Boundary pixels to mask


if troubleshooting
    imshowbig(LabNBig)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LabNBigDil=imdilate(LabNBig,strel('disk',2));                     % Dilate region to make sure region touches edges if it's close
RatioPerm = sum(LabNBigDil(framepixels))./length(framepixels);    % Fraction of image edge made up of region
RatioArea = sum(LabNBig(:))/(r*c);                                % Fraction of image area made up of region


if RatioPerm>maxperm || (RatioArea<minarea) || (RatioArea>maxarea)

    %%%%%%%%%%%%%% If area still too big or small or too much of periphery made
    %%%%%%%%%%%%%% up of region third pass. This is meant to sepearte only
    %%%%%%%%%%%%%% the images where embryos are touching very closely.
    %%%%%%%%%%%%%% Need to resolve fine structure to do this effectively
    %%%%%%%%%%%%%% which is why the strategies and parametgers are not
    %%%%%%%%%%%%%% effective for easier to seperate embryos.



    H=fspecial('log',ceil(logsize*scalingfactor),ceil(logradius*scalingfactor*factor1));  % Laplacian of gaussian filter kernal to find edges, large scale
    NucLogFine = imfilter(single(IE),H,'replicate');                                  % Filtering of image
    BW = im2bw(NucLogFine,graythresh(NucLogFine(NucLogFine>0)));                      %
    BWO=BW;


    if troubleshooting % heat map of log filtered images
        imshowbig(label2rgb(round((NucLogFine-min(NucLogFine(:)))/(max(NucLogFine(:))-min(NucLogFine(:)))*10000),'jet',[0,0,0]))
    end


    if troubleshooting
        imshowbig(BW)                                                      %Display
    end

    BW = bwdist(BW)<dil1;

    BW=bwmorph(BW,'thin',dil1)|BWO;

    if troubleshooting
        imshowbig(BW)                                                      %Display
    end

    if troubleshooting
        imshowbig(BW)                                                      %Display
    end

    BWL=bwlabeln(BW);                                                  %Label black and white image


    LabBWLBig=BWL==mode(double(nonzeros(BWL)));                        %Selecting largest boundary element

    LabBWLBig2=ismember(BWL,nonzeros(unique(BWL(framepixels))));       %Selecting all elements on image boundary

    LabBWLBig=LabBWLBig|LabBWLBig2;

    %
    if troubleshooting
        imshowbig(LabBWLBig)                                               %Display
    end



    BW=(~LabBWLBig);                                                   %Inverse of boundary to find solid objects

    if troubleshooting
        imshowbig(BW)                                                      %Display
    end

    LabN=bwlabeln(BW,4);                                                 %Labelling of solid objects

 %%%%%%%%%%%%%%%%%% Finding biggest and brightest object %%%%%%%%%%%%

SizeAndIntensitySort=zeros(max(LabN(:)),3);    % Initialize matrix to store

for i=1:max(LabN(:));                          % Determine size and total intensity of all objects

    SizeAndIntensitySort(i,1)=i;                                    % Index
    SizeAndIntensitySort(i,2)=sum(IE(find(ismember(LabN,i))));      % Total intensity
    SizeAndIntensitySort(i,3)=length(find(ismember(LabN,i)));       % Size

end

SizeAndIntensitySort=sortrows(SizeAndIntensitySort,-2);   % Sort ascendingly according to total intensity

[dummy,Bright] = max(SizeAndIntensitySort(:,2));           % Find index that corresponds to object with highes intensity

LabNBig =LabN==SizeAndIntensitySort(Bright,1);            % bw of brightest object

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Boundaries =  cell2mat(bwboundaries(LabNBig));                     %Boundaries of solid object

    ConInd = convhull(Boundaries(:,1),Boundaries(:,2));                %Pixels that correspond to convex hull
    XX = Boundaries(ConInd,1);                                         %Boundaries
    YY = Boundaries(ConInd,2);
    LabNBig = poly2mask(YY,XX,r,c);



    if troubleshooting
        imshowbig(LabNBig)
    end

    
    LabNBigDil=imdilate(LabNBig,strel('disk',2));                     % Dilate region to make sure region touches edges if it's close
    RatioPerm = sum(LabNBigDil(framepixels))./length(framepixels);    % Fraction of image edge made up of region
    RatioArea = sum(LabNBig(:))/(r*c);            

    if RatioPerm>maxperm || (RatioArea<minarea) || (RatioArea>maxarea)

        %%%%%%%%%%%%%% If area still too big or small or too much of periphery made
        %%%%%%%%%%%%%% up of region third pass. This is meant to sepearte only
        %%%%%%%%%%%%%% the images where embryos are touching very closely.
        %%%%%%%%%%%%%% Need to resolve fine structure to do this effectively
        %%%%%%%%%%%%%% which is why the strategies and parametgers are not
        %%%%%%%%%%%%%% effective for easier to seperate embryos.





        H=fspecial('log',ceil(logsize*scalingfactor),ceil(logradius*scalingfactor*factor2));  % Laplacian of gaussian filter kernal to find edges, large scale
        NucLogFine = imfilter(single(IE),H,'replicate');                                  % Filtering of image
        BW = im2bw(NucLogFine,graythresh(NucLogFine(NucLogFine>0)));                      %
        BWO=BW;


        if troubleshooting % heat map of log filtered images
            imshowbig(label2rgb(round((NucLogFine-min(NucLogFine(:)))/(max(NucLogFine(:))-min(NucLogFine(:)))*10000),'jet',[0,0,0]))
        end


        if troubleshooting
            imshowbig(BW)                                                      %Display
        end


        BW = bwdist(BW)<dil2;

        BW=bwmorph(BW,'thin',dil2)|BWO;

        if troubleshooting
            imshowbig(BW)                                                      %Display
        end

        %        REG=BW;

        if troubleshooting
            imshowbig(BW)                                                      %Display
        end

        BWL=bwlabeln(BW);                                                  %Label black and white image


        LabBWLBig=BWL==mode(double(nonzeros(BWL)));                        %Selecting largest boundary element

        LabBWLBig2=ismember(BWL,nonzeros(unique(BWL(framepixels))));       %Selecting all elements on image boundary

        LabBWLBig=LabBWLBig|LabBWLBig2;

        %
        if troubleshooting
            imshowbig(LabBWLBig)                                               %Display
        end



        BW=(~LabBWLBig);                                                   %Inverse of boundary to find solid objects

        if troubleshooting
            imshowbig(BW)                                                      %Display
        end

        LabN=bwlabeln(BW,4);                                                 %Labelling of solid objects

%%%%%%%%%%%%%%%%%% Finding biggest and brightest object %%%%%%%%%%%%

SizeAndIntensitySort=zeros(max(LabN(:)),3);    % Initialize matrix to store

for i=1:max(LabN(:));                          % Determine size and total intensity of all objects

    SizeAndIntensitySort(i,1)=i;                                    % Index
    SizeAndIntensitySort(i,2)=sum(IE(find(ismember(LabN,i))));      % Total intensity
    SizeAndIntensitySort(i,3)=length(find(ismember(LabN,i)));       % Size

end

SizeAndIntensitySort=sortrows(SizeAndIntensitySort,-2);   % Sort ascendingly according to total intensity

[dummy,Bright] = max(SizeAndIntensitySort(:,2));           % Find index that corresponds to object with highes intensity

LabNBig =LabN==SizeAndIntensitySort(Bright,1);            % bw of brightest object

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Boundaries =  cell2mat(bwboundaries(LabNBig));                     %Boundaries of solid object

        ConInd = convhull(Boundaries(:,1),Boundaries(:,2));                %Pixels that correspond to convex hull
        XX = Boundaries(ConInd,1);                                         %Boundaries
        YY = Boundaries(ConInd,2);
        LabNBig = poly2mask(YY,XX,r,c);



        if troubleshooting
            imshowbig(LabNBig)
        end

    end


end


LabNBig=imresize(LabNBig,[ro,co]);                        %Resize image to original size

EmbryoRegion=LabNBig;

