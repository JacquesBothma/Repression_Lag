%%                              segmentnucleiWatershed.m
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
% This code segments the nuclei in a confocal projection of an early
% drosophila embryo. The cores of the nuclei is one of the inputs.
% The original image is then slightly blurred and the pixels corresponding
% to the nuclear cores replaced with maximum intesity pixels to guide the
% next watershed segmenting to reduce oversegmentation. The edgepixels from
% the watershed algorithm are then assigned to nuclei based on which direction
% the spatially averaged gradient filter points at that unassigned pixel.
%
% The inputs are the image of the embryo and then filter parameters. The
% output is the lablel matrix of the nuclei. For input
% details see below
%
% Inputs:
%
% I -  fluorescent image of nuclei
% NucCore - bw image of the cores of nuclei (from segmentnucleiCore)
% sizegaus  - size of the gaussian filter to use (10-30)
% radiusgaus - radius of the gaussian filter to use (2-4)
% determining nuclear mask. 
% averagefiltsize - size of square averaging filter used to average image
% gradient used in determing edge pixels.
%
% Outputs:
%
% LabNucNew - Matrix with labelled nuclei
%
% Comments:
% NucCore input determined using segmentnucleiCore
% Compare with segmentnucleiExpand that expands cores and dilates instead of
% using watershed



function LabNucNew=segmentnucleiWatershed(I,NucCore,sizegaus,radiusgaus,averagefiltsize)

troubleshooting=0;

NucCourseFiltT= NucCore;
II=single(I);


%%%%%%%% Fine filtering of image to be used for watershed after core
%%%%%%%% nuclear pixles have been maxed out

H = fspecial('gaussian',sizegaus,radiusgaus);                  % filter kernal for fine gaussian filter
NucFineFilt = imfilter(II,H,'replicate');                      % filter image with fine gaussian filter
NucFineFilt = ((NucFineFilt./max(NucFineFilt(:))));            % scaling

if troubleshooting
imshowbig(NucFineFilt)
end

ModNuc=NucFineFilt.*(~NucCourseFiltT)+NucCourseFiltT;          % redefining core pixel values as maximum intensity
FinalWatershed=watershed(-ModNuc,4);                           % watershed algorithm

if troubleshooting
imshowbig(ModNuc)
end



%%%%% Determining the approximate limits of the embryo %%%%%%%

Backdrop = EmbryoFindPostProc(I,0.5,30,100,40);


%%%%% clean up the net of boundary pixels

[BorderPixels,NUM] = bwlabeln(Backdrop.*~(FinalWatershed>0)); % label the edge pixels
BorderPixels=(BorderPixels==mode(nonzeros(BorderPixels(:))));  % select border pixels that form the largest net
[LabNuc,NUM] = bwlabeln(Backdrop-BorderPixels,4);          %Code that labels unique nuclear segments 

%%%%%% remove pixles that are on the edge of the image so as not to run
%%%%%% into boundary issues

[a,b]=size(I);
E1=[1:a];            
E2=[a*b-a+1:a*b];
E3=[1+a:a:a*b];
E4=[2*a:a:a*b];

EDGE=[E1,E2,E3,E4]; %linear indices of pixels on frame edge

A=find(BorderPixels==1);
Borderpixelsindex=setdiff(A,EDGE); 

%%%%%% Assigning the edge pixels to nuclei based on direction of smoothed image
%%%%%% gradient at the edge pixel.

N=length(Borderpixelsindex);  % Numberof pixels that need to be assigned to nuclei
ImYSize=size(I,1); % Y length of image, used in linear indices.
NeighLinearInd=[ImYSize,ImYSize-1,-1,-ImYSize-1,-ImYSize,-ImYSize+1,1,ImYSize+1]; %linear indices of the 8th closest neighbours of the 0th index.
LabNucNew=LabNuc; % New label matrix for nuclei to use to assign unassigned pixels

[FX,FY] = gradient(-NucFineFilt);  %Gradient of the filtered image

Gfield=FX+sqrt(-1)*FY;             %Combining the X and Y gradietn inforation ito a single number by representing as complex number
H=fspecial('average',averagefiltsize); %Spatial average kernal 
GfieldFilt = imfilter(Gfield,H,'replicate'); %Spatially averaged gradient, the use of complex numbers allows x and y to be automatically seperated

%%%%%% Defining unassigned pixels done by looking at what the nuclear label
%%%%%% is of the nearest pixel that the gradient points to. Done quickly by
%%%%%% just looking at the angle of the gradient vector.

Thet=angle(GfieldFilt)+2*pi;  %Angle of the complex number which corresponds to direction of grad
Thet=mod(Thet,2*pi);      %Modulo 2pi
Thet=Thet/(2*pi)*8;       % 8 nearest neightbours, only relevant which group angle falls into
RefA=mod(round(Thet),8)+1;% Mudulo 8


% Run over all unassigned nuclei and assign (might still be room for
% improvement on method but ok for now)

for i=1:N
     sindex=Borderpixelsindex(i);   %linear index of unassigned pixel
     idx= NeighLinearInd(RefA(sindex))+sindex; %linear index of pixel that grad points to
     LabNucNew(sindex)=LabNuc(idx);    %defining the the value of the pixel to be equal to pixel that the grad points to
      
      
%      if LabNucNew(sindex)==0               %%%%%% If the pixel that gets pointed to has not been assigned then the pixel it points to isused for the assignment
%      idx2=NeighLinearInd(RefA(idx))+idx;
%      LabNucNew(sindex)=LabNuc(idx2);
%      end
      
      if LabNucNew(sindex)==0 %If pixel still not assigned assignment goes anticlockwise along neightbouring pixels untill pixel is assgned.
      ct=1;
        while LabNucNew(sindex)==0 && ct < 9;
        idx= NeighLinearInd(mod(RefA(sindex)-(1+ct),8)+1)+sindex;
        LabNucNew(sindex)=LabNuc(idx);
        ct=ct+1;
        end
      end
   
 end


%plotting 
CM=label2rgb(LabNucNew, 'jet', [0,0,0],'shuffle');
DI = cast(bsxfun(@times,single(CM)/255,single(I)),class(I));
imshowbig(DI)
   





