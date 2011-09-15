
%%                              mRNAthreshold.m
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
% This code determines the threshold to be chosen for segmenting intronic
% mRNA signals from a filtered image of the mRNA. The method is based on
% the idea that when the fraction of objects that consist of more than one
% pixel peak as a function of the threshold  that you have come to 
% a good approximation of the optimal threshold. The
% code finds the number of objects that consist of one pixel as a function
% of the threshold and then finds the maximum of that for a smoothed
% version of the numbers. If there is a problem or what seems like an
% excessive number of objects the code defaults the threshold to 1 which
% excludes all objects.
% 
% 
% Inputs & outputs:
% The input is a filtered image of the mRNA and some parameters. The output
% is a single number, the optimal threshold.
%
% I - Filtered mRNA image
% N - Number of points in threshold range to use  (~100)
% smo- number of points to use in smoothing data (~4)
% Thresmin - Minimum threshold value to use
% Thresmax - Maximum threshold value to use


function f=mRNAthreshold(I,N,smo,Thresmin,Thresmax)

try
    
troubleshooting=0;


        I=single(I)/single(max(I(:))); %Scaling of I so pixels are all less than 1
        DataArray=zeros([N,2]); %Initialize array to store calculated parameters
        
        ThresholdArray=linspace(Thresmin,Thresmax,N); % Define array of threshold values

        for i=1:N  % Loop over possible threshold values
    
        ThreshImage=I>ThresholdArray(i); % Threshold image
        [LabThreshImage,NumObj]=bwlabel(ThreshImage,4);   % Label objects with 4 nearest neighbours
    
        LabPixels = LabThreshImage(ThreshImage);  %Array containing all lableled pixels. 
        [bb, m1, n1] = unique(LabPixels,'first'); % Arrays containing the ordered unique pixels with indices (first occurenc)
        [bb, m2, n2] = unique(LabPixels,'last'); % Arrays containing the ordered unique pixels with indices (last occurenc)


         DataArray(i,1) =sum(m1==m2);    %Find number of objects that consis of only 1 pixel
         DataArray(i,2)=NumObj;          %Find total number of objects
        
        end
        
FracObjGrOnePixel = (DataArray(:,2)-DataArray(:,1))./(1+DataArray(:,1)); %Calculates the fraction of objects in image that consist of more than one pixel
FracObjGrOnePixelSmooth=smooth(FracObjGrOnePixel,smo);%Smoothes this as a function of threshold

if troubleshooting

figure, plot(ThresholdArray,FracObjGrOnePixelSmooth);
hold on, plot(ThresholdArray,FracObjGrOnePixel,'r'); hold off

title('Frac Obj Greater than one pixel')

figure, plot(ThresholdArray, DataArray(:,2),'g.')

title('Number of objects')

figure, hist(DataArray(:,2),15)

end

IndexMax = find(FracObjGrOnePixelSmooth==max(FracObjGrOnePixelSmooth)); %Linear index of threshold maximum


f=ThresholdArray(IndexMax(end)); %Determine the threshold value at which the fraction of objects greater than one pixel peaks


%%%%%%%%% Threshold to 1 if there are too many objects

if max(max(bwlabeln(I>f)))>=3000 & f < 0.05
    f=1;
end

%%%%%%%%%%%%%%

catch
    f=1;
end
    
