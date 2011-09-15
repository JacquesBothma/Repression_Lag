
%%                              conlabelobj.m
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
% This code determines the connectivity matrix of labelled objects in a 2D
% matrix. Works by shifting label matrix in eather 4 or 8 directions and
% looks at difference, pixels that have a non-zero difference touch other
% objects. Each such touch is counted and used to generate connectivity
% matrix where the size of the entry corresponds to twice the number of
% touches. Can specify whether connectivity is through sides only (pixels have 4
% connected neighbouring pixles) or through diagonals as well (pixels have
% 8 connected neighbouring pixles).
%
% Inputs & outputs: 
% The inputs are the label matrix and the number of neighbours. The
% output is the connectivity matrix.
%
% Lab -  Label Matrix
% Neighb - Number of neighbours 


function ACON=conlabelobj(Lab,Neighb)

if ~(Neighb==4 || Neighb==8)
    error('Error, number of neighbours needs to be either 4 or 8')
end

%%% Removing boundary pixels to avoid breakdown

[l,w]=size(Lab);

framepixels=[1:l,w*l:-1:w*(l-1),2*l:l:l*(w-1),l+1:l:l*(w-2)+1];

Lab(framepixels)=0;

%%%%%%%%%%


UU = max(Lab(:))+1; %Determing the number of objects in lablel matrix

if UU<=2^7-1 %Loop that reassigns variabletype of labelmatrix to speed up calculations
    Lab=cast(Lab,'int8');
elseif UU<=2^15-1
    Lab=cast(Lab,'int16');
elseif  UU<=2^31-1
    Lab=cast(Lab,'int32');
else
    Lab=cast(Lab,'single');
end


ImYSize=size(Lab,1); %Y length of image, used in linear indices.
ShiftVal=[0 1; 1 0; -1 0; 0 -1; 1 1; -1 1; 1 -1; -1 -1]; %Shift dimensions
NeighLinearInd=[ImYSize,1,-1,-ImYSize,ImYSize+1,+ImYSize-1,-ImYSize+1,-ImYSize-1]; %linear indices of the 8th closest neighbours of the 0th index.

NA=[]; %initiate arrays to store linear indices of pixels that have non-zero difference
NB=[];

LabNew=Lab+1; %Define new label matrix where labels increased by 1 to let background be labeled as 1



for i=1:Neighb %loops over different neigbour orientations to determine connectivity
    NATT=(LabNew-circshift(LabNew,ShiftVal(i,:))); %difference between normal and shifted matrix
    NAT=find(NATT); %finds linear indices of pixels where difference is non-zero
    NA=[NA;NAT];                     %Stores index of pixel that is made non-zero by subtraction
    NB=[NB;NAT-NeighLinearInd(i)];   %Stroes index of pixel that was subtracted away
end

NAV=LabNew(NA); %index of objects at the edgepixels that are connected
NBV=LabNew(NB); 

NN=[NAV,NBV]; %array that stores index pairs that represents a connection

ACONN=zeros(max(Lab(:))); %Zero matrix for connectivity array
ACONT = accumarray(NN,1); %Command that builds connectivity array from connections
ACONN(1:length(ACONT)-1,1:length(ACONT)-1)=ACONT(2:end,2:end); %store connectivity array in intialized matrix and remove connections to backround

ACON=sparse(ACONN);
