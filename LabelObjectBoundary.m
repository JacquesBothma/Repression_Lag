%%                             LabelObjectBoundary.m
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
%
% This code determines the boudaries of objects in a labelled image. The
% thickness of the boundary can be arbitrarily selected.
%
% Input:
%
% LabelMatrix -  Label matrix of objects
% thicknes -  thicknes of boundary in pixels
%
% Output:
%
% LabBoundary - Labeled boundary pixels 
%
% Comments:

function  LabBoundary = LabelObjectBoundary(LabelMatrix,thicknes)

LabelMatrix=cast(LabelMatrix,'uint16');  % Change variable type to uint16 to speed up calculations

    Boundary=ones(size(LabelMatrix))>0;  % Initialize array containing edges of nuclei
    S=[0 1; 1 0; -1 0; 0 -1];            % Shift matrix which covers four nearest neighbours

    for i=1:4                           
    TEMP=(LabelMatrix==circshift(LabelMatrix,S(i,:)));
    Boundary=Boundary&TEMP;
    end

    if thicknes==1
    LabBoundary=LabelMatrix.*uint16(~Boundary); % Pixels that make up the edges of objects
    else
     Bn=(~Boundary);
        for i=1:thicknes-1
        Bn=imdilate(Bn,strel('square',3));
        end
    LabBoundary = LabelMatrix.*uint16(Bn);
    end
    
  LabBoundary=cast(LabBoundary,class(LabelMatrix));
    