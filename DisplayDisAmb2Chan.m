%%                              DisplayDisAmb2Chan.m
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
%
% This code dipslays output of AmbigiousmRNAsDiff in a representation that
% makes it easier to see what it did. This includes labelling where nuclei
% have lost mRNAs and where they have been moved to
%








function DisplayDisAmb2Chan(LmRNA,mRNA1oInd, mRNA2oInd, mRNApInd, AllmRNAinNucPre, AllmRNAinNucPost,NucmRNA1oInd,NucmRNA2oInd,NucmRNApInd,handles)

RegionPropNuc=regionprops(handles.LabelNucSeg,'PixelIdxList');

LabBoundary = LabelObjectBoundary(handles.LabelNucSeg,1);

disambmRNA=find((AllmRNAinNucPre-AllmRNAinNucPost)~=0);
RNANucBefore=zeros(1,max(handles.LabelNucSeg(:)))';
RNANucAfter=zeros(1,max(handles.LabelNucSeg(:)))';
 
 
for i=1:length(AllmRNAinNucPre);
    
    if ismember(i,mRNApInd)
    RNANucBefore(AllmRNAinNucPre(i))=RNANucBefore(AllmRNAinNucPre(i))+2;
    RNANucAfter(AllmRNAinNucPost(i))=RNANucAfter(AllmRNAinNucPost(i))+2;
    else
    RNANucBefore(AllmRNAinNucPre(i))=RNANucBefore(AllmRNAinNucPre(i))+1;
    RNANucAfter(AllmRNAinNucPost(i))=RNANucAfter(AllmRNAinNucPost(i))+1;    
    end
    
    
end


Increase=RNANucAfter-RNANucBefore>0;
Decrease=RNANucAfter-RNANucBefore<0;
IncDec=ismember(LabBoundary,find(Increase))*2+ismember(LabBoundary,find(Decrease));

blank = (IncDec>0) + (LmRNA>0); 

ZZ2=uint8(handles.LabelNucSeg>0);

for i=1:max(handles.LabelNucSeg(:))
    ZZ2(RegionPropNuc(i).PixelIdxList)=RNANucAfter(i)+ZZ2(RegionPropNuc(i).PixelIdxList);
end

NucDis = label2rgbBackdrop(ZZ2,'jet',[0,0,0],handles.NucChan.*uint8(~blank));
EdgeDis= label2rgb(IncDec.*(~LmRNA>0),[1,0,1;0,1,1],[0,0,0]);
mRNADis= label2rgb(ismember(LmRNA,mRNA1oInd),[1,0,0],[0,0,0]) + label2rgb(ismember(LmRNA,mRNA2oInd),[0,1,0],[0,0,0])+label2rgb(ismember(LmRNA,mRNApInd),[1,1,0],[0,0,0]);

imshowbig(NucDis+EdgeDis+mRNADis)

hold on,

 [x,y] = ind2sub(size(handles.NucChan), find(ismember(LmRNA,intersect(disambmRNA,mRNA1oInd))));
 scatter(y,x,50,'r');
 
 
 [x,y] = ind2sub(size(handles.NucChan), find(ismember(LmRNA,intersect(disambmRNA, mRNA2oInd))));
 scatter(y,x,50,'g');
 
 
 [x,y] = ind2sub(size(handles.NucChan), find(ismember(LmRNA,intersect(disambmRNA,mRNApInd))));
 scatter(y,x,50,'y');