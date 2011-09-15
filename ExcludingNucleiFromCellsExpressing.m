%%                ExcludingNucleiFromCellsExpressing.m
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
% This code excludes nuclei from and mRNA after being processed to remove
% ambigious mRNAs.
%
%
%
%


function [LmRNANew,mRNA1oIndNew, mRNA2oIndNew, mRNApIndNew, AllmRNAinNucPostNew,NucmRNA1oIndNew,NucmRNA2oIndNew,NucmRNApIndNew] = ExcludingNucleiFromCellsExpressing(LmRNA,mRNA1oInd, mRNA2oInd, mRNApInd, AllmRNAinNucPost,NucInConfReg,NucLab)


mRNAsInConfReg = find(ismember(AllmRNAinNucPost,NucInConfReg));

mRNA1oIndLeft = mRNA1oInd(ismember(mRNA1oInd,mRNAsInConfReg))';
mRNA2oIndLeft = mRNA2oInd(ismember(mRNA2oInd,mRNAsInConfReg))';
mRNApIndLeft = mRNApInd(ismember(mRNApInd,mRNAsInConfReg))';

AllOld=[mRNA1oIndLeft;mRNA2oIndLeft;mRNApIndLeft];
AllNew=[1:length(AllOld)]';

mRNA1oIndNew=[1:length(mRNA1oIndLeft)];
mRNA2oIndNew=[1+length(mRNA1oIndLeft):length(mRNA2oIndLeft)+length(mRNA1oIndLeft)];
mRNApIndNew=[1+length(mRNA1oIndLeft)+length(mRNA2oIndLeft):length(mRNA2oIndLeft)+length(mRNA1oIndLeft)+length(mRNApIndLeft)];


LmRNANew=zeros(size(LmRNA));


AllmRNAinNucPostNew=zeros([length(AllOld),1]);

LmRNA=sparse(LmRNA);

for i=1:length(AllOld)
    
    LmRNANew(find(LmRNA==AllOld(i)))=AllNew(i);
    AllmRNAinNucPostNew(i)=AllmRNAinNucPost(AllOld(i));
    
end


NucmRNA1oIndNew=zeros([max(NucLab(:))],1);
NucmRNA2oIndNew=zeros([max(NucLab(:))],1);
NucmRNApIndNew=zeros([max(NucLab(:))],1);
 
 
for i=1:length(mRNA1oIndNew)
NucmRNA1oIndNew(AllmRNAinNucPostNew(mRNA1oIndNew(i)))=NucmRNA1oIndNew(AllmRNAinNucPostNew(mRNA1oIndNew(i)))+1;
end
 

for i=1:length(mRNA2oIndNew)
NucmRNA2oIndNew(AllmRNAinNucPostNew(mRNA2oIndNew(i)))=NucmRNA2oIndNew(AllmRNAinNucPostNew(mRNA2oIndNew(i)))+1;
end 


for i=1:length(mRNApIndNew)
NucmRNApIndNew(AllmRNAinNucPostNew(mRNApIndNew(i)))=NucmRNApIndNew(AllmRNAinNucPostNew(mRNApIndNew(i)))+1;
end 
 


