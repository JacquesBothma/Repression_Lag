%%                               AmbigiousmRNAsDiff.m
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
% This code reassigns ambigious mRNA on nuclear boundaries for cells that have too many mRNAs. The region of focus is centered on cells that have too many mRNAs.
% Ambigious mRNAs connected by shared nuclei are identified as clusters and all possible mRNA assignment schemes are
% calculated and scored. This method only works efficiently for clusters
% consisting of less than 10 mRNAs and so big clusters need to be
% artificially severed. This is done in the first section of the code.
%
% The inputs:
%
% mRNA1o - logical array of mRNA1o centers
% mRNA2o - logical array of mRNA2o centers
% mRNAp - logical array of mRNAp centers
% NucLab - Segmented and labelled Nuclei
% Nuc - Raw nuclear image
% DilR - Radius of disk used to dilate nuclear edges to define ambigious
% region (1-3)
% radarea - Radius of disk used around mRNA to define neighbouring nuclei
% (1-5)
% MaxClusterSize - Largest number of mRNAs in a cluster  (<10)
% fracsplit - fraction by which ambigious mRNA clusters are split (0.5)
% genotype - 1 for heterozygous, i.e maximum of one site of transcription per
% nucleus, or 2 for homozygous. 0 for autodeterm based on median number in
% on cells.
% ACON - connectivity array of nuclei
% SrchDepth - maximum number of connections that a nucleus is considered
% within the ambigious region.

% The outputs:
%
% LmRNA - Label matrix of all mRNAs inside nulcear region
% mRNA1oInd - Indices of LmRNA that are mrna1o
% mRNA2oInd - Indices of LmRNA that are mrna2o
% mRNApInd - Indices of LmRNA that are mrnap
% AllmRNAinNucPre - Allocation of mRNA in LmRNA to nuclei, rows correpond
% to mRNA indices and entries to the nuclei that the mRNA is associated
% with.
% AllmRNAinNucPost - Allocation of mRNA in LmRNA to nuclei after processing
%
%Notes :
%
%11/12/10 - Added important modification  "idx(idx<1 |
%idx>(SIZ(1)*SIZ(2)))=[];" which partially deals with the boundary issue
%that has caused problems in some images. Makes sure that if indices are
%less than 1 and more than maximum number of pixles they are disgarded. It
%does mean though that there could be wrappin on the top and bottom of the
%image. not considered important since for good images boundary should take
%care of that. 
%


function [LmRNA,mRNA1oInd, mRNA2oInd, mRNApInd, AllmRNAinNucPre, AllmRNAinNucPost,NucmRNA1oInd,NucmRNA2oInd,NucmRNApInd,genotype]=AmbigiousmRNAsDiff(mRNA1o,mRNA2o,mRNAp,NucLab,Nuc,DilR,radarea,MaxClusterSize,fracsplit,genotype,ACON,SrchDepth)

troubleshooting=1;

if genotype == 0                          % Determining the genotype of cells by looking at the median number of mRNAs in all
    
Np= hist(double(NucLab(mRNAp>0)),double(max(NucLab(:))));
N1= hist(double(NucLab(mRNA1o>0)),double(max(NucLab(:))));
N2= hist(double(NucLab(mRNA2o>0)),double(max(NucLab(:))));

[All,AllN] =hist(Np+N1+N2,[0:max(Np+N1+N2)]);

genotype = median(nonzeros(Np+N1+N2))


All(3)/All(2)

if genotype==1 && All(3)/All(2)>0.3
    genotype=2
end
end

%%%%% Initialzing nuclear arrays and defining ambigious mRNA region as well as
%%%%% search area for nuclei

NucLab_PixelIdxList = regionprops(NucLab, {'PixelIdxList'}); % Structure containing linear indices of nuclear pixels
totalnuc=max(NucLab(:));                                     % Total number of nuclei
NucLabIndex=[1:totalnuc]';                                   % Array containing indices of nuclei
SIZ = size(NucLab);                                          % Dimensions of nuclear array

%%%%% Finding nuclei that have too many mRNAs

testnumbermRNA1o=histc(NucLab(mRNA1o>0),NucLabIndex);   
testnumbermRNA2o=histc(NucLab(mRNA2o>0),NucLabIndex);
testnumbermRNAp=histc(NucLab(mRNAp>0),NucLabIndex);

 if  genotype == 2
     NucOver=find(((testnumbermRNA1o+testnumbermRNA2o+ 2*testnumbermRNAp)>4 | (testnumbermRNA1o)>2 | (testnumbermRNA2o)>2) | (testnumbermRNA2o+testnumbermRNA1o+testnumbermRNAp)>2);      %Penalty for total of more than 4 mrnas per nucleus 
 else
     NucOver=find(((testnumbermRNA1o+testnumbermRNA2o+ 2*testnumbermRNAp)>2 | (testnumbermRNA1o)>1 | (testnumbermRNA2o)>1));
 end
 
%%%%%

NucIntensityCent=zeros([totalnuc,2]);                        % Matrix to store intensity centroid of nuclei

for i=1:totalnuc                                             % Loop to determine intensity weighted centroids of nuclei
     [r,c]=ind2sub(SIZ,NucLab_PixelIdxList(i).PixelIdxList);
     NucIntensityCent(i,1)=sum(r.*single(Nuc(NucLab_PixelIdxList(i).PixelIdxList)))/sum(single(Nuc(NucLab_PixelIdxList(i).PixelIdxList)));
     NucIntensityCent(i,2)=sum(c.*single(Nuc(NucLab_PixelIdxList(i).PixelIdxList)))/sum(single(Nuc(NucLab_PixelIdxList(i).PixelIdxList))); 
end  


    %%%%% Defining Region in which mRNAs are considered ambigious

    mebers=[];
    
    for oo=1:length(NucOver)
    
         mebers = [mebers, graphtraverse(ACON>0,NucOver(oo),'DEPTH',SrchDepth, ...
        'method','BFS','directed',false)]; %mRNAs to be split from cluster    
    end
    
    mebers=unique(mebers);
    
    LabBoundary = LabelObjectBoundary(NucLab,DilR);
    
    boundary=ismember(LabBoundary,mebers);                       %Pixels that make up the edges of the nuclei
    
    NucLabOutBndry=(~boundary).*NucLab;                %labelled nuclei without edges
    
    
    
    

if troubleshooting
    imshowbig(boundary)
end
    
    %%%%%%%%%%%%%%%%%%%%%
    
    %%%%% Determining linear indices of circular nearest neighbours
 

     T=false(2*radarea+1);                       % blank matrix of appropriate size
     T(radarea+1,radarea+1) = 1;                 % seed the blank with a 1 at the center
     T2 = imdilate(T,strel('disk',radarea,0));   % use imdilate to create disk
     [row,col] = ind2sub(size(T2),find(T2));     % get the row and column value of the non-zero indices
     CircNeib = sub2ind(SIZ,row,col) ...
         -sub2ind(SIZ,radarea+1,radarea+1);      % linear indices of neighbours
     %%%%%%%%%%%

%%%%%%%%%%%%%%%
                                   

%%%%%%  Initializing mRNA matrices by labelling individual mrna types,
%%%%%%  combining into a single mRNA matrix and then tallying how many
%%%%%%  mRNAs per nucleus

[LmRNA1o, numRNA1o] = bwlabeln(mRNA1o.*(NucLab>0));         % labeled matrix with only mRNA1 in region
[LmRNA2o, numRNA2o] = bwlabeln(mRNA2o.*(NucLab>0));         % labeled matrix with only mRNA2 in region
[LmRNAp, numRNAp] = bwlabeln(mRNAp.*(NucLab>0));            % labeled matrix with only mRNA pairs in region


LmRNA2o(find(LmRNA2o))=LmRNA2o(find(LmRNA2o))+numRNA1o;     %Incresaing labels of mRNA2 for big matrix
LmRNAp(find(LmRNAp))=LmRNAp(find(LmRNAp))+numRNA1o+numRNA2o;%Incresaing labels of mRNAp for big matrix
LmRNA=LmRNA1o+LmRNA2o+LmRNAp;                               % Label matrix that includes all mRNA types

AllmRNAinNucPre = [LmRNA(find(LmRNA)),NucLab(find(LmRNA))]; % Array specifying which mRNAs are assigned to which nuclei before 
AllmRNAinNucPre =sortrows(AllmRNAinNucPre,1);               % ambigioous mRNAs are reassigned
AllmRNAinNucPre=AllmRNAinNucPre(:,2);


mRNA1oInd=[1:numRNA1o];                                   %Big matrix lables for mRNA1
mRNA2oInd=[numRNA1o+1:numRNA1o+numRNA2o];                 %Big matrix lables for mRNA2
mRNApInd=[numRNA1o+numRNA2o+1:numRNAp+numRNA1o+numRNA2o]; %Big matrix lables for mRNAp

mRNAOnBndry=boundary.*(LmRNA>0);                       % Defining mRNAs in ambigious region
mRNAOnBndry_linearindices=find(mRNAOnBndry);           % Vector containing linear indices of ambigious mRNAs
[Rmrna,Cmrna]=ind2sub(SIZ,mRNAOnBndry_linearindices);  % Vector containing row and column locations of the ambigious mRNAs
mRNAOnBndry_indices=LmRNA(mRNAOnBndry_linearindices);  % Vector containing indices of ambigious mRNAs
Number_Ambig_mRNAs=length(mRNAOnBndry_indices);        % Number of ambigious mRNAs

AllmRNAinNucPost=AllmRNAinNucPre;                      %Initializing array specifying which mRNAs are assigned to which nuclei postproc
AllmRNAinNucPost(mRNAOnBndry_indices)=0;


NumRNA1oPerNuc_OutBndry = histc(NucLabOutBndry(find(LmRNA1o)),NucLabIndex); %Number mRNA1o per nuc not including ambigious region

NumRNA2oPerNuc_OutBndry = histc(NucLabOutBndry(find(LmRNA2o)),NucLabIndex); %Number mRNA2o per nuc not including ambigious region

NumRNApPerNuc_OutBndry = histc(NucLabOutBndry(find(LmRNAp)),NucLabIndex);  %Number mRNAp per nuc not including ambigious region

NumRNA1oPerNuc_Assign = NumRNA1oPerNuc_OutBndry;              %Initialize arrays where processed mRNAs are allocated
NumRNA2oPerNuc_Assign = NumRNA2oPerNuc_OutBndry;
NumRNApPerNuc_Assign = NumRNApPerNuc_OutBndry;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding indices of mRNAs in overpop nuclei  


OverfullmRNA=nonzeros(LmRNA(ismember(NucLab,NucOver)>0));
OverfullmRNAonbnd=intersect(mRNAOnBndry_indices,OverfullmRNA);
AmbigmRNAindex = find(ismember(mRNAOnBndry_indices,OverfullmRNAonbnd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Determining connectivity of ambigious mRNA to nuclei, ambigious mRNA
%%%% to ambigious mRNA through shared nuclei

NeigNucAmbmRNALab=[];                    %Initialize array
AmbigmRNALab = [];
 
for i=1:Number_Ambig_mRNAs
     
sindex=mRNAOnBndry_linearindices(i);             %linear index of the ith ambigious mRNAs

idx= CircNeib+sindex;                           %linear indices of neighbouring pixels

idx(idx<1 | idx>(SIZ(1)*SIZ(2)))=[];            % Line that sorts out boundary issues on left and rigth edges, still get wrapping on top and bottom

bordernuc=nonzeros(unique(NucLab(idx)));         %indices of unique nuclei that fall within neighbour radius
mrnalabel=i*ones(length(bordernuc),1);           %Vector with ambigious mRNA label 

AmbigmRNALab=[mrnalabel;AmbigmRNALab];           % Matrix to store ambig mrna label
NeigNucAmbmRNALab=[bordernuc;NeigNucAmbmRNALab]; % Matrix to store neighbouring nuclei label
 
end


 
ConMatARNANuc = sparse(accumarray([AmbigmRNALab,NeigNucAmbmRNALab],1)); %Connectivity matrix of ambigious mRNAs to nuclei
ConMatOnlyARNA = sparse(ConMatARNANuc*ConMatARNANuc'>0); %Connectivity matrix of ambigious mRNAs to each other

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%%%%%%%%%%%% Finding and breaking up of connected ambigious mRNA clusters

[numberofmRNAclusters, clusterlabel] = graphconncomp(ConMatOnlyARNA,'Directed', 'false'); %Labelling of clusters of mRNAs
numberofmRNAspercluster=histc(clusterlabel,[1:max(clusterlabel)]);                                  %How many mrnas in a cluster


%%%% Displaying clusters of ambigious mRNAs

if troubleshooting

NucDisp=zeros(size(NucLab));
cntr=1;

for i=1:numberofmRNAclusters

    mRNAinithcluster = find(clusterlabel==i);

    NucInCluster=[];

    for j=1:length(mRNAinithcluster)

        NucInCluster=[NucInCluster, find(ConMatARNANuc(mRNAinithcluster(j),:))];

    end

    NucInCluster = unique(NucInCluster);

    for ii=1:length(NucInCluster)

        NucDisp(NucLab_PixelIdxList(NucInCluster(ii)).PixelIdxList)=cntr;

    end
    
    cntr=cntr+1;

end

imshowbig(label2rgb(NucDisp,'jet',[0,0,0],'shuffle'));
hold on
[x,y] = ind2sub(size(NucDisp), find(ismember(LmRNA,mRNAOnBndry_indices(AmbigmRNAindex))));
scatter(y,x,50,'w');
hold off

end

%%%%

    %%% Loop that breaks up mRNA clusters. Finds two furthest mRNAs in
    %%% cluster breaks cluster in two pieces, size of is some fraction of
    %%% this distance


        
    while max(numberofmRNAspercluster)>MaxClusterSize

    Distancebetween=graphallshortestpaths(ConMatOnlyARNA);              %All the shortest paths
    Distancebetween(Distancebetween==Inf)=0;                            %Setting not connected nodes to 0
    clustersgreatNCmRNAs=find(numberofmRNAspercluster>MaxClusterSize);  %indices of clusters that have more than MaxClusterSize mRNAs
    numbigclust=length(clustersgreatNCmRNAs);                           %number of clusters that have more than MaxClusterSize mRNAs

    
    for jj=1:numbigclust                                            %Loop that breaks down clusters

            mRNAclust = find(clusterlabel==clustersgreatNCmRNAs(jj));   % indices of mrnas in jj'th cluster
            Allpairs = nchoosek(mRNAclust,2);                           % All possible pairs of mRNA in cluster
            alldist=Distancebetween(sub2ind(size(ConMatOnlyARNA), Allpairs(:,1), Allpairs(:,2))); %number of connections between all pairs
            maxdist=max(alldist);                                       % maximum number of connections between two points
            ind = find(alldist==maxdist,1);                             % Linear index of first member of pair of mRNAs sepearted by largest distance
            mRNAtar=Allpairs(ind,1);                                    % index of ambigious mRNA to use as starting point
            distsplit=floor(maxdist*fracsplit);                         % Number of connections that correspond to this fraction

            mRNAinA = graphtraverse(ConMatOnlyARNA,mRNAtar,'DEPTH',distsplit, ...
                'method','BFS','directed',false); %mRNAs to be split from cluster
            mRNAinB=setdiff(mRNAclust,mRNAinA); %mRNAs left in cluster

            for i=1:length(mRNAinA) %Loop that severs connections between two groups
                ConMatOnlyARNA(mRNAinA(i),mRNAinB)=0;
                ConMatOnlyARNA(mRNAinB,mRNAinA(i))=0;
            end



    end


 
    [numberofmRNAclusters, clusterlabel] = graphconncomp(ConMatOnlyARNA,'Directed', 'false'); % Relabeling of clusters after split
    numberofmRNAspercluster=histc(clusterlabel,[1:max(clusterlabel)]);                        % Number of mRNAs in cluster

    max(numberofmRNAspercluster)

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Displaying clusters of ambigious mRNAs
if troubleshooting

NucDisp=zeros(size(NucLab));
cntr=1;
for i=1:numberofmRNAclusters

    mRNAinithcluster = find(clusterlabel==i);

    NucInCluster=[];

    for j=1:length(mRNAinithcluster)

        NucInCluster=[NucInCluster, find(ConMatARNANuc(mRNAinithcluster(j),:))];

    end

    NucInCluster = unique(NucInCluster);

    for ii=1:length(NucInCluster)

        NucDisp(NucLab_PixelIdxList(NucInCluster(ii)).PixelIdxList)=cntr;

    end
    
    cntr=cntr+1;

end



imshowbig(label2rgb(NucDisp,'jet',[0,0,0],'shuffle'));
hold on
[x,y] = ind2sub(size(NucDisp), find(ismember(LmRNA,mRNAOnBndry_indices(AmbigmRNAindex))));
scatter(y,x,50,'w');
hold off

end

%%%%


%%%%%% Processing of clusters
    
 for i=1:numberofmRNAclusters
 
    mrnalabelincluster = find(clusterlabel==i); %labels of mrna in the i'th clustercluster
    
    n=length(mrnalabelincluster); % number of mrnas in the i'th cluster 
    
    
    AmbigRNAtype=ismember(mRNAOnBndry_indices(mrnalabelincluster),mRNA1oInd)+ ...
        2*ismember(mRNAOnBndry_indices(mrnalabelincluster),mRNA2oInd)+ ...
        3*ismember(mRNAOnBndry_indices(mrnalabelincluster),mRNApInd); % Classifies the mRNAs in 
    % the clusters as different types

    mult=1; %initializing number of total different number of nuclear combinations
    
    nucs=[];
    
    for ii=1:n  % loop that determines the nuclei that are associated with each mrna and with the whole cluster
        nuc = find(ConMatARNANuc(mrnalabelincluster(ii),:));
        mult= mult*length(nuc);
        nucs= [nucs,nuc];
    end
    
    Nucs=unique(nucs);             % Array of all the nuclei in the cluster without repetition
    Score=zeros([1,mult]);         % Defining array that will be used to store the score of different combinations
    DisFromNucCen=zeros([1,mult]); % Defining array that will store the sum of differences between nuclear centers and mRNA
    
    
        for ii=1:n % Loop creates matrix where the columns represent all possible assignment of mRNAs to nearest nuclei,
             ... assigned nuclei's indexes appear in columns

            if ii==1;

            AllPosmRNAAss = find(ConMatARNANuc(mrnalabelincluster(ii),:)); %Nuclei associated with ii'th ambigious mRNA

            else

            [a,b] = size(AllPosmRNAAss); 

            nuc=find(ConMatARNANuc(mrnalabelincluster(ii),:));
            
            ln = length(nuc);

            AA=[];

            for yy=1:ln
                AA=[AA,repmat(nuc(yy),1,b)];
            end

            AllPosmRNAAss=[repmat(AllPosmRNAAss,1,ln);AA];

            end

        end
    
       
  
    
    for ii=1:mult % loop that calculates the distance from mRNA to intensity weighted centroid of sharing nuclei
        
        rmrna=Rmrna(mrnalabelincluster); %row of mRNA location
        cmrna=Cmrna(mrnalabelincluster); %column of mRNA location
        rnuc = NucIntensityCent(AllPosmRNAAss(:,ii),1); %row(s) of nuclei location(s)
        cnuc = NucIntensityCent(AllPosmRNAAss(:,ii),2); %column(s) of nuclei location(s)
        dis=sqrt((bsxfun(@minus,rnuc,rmrna)).^2+  (bsxfun(@minus,cnuc,cmrna)).^2); % Distance between mRNA and each nucleus
        DisFromNucCen(1,ii)=sum(dis);           %Sum of differences in distance
   
    end
    
    for kk=1:mult % Loop that temporarily adds mRNAs to nuclei and scores different combinations
        
 
        numtem1=NumRNA1oPerNuc_Assign; %create temporary matrix to store number of mRNAs
        numtem2=NumRNA2oPerNuc_Assign;
        numtemp=NumRNApPerNuc_Assign;
        
        for ii=1:n %Assign mRNAs to nuclei
            if AmbigRNAtype(ii)==1
            numtem1(AllPosmRNAAss(ii,kk))=numtem1(AllPosmRNAAss(ii,kk))+1;
            elseif AmbigRNAtype(ii)==2
            numtem2(AllPosmRNAAss(ii,kk))=numtem2(AllPosmRNAAss(ii,kk))+1;
            else
            numtemp(AllPosmRNAAss(ii,kk))=numtemp(AllPosmRNAAss(ii,kk))+1;
            end
        end
        
 
   
        

 
     ntem1N=numtem1(Nucs);                        %Number of mRNA1s in appropriate nuclei
     ntem2N=numtem2(Nucs);                        %Number of mRNA2s in appropriate nuclei
     ntempN=numtemp(Nucs);                        %Number of mRNAps in appropriate nuclei


     S1=zeros(size(length(ntempN)));

 if  genotype == 2
     
     S1=S1 + 100*sum((ntem1N+ntem2N+2*ntempN)>4);      %Penalty for total of more than 4 mrnas per nucleus 
     S1=S1 + 10*sum(ntem1N>2);                         %Penalty for more than 2 mrna1s per nucleus
     S1=S1 + 10*sum(ntem2N>2);                         %Penalty for more than 2 mrna2s per nucleus
     S1=S1 + 10*sum(ntempN>2);                         %Penalty for more than 2 mrnaps per nucleus
     S1=S1 + 100*sum((ntem1N+ntem2N+ntempN)>2);

 
 else
     
     S1=S1 + 100*sum((ntem1N+ntem2N+2*ntempN)>2);      %Penalty for total of more than 2 mrnas per nucleus 
     S1=S1 + 10*sum(ntem1N>1);                         %Penalty for more than 1 mrna1s per nucleus
     S1=S1 + 10*sum(ntem2N>1);                         %Penalty for more than 1 mrna2s per nucleus
     S1=S1 + 10*sum(ntempN>1);                         %Penalty for more than 1 mrnaps per nucleus
     
 end

     Score(:,kk)=S1;

    end
  
 
%%%%%%%% Finding best configuration

    minscore = find(Score==min(Score)); %lowest score
    colval = find(DisFromNucCen==min(DisFromNucCen(minscore)),1); %lowest distance from center

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
%%%%% Assigning ambigious mRNAs to nuclei

      for ii=1:n 

          AllmRNAinNucPost(mRNAOnBndry_indices(mrnalabelincluster(ii))) = AllPosmRNAAss(ii,colval);
          
            if AmbigRNAtype(ii)==1
            NumRNA1oPerNuc_Assign(AllPosmRNAAss(ii,colval))=NumRNA1oPerNuc_Assign(AllPosmRNAAss(ii,colval))+1;
            elseif AmbigRNAtype(ii)==2
            NumRNA2oPerNuc_Assign(AllPosmRNAAss(ii,colval))=NumRNA2oPerNuc_Assign(AllPosmRNAAss(ii,colval))+1;
            else
            NumRNApPerNuc_Assign(AllPosmRNAAss(ii,colval))=NumRNApPerNuc_Assign(AllPosmRNAAss(ii,colval))+1;
            end
      end
 
     
 end

NucmRNA1oInd=zeros([max(NucLab(:))],1);
NucmRNA2oInd=zeros([max(NucLab(:))],1);
NucmRNApInd=zeros([max(NucLab(:))],1);
 
 
for i=1:length(mRNA1oInd)
NucmRNA1oInd(AllmRNAinNucPost(mRNA1oInd(i)))=NucmRNA1oInd(AllmRNAinNucPost(mRNA1oInd(i)))+1;
end
 

for i=1:length(mRNA2oInd)
NucmRNA2oInd(AllmRNAinNucPost(mRNA2oInd(i)))=NucmRNA2oInd(AllmRNAinNucPost(mRNA2oInd(i)))+1;
end 


for i=1:length(mRNApInd)
NucmRNApInd(AllmRNAinNucPost(mRNApInd(i)))=NucmRNApInd(AllmRNAinNucPost(mRNApInd(i)))+1;
end 
 
 
 
 