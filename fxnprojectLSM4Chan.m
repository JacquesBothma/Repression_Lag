%%                             fxnprojectLSM4Chan.m
%
% Jacques Bothma & Alistair Boettiger                     Last Modified: 09/13/11     
% Levine Lab, UC Berkeley                       
% Functionally complete                          
%
%
%% Attribution:
%  Feel free to use, modify and distribute this code provided that you
%  attribute Jacques Bothma and Alistair Boettiger for development.
%  This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
%  To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.
%%
% Overview:
%
%  This code takes a stack of fluorescence confocal images of a Drosophila embryo
%  that has been stained for DNA and up to three different types of mRNA
%  and projects this onto a 2D image. The code allows for the user to
%  select whether it is a simple max project or whether the mRNA layers are
%  filtered before max projetcing to slect for only the cores of bright
%  pixels. It also finds the region that corresponds to main embryo on in
%  the field of view.

% Input:
%
% fin - folder name
% name - mat file with parsed info for the lsm file
% chns - Channels that have data in them and how they need to be ordered in
% final file.
% subsets - Subsets of channels to be used in max projections (between 0
% and 1)
% maxProj - Channels to be max projected
% LogFilter - Channels to be filtered and then max projected
% FilterParm - Paramters to use in filtering mRNA
% nucchan - The channel that contains the nuclear stain
% objT - Number between 0 and 1 that defines how embryo region is defined.
% imnum - Embryo (stack) number to be processed
%
%
%
%
% Output:
% I - An array with the folowing sustructures:
% I.Io - The projected image
% I.IoF - The first image of the stack
% I.IoL - The last image of the stack
% I.NucBright - The brightest layer of nuclei
% I.bw - A black and white image showing the location of the main embryo
%
%
% Comments:
%
%
%
%
%



function I = fxnprojectLSM4Chan(fin,name,chns,subsets,maxProj,LogFilter,FilterParm,nucchan,objT,imnum)


%Read in image stack
stack = loadlsm([fin '\' name],imnum);

%Determine number of layers in stack

frames= length(stack);

subsetslayer=ceil(frames*subsets);

disp(['fxnproject: ' num2str(imnum)]);

watcher = 0; clear Im;


%%%%%% Store image data in structure


    for f = 1:frames        
        
        for i=1:length(chns)
            
            try  % try in case the file does not exist
                
                 Im{i}{f} = cell2mat(stack{1,f}(i));
                 
            catch     % stop iterating if exhausted all files.
                disp(['error, cannot find file ', file{i,f}]);
                watcher = 1;
            end
        end
          if watcher == 1;
          disp('Loop broken to populate m!!!')
              break; 
          end;      
  
    end

%%%%%% Store image data in structure    


    [imsize1, imsize2] = size(Im{1}{1}); % get image size
    
        Io = zeros(imsize1,imsize2,length(find(chns)),class(Im{1}{1}));  % blank to store projected data


%%%%%% Find brightest nuclear layer in order to find embryo

IntensityNuc=zeros(frames,1);

disp([' Number of frames: ' num2str(frames)])

for i=1:frames
    IntensityNuc(i)=mean(mean(Im{nucchan}{i}));
end


[dummy,layerbrigth] = max(IntensityNuc);        

NucBright = Im{nucchan}{layerbrigth};

%%%%%% Find embryo region. This can be done using the EmbryoFindNew
%%%%%% function (objT=0), with a straight threshold (1>objT>0) on the nuclear channel or 
%%%%%% by having the user seelct it (objT=1).
       
if objT == 0 % loop that uses EmbryoFind to select region

    bw = EmbryoFindNew(NucBright,0.5,150,30,30,0.5,0.5*0.8,6,6,1/8,0.1,0.5);

elseif objT<1

    outims=EveningFilter(NucBright,100,60);

    bw = im2bw(outims,graythresh(nonzeros(outims))); %Apply user chosen threshold
    bw=imfill(bw,'holes');

else
    %allow manual area select for objT

    bw=ManualAreaSelect(1,NucBright);

end


EmbPixels=find(bw);
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert into stack and max project stack


chnswithdata=find(chns);

for i=1:length(chnswithdata);

Intensity.(['Ch' num2str(chnswithdata(i))])=zeros(frames,1);

end

ChanNotNuc=setdiff(chnswithdata,nucchan);

for i=1:frames

    Intensity.(['Ch' num2str(nucchan)])(i)=mean(Im{nucchan}{i}(EmbPixels));

end


[dummy,layerbrigthCorr] = max(Intensity.(['Ch' num2str(nucchan)]));    


if  subsetslayer(nucchan,1)-subsetslayer(nucchan,2)==0

    top = layerbrigthCorr-5;
    bottom = layerbrigthCorr;
    
     if top<=0
         top=1;
     end
     
     
    subsetslayer(nucchan,1)=top;
    subsetslayer(nucchan,2)=bottom;

end 


%%%%%%%%% First and last images %%%%%%%%%%%%%%


IoF = zeros(imsize1,imsize2,length(find(chns)),class(Im{1}{1}));
IoL = zeros(imsize1,imsize2,length(find(chns)),class(Im{1}{1}));

for i=1:length(chnswithdata)
    
    IoF(:,:,i)=Im{chnswithdata(i)}{1};
    IoL(:,:,i)=Im{chnswithdata(i)}{frames};
    
end

IoF=imresize(IoF,0.5,'nearest');
IoL=imresize(IoL,0.5,'nearest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
   for i=1:length(chnswithdata)
    
          if ismember(chnswithdata(i),find(LogFilter));
       
           
           Blank=zeros(imsize1,imsize2,2,'single');
           
           H = -fspecial('log',FilterParm(1),FilterParm(2));

           
           for j=subsetslayer(chnswithdata(i),1):subsetslayer(chnswithdata(i),2)
               
               Blank(:,:,1) = imfilter(single(Im{chnswithdata(i)}{j})-mean(nonzeros(single(Im{chnswithdata(i)}{j}))),H,'replicate');
               Blank(:,:,2) = max(Blank,[],3);

           end
           
           
           
           MaxProjChanResFilt = max(Blank,[],3);

           MaxProjChanResFilt= MaxProjChanResFilt./mean(MaxProjChanResFilt(MaxProjChanResFilt>(max(MaxProjChanResFilt(:))*FilterParm(3))));

           if strcmp(class(Im{1}{1}),'uint8')
           MaxProjChanResFilt=MaxProjChanResFilt*255;    
           else strcmp(class(Im{1}{1}),'uint16')
           MaxProjChanResFilt=MaxProjChanResFilt*(2^16-1);               
           end
               
           Itm{i}=cast(MaxProjChanResFilt,class(Im{1}{1}));
                          
           elseif ismember(chnswithdata(i),find(maxProj))
           
           subsetslayer(chnswithdata(i),1):subsetslayer(chnswithdata(i),2);    
               
           maxreshape=cell2mat(Im{chnswithdata(i)}(subsetslayer(chnswithdata(i),1):subsetslayer(chnswithdata(i),2)));
           
           maxreshape = reshape(maxreshape,imsize1,imsize2,length(Im{chnswithdata(i)}(subsetslayer(chnswithdata(i),1):subsetslayer(chnswithdata(i),2))));
       
           Itm{i} = max(maxreshape,[],3);

           maxreshape=[];
       

       end
                                     
    Io(:,:,i) = Itm{i};
   end

I.Io=Io;
I.IoF=IoF;
I.IoL=IoL;
I.NucBright=NucBright;
I.bw=bw;


end






