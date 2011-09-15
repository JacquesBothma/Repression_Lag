%%                             fxnclean4Chan.m
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
% This code takes a projected confocal stack, the image that defines where
% the embryo is and then re-orients the embryo such that the major axis is
% horizontal and crops excess area. The embryo is then oriented
% automatically so that if descernible ventral is down and anterior is to
% the left, this works ok, but there is an option to do this manually. The edge of the cropped area is
% smoothed with a gaussian blur to avoid sharp edges that will interfere with
% segmentation later.
%
% Input:
%
% I - projected image stack from fxnprojectLSM4chan
% emb_chn - channel that contains nuclear staining
% bw - black and white image 
% NucBright - Brightest layer of nuclei, used in automatic orientation
% determination
% ManualOrientation - Logical one will allow for manual orientation
%
%
% Output:
%
%
%
% Comments:
%
%
%


function I3 = fxnclean4Chan(I,emb_chn,bw,NucBright,ManualOrientation)    
           
       
      % choose largest object, remove minors, then orient
       L = bwlabel(bw,4);   % create label matrix of objects in image
  
       imdata = regionprops(L,'Area', 'Orientation'); % measure properties of objects
       keep = find([imdata.Area]==max([imdata.Area])); % keep largest object
       
       theta = [imdata.Orientation]; % orientation of objects
       thetaK=theta(keep);
       
       obj_loc = ismember(L,keep); 
               
            H = fspecial('gaussian',100,15); % Filter Kernel
            outims = imfilter(single(obj_loc),H,'replicate'); %Apply Filter
            
            
            
        if strcmp(class(I(:,:,emb_chn)),'uint8')
            
             I2=uint8(bsxfun(@times,single(I),outims));       
             
        elseif strcmp(class(I(:,:,emb_chn)),'uint16')

            I2=uint16(bsxfun(@times,single(I),outims));

           else
            disp('Error! unknown input format!')
        end
         
        [w,l,h]=size(I);

        Sp=zeros([w,l,h+2],class(I2));
         
        Sp(:,:,1:h)=I2;
        Sp(:,:,h+1)=obj_loc;
        Sp(:,:,h+2)=NucBright;
        
        Sp=imrotate(Sp,360-thetaK,'bilinear');


        
        I2=Sp(:,:,1:h);
        obj_loc=Sp(:,:,h+1);
        NucBright=Sp(:,:,h+2);

         
         L = bwlabel(obj_loc,8); 
         imdata = regionprops(L,'BoundingBox'); % compute boundaries
         c = round([imdata.BoundingBox]);
         
         try
         I3 = I2(c(2):c(4)+c(2),c(1):c(3)+c(1),:);   
         SpR= obj_loc(c(2):c(4)+c(2),c(1):c(3)+c(1),:);   
         NucBright= NucBright(c(2):c(4)+c(2),c(1):c(3)+c(1),:);   
         
         catch
             I3 = I2(c(2):c(4)+c(2)-1,c(1):c(3)+c(1)-1,:);
             SpR = obj_loc(c(2):c(4)+c(2)-1,c(1):c(3)+c(1)-1,:); 
             NucBright = NucBright(c(2):c(4)+c(2)-1,c(1):c(3)+c(1)-1,:);    
         end

         imshowbig(NucBright)

         
I3=EmbryoOrientation(SpR,I3,NucBright,ManualOrientation); %Reorient image of the embryo so the anterior part
% of the embryo (the head) is to the left and the ventral part of the
% embryo (the stomach) is to the bottom.
            
            
        
            