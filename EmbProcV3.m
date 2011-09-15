
%%                              EmbProcV3.m
% Jacques Bothma & Joe Magliocco                  Last Modified: 09/13/11     
% Levine Lab, UC Berkeley                       
% Functionally complete                          
%
%
%% Attribution:
% Feel free to use, modify and distribute this code provided that you
% attribute Jacques Bothma & Joe Magliocco for development.
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.
%
%
%
%  Important Notes:
% % This version written for Windows.
%
% Overview:
%  This gui loads double stained images of nuclei and mRNA nascent
%  transcripts and then using semi-automated segmentation techniques
%  locates nuclei and sites of nascent transcription.
%
%  It is a step based approach where all successive filtering and
%  segmentation steps are grouped into individual steps that allow for user
%  input. There is also a an autocycle mode in which the standard
%  parameters are applied blindly. Some parameter vlaues are stored in
%  external files while most of the image data is stored in the hadles structure of
%  the gui figure.
%
% Inputs can be chosen to be either a tif image or a mat file.
%
%
%

function varargout = EmbProcV3(varargin)

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EmbProcV3_OpeningFcn, ...
                   'gui_OutputFcn',  @EmbProcV3_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



function EmbProcV3_OpeningFcn(hObject, eventdata, handles, varargin)

    % This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




handles.output = hObject; 
  
  % Some initial setup 
    handles.step = 1;  % Starting step
   
    handles.datafolder = 'C:\Users\bothma\Documents\MATLAB\Repression_Lag\ParameterFiles'; % Specifies the folder where the Paramater files are kept
   
    handles.stepsorder=['a';'b';'c';'d';'e';'f';'g';'h';'i';'s']; %Specifies the steps that are executed and their order. 
    % Defining steps in this way allows the relative order of steps to be easily changed and
    % s
    
    %%%% Description of steps
    
    % a - Load Data into script
    % b - Segment Cores of Nuclei
    % c - Define Nulcear regions
    % d - Allows region of high confidence in segmentation to be defined
    % e - Finding mRNA of type 1
    % f - Finding mRNA of type 2
    % g - finding mRNA pairs
    % h - Disambiguate mRNA
    % i - Display segmented image
    % s - Save data

    
    
    %%%%%% Structure that contains the labels of the inputs at varios
    %%%%%% stages.
    
handles.steplabels.a =[{'Nuclei channel'},{'mRNA1 channel'},{'mRNA2 channel'},{'mRNA3 channel'},{'mat of tif'},{''},{''},{''}];
handles.steplabels.dir.a = 'Loading of data file.';

handles.steplabels.b =[{'Size of Filter Kernel'},{'Radius of filter'},{'Threshold'},{''},{''},{''},{'Evening Filter Size'},{'Radius of Evening Filter'}];
handles.steplabels.dir.b = 'Core segmentation of nuclei.';    

handles.steplabels.c =[{'Size of Blurring Kernel'},{'Radius of Blurring kernel'},{'Complex average'},{'Expand Pixels'},{'Contract Pixels'},{'Expand?'},{''},{''}];
handles.steplabels.dir.c = 'Defining Nuclear regions ';        

handles.steplabels.d =[{'Size Errodion Disk'},{''},{''},{''},{''},{''},{''},{'Manual Area Select'}];
handles.steplabels.dir.d = 'Defining Region of high confidence in nuclear segmentation';        

handles.steplabels.e =[{'Size of LOG filter'},{'Radius of LOG filter'},{'Threshold'},{'Number points for opt'},{'Smooth'},{'ThreshMin'},{'ThreshMax'},{''}];
handles.steplabels.dir.e = 'Segmenting mRNA1';     

handles.steplabels.f =[{'Size of LOG filter'},{'Radius of LOG filter'},{'Threshold'},{'Number points for opt'},{'Smooth'},{'ThreshMin'},{'ThreshMax'},{''}];
handles.steplabels.dir.f = 'Segmenting mRNA2';     

handles.steplabels.g =[{'Shell Radius'},{''},{''},{''},{''},{''},{''},{''}];
handles.steplabels.dir.g = 'Finding mRNA pairs';

handles.steplabels.h =[{'Radius of ambigious region'},{'Radius around mRNA'},{'Max Number mRNA in cluster'},{'Frac split'},{'Genotype'},{'Search Depth'},{''},{''}];
handles.steplabels.dir.h = 'Disambiguiting mRNA';

handles.steplabels.i =[{''},{''},{''},{''},{''},{''},{''},{''}];
handles.steplabels.dir.i = 'Display Nuclei';

handles.steplabels.s =[{''},{''},{''},{''},{''},{''},{''},{''}];
handles.steplabels.dir.s = 'Saving Data';

warning('off','Images:initSize:adjustingMag') % Suppress warning about magnification

    handles.totstepnum = length(handles.stepsorder); %stores the total number of steps
    

    set(handles.stepnum,'String',handles.step); % change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
    setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels
    
     try %tries to load the name of the last processed image if it exists. this gets saved in step 1. 
     
     load([handles.datafolder '\LastNamesEmbProcV3.mat']);
     
     set(handles.source,'String',namestore.source);
     set(handles.froot,'String',namestore.froot);
     set(handles.DataFile,'String',namestore.savename);
     
     catch
         
     end
     

% Update handles structure
guidata(hObject, handles);



%=========================================================================%
%                          Primary Analysis Section                       %      
%=========================================================================%
% All of the functional processing script is called in this section.
% Where possible the actual processing is done by subroutines in other
% matlab m files

function run_Callback(hObject, eventdata)
    
handles = guidata(hObject);
step = handles.step;


% Step 1 (a): Load Image Data

if handles.stepsorder(step) == 'a';



SaveParm(hObject, eventdata, handles); % Function below that saves input parameters to handles.
    
NucChanInd = str2double(get(handles.in1,'String'));    %  Index of nuclear channel
mRNA1ChanInd = str2double(get(handles.in2,'String'));  %  Index of first mRNA channel
mRNA2ChanInd = str2double(get(handles.in3,'String'));  %  Index of second mRNA channel
mRNA3ChanInd = str2double(get(handles.in4,'String'));  %  Index of third mRNA channel

MatorTif = get(handles.in5,'String');                  %  Specifies input file type
    
handles.fin = get(handles.source,'String');            % input folder
handles.fname = get(handles.froot,'String');           % embryo name
handles.emb = get(handles.embin,'String');             % embryo number


      % full file name
      handles.filename = [handles.fin,'\',handles.fname,'_',handles.emb,'.mat'];
      handles.savedata.filename= handles.filename;
      
      
    try % Attempt to load input file
       
        
        if strcmp(MatorTif,'mat') % Loop that loads either tif or mat
            
        load(handles.filename);
       
        handles.It = Data.Image;
        
        handles.NucChanO=handles.It(:,:,NucChanInd); % Oringinal Nuclear channel that won't be modified
        handles.NucChan=handles.It(:,:,NucChanInd);
        handles.mRNA1Chan=handles.It(:,:,mRNA1ChanInd);
        handles.mRNA2Chan=handles.It(:,:,mRNA2ChanInd);
       
        else
            
           I= imread([handles.filename(1:end-3), 'tif']);
           
        handles.NucChanO=I(:,:,NucChanInd);
        handles.NucChan=I(:,:,NucChanInd);
        handles.mRNA1Chan=I(:,:,mRNA1ChanInd);
        handles.mRNA2Chan=I(:,:,mRNA2ChanInd);
        
        end
        
       
       
       
      
       try % Loads the size data that is normally stored in the mat file
           
       handles.VoxelsizeX=Data.VoxelsizeX;
       handles.savedata.VoxelsizeX=handles.VoxelsizeX;
       handles.VoxelsizeY=Data.VoxelsizeY;
       handles.savedata.VoxelsizeY=handles.VoxelsizeY;
       
       catch
       
       % It there is a problem with reading in the actual value default to
       % 1
           
       handles.VoxelsizeX=1;
       handles.savedata.VoxelsizeX=1;
       
       handles.VoxelsizeY=1;
       handles.savedata.VoxelsizeY=1;
       
       disp('Voxel size set to 1')
           
       end
       
       
       %%%% Formatting display image. Default is to dislpay nuclei as blue,
       %%%% mRNA1 red, mRNA 2 green and mRNA 3 in white, if it exists.
       
       if ~isnan(mRNA3ChanInd)
       handles.mRNA3Chan=handles.It(:,:,mRNA3ChanInd);
       end
       
       DISPT=zeros(size(handles.NucChan),class(handles.NucChan));
       
       DISPT(:,:,1)=handles.mRNA1Chan;
       DISPT(:,:,2)=handles.mRNA2Chan;
       DISPT(:,:,3)=handles.NucChan;
       
       if ~isnan(mRNA3ChanInd)
           
       Dis1=label2rgbBackdrop([],[1,1,1],[0,0,0],handles.mRNA3Chan);
       imshowbig(DISPT+Dis1)
       
       else
              imshowbig(DISPT)       
       end


       
       
       title(handles.emb);
       
    catch   % display error message if there is a problem loading the file
        dir1 = ['Error loading file ', handles.filename];
        set(handles.directions,'String',dir1);
        
    end
 
% Section that reloads saved data if image has laready been processed.    

       savename = get(handles.DataFile,'String');
       savename = ([savename '_' handles.emb '_EmbProcV3Data']);
      
     if exist([handles.fin '\' savename '.mat'],'file')>0
         

         
     load([handles.fin '\' savename]);
     han=fieldnames(data);
        
     for i=1:length(han)
          handles.(char(han(i)))=data.(char(han(i))); %loading of old data into handles
          handles.savedata.(char(han(i)))=data.(char(han(i))); %loading of
        %  old data into savedata so it is carried through to next saving
        %  unless overwritten
     end

        dir1 = 'Image has already been processed and data has been loaded';
        set(handles.directions,'String',dir1);
        
     end 

%%%%% This saves the name information of the file currently being loaded so nect time the gui is started it defaults to that name.     
namestore=[];
namestore.source=get(handles.source,'String');
namestore.froot=get(handles.froot,'String');
namestore.savename = get(handles.DataFile,'String'); 
save([handles.datafolder '\LastNamesEmbProcV3.mat'],'namestore')
     
end


% Step b : Segmenting Nuclear cores
if handles.stepsorder(step) == 'b';

    
SaveParm(hObject, eventdata, handles); % Saves input paramteres to handles
    
sizefilt = str2double(get(handles.in1,'String'));  %  Size of the matrix with laplacian of gaussian kernel
radiusfilt = str2double(get(handles.in2,'String'));  %  Radius of the laplacian of gaussian kernel, if chosen to be zero uses automatically determind value
coursefilterthreshold = str2double(get(handles.in3,'String'));  % Threshold value for nuclear segmentation 

SizGaus = str2double(get(handles.in7,'String'));  % Size of gaussian kernel to use in evening filter  
RadGaus = str2double(get(handles.in8,'String'));  % Radius of gaussian kernel to use in evening filter  

%%%% Finding region of image that contains embryo %%%%%%%%
handles.EmbryoReg=EmbryoFindPostProc(handles.NucChanO,0.5,30,100,40);
handles.savedata.EmbryoReg=handles.EmbryoReg;
imshowbig(handles.EmbryoReg)
title('Embryo Region')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Performing Evening Filter %%%%%
handles.NucChan=EveningFilter(handles.NucChanO,SizGaus,RadGaus);
handles.savedata.NucChanFilt = handles.NucChan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Automatically determine best size of filter for nuclear segmentation  %%%%%
if sizefilt==0
radiusfilt = OptimalFilterRadiusNuclearSeg(handles.NucChan,1,10,100,0.25,80,1.5);
sizefilt = round(8*radiusfilt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


imshowbig(handles.NucChan)

%%%% Perform nuclear segmentation  %%%%%
handles.NucCoreSeg=segmentnucleiCore(handles.NucChan,sizefilt,radiusfilt,coursefilterthreshold);
%%%%%%%%%%%%%%%%%%%

handles.NucCoreSeg=sparse(handles.NucCoreSeg);

handles.savedata.NucCoreSeg = handles.NucCoreSeg;

%%%%%%%% Label Nuclei %%%%%%%%%%%%%%%%%%%%%%
LabelNucCoreSeg=bwlabeln(full(handles.NucCoreSeg));
%%%%%%%%%%%%%%%%%%%%%%%

imshowbig(label2rgb(LabelNucCoreSeg,'jet',[0,0,0],'shuffle'))

handles.savedata.LabelNucCoreSeg=uint16(LabelNucCoreSeg);

end

% Step c : Defining Nuclear Regions and connectivity of these regions

if handles.stepsorder(step) == 'c';

SaveParm(hObject, eventdata, handles);


sizegaus = str2double(get(handles.in1,'String'));  %
radiusgaus = str2double(get(handles.in2,'String'));  %  
averagefiltsize = str2double(get(handles.in3,'String'));  %  

Mthick = str2double(get(handles.in4,'String'));  %  
Mthin = str2double(get(handles.in5,'String'));  %  

ExpandorWater = str2double(get(handles.in6,'String'));  %  

%%%%%% Define nuclear regions
if ExpandorWater
handles.LabelNucSeg = segmentnucleiExpand(handles.NucChan,full(handles.NucCoreSeg),sizegaus,radiusgaus,handles.EmbryoReg,averagefiltsize,Mthick,Mthin);
else
handles.LabelNucSeg=segmentnucleiWatershed(handles.NucChan,full(handles.NucCoreSeg),sizegaus,radiusgaus,averagefiltsize);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Determine connectivity matrix (which nuclei are next to which
%%%%%%%%%% other nuclei) of all nuclei
handles.ACON=conlabelobj(handles.LabelNucSeg,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


handles.RegionPropNuc=regionprops(handles.LabelNucSeg, {'PixelIdxList'}, {'Centroid'});

handles.savedata.LabelNucSeg=uint16(handles.LabelNucSeg);
handles.savedata.ACON=sparse(handles.ACON);
handles.savedata.RegionPropNuc=handles.RegionPropNuc;


end


% Step d : Defining Region of Confidence in segmentation and boundary
% region as well as log 2 of nuclear density

if handles.stepsorder(step) == 'd';
    
SaveParm(hObject, eventdata, handles);
    
SizeFilt = str2double(get(handles.in1,'String'));  %  
Manarea = str2double(get(handles.in8,'String'));  % 


if Manarea==1 %Allows user to manually select area of high confidence
    imshowbig(label2rgb(handles.NucChan,'jet',[0,0,0]))
    handles.ConfReg = ManualAreaSelect(1,handles.NucChan,50);
else
    handles.ConfReg=imerode(handles.EmbryoReg,strel('disk',SizeFilt)); % Construct region of confidence by eroding edges of embryo
end

handles.savedata.ConfReg = handles.ConfReg;

imshowbig(label2rgbBackdrop(handles.ConfReg,[1,0,0],[0,1,0],handles.NucChan)) % Display showing high confidence region in red and low confidence in green

NucIndConfReg = nonzeros(unique(handles.LabelNucSeg(handles.ConfReg>0))); % The nuclear indices of nuclei within the selected region

NucIndConRegProp=zeros([length(NucIndConfReg),1]); %Initializing array to store indices

for i=1:length(NucIndConfReg) %Loop to exclude nuclei that have more than some fraction of their area outside of region

    L=mean(handles.ConfReg(handles.RegionPropNuc(NucIndConfReg(i)).PixelIdxList));

    if L>=0.5
        NucIndConRegProp(i)=NucIndConfReg(i);
    end
end

NucIndConRegProp=nonzeros(NucIndConRegProp);

handles.NucIndConfReg = NucIndConfReg;                          % Nuclei in original region of confidence
handles.savedata.NucIndConfReg = NucIndConfReg;

handles.NucIndConRegProp=NucIndConRegProp;                      % Nuclei that have more than some fraction of their area in the region of confidence
handles.savedata.NucIndConRegProp=NucIndConRegProp;

handles.bwConfRegProp=ismember(handles.LabelNucSeg,NucIndConRegProp);  
handles.savedata.bwConfRegProp = handles.bwConfRegProp;

NucIndConRegBnd = nonzeros(unique(handles.LabelNucSeg((imdilate(handles.bwConfRegProp,strel('square',5))-handles.bwConfRegProp)>0)));

handles.NucIndConRegBnd=NucIndConRegBnd;

handles.savedata.NucIndConRegBnd=NucIndConRegBnd;

handles.bwConfRegIncBnd=ismember(handles.LabelNucSeg,[NucIndConRegBnd;NucIndConRegProp]);

handles.savedata.bwConfRegIncBnd=handles.bwConfRegIncBnd;

DISPconf=ismember(handles.LabelNucSeg,NucIndConRegProp)*2+ismember(handles.LabelNucSeg,NucIndConRegBnd)+(handles.LabelNucSeg>0);

imshowbig(label2rgbBackdrop(DISPconf,[0,0,1;0,1,0;1,0,0],[0,0,0], handles.NucChan))

%%%%%% Calculating nuclear density to determin staging
%%%%%%

PixelArea = handles.VoxelsizeX*handles.VoxelsizeX*10^12; % To convert from m^2 tu micro m^2;

Log2NuclearDensity =  NuclearDensityEmbProc(handles.LabelNucSeg,handles.bwConfRegProp,PixelArea,20);

handles.Log2NuclearDensity=Log2NuclearDensity;
handles.savedata.Log2NuclearDensity=Log2NuclearDensity;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% Step e : Segment mRNA1

if handles.stepsorder(step) == 'e';

SaveParm(hObject, eventdata, handles);    
    
SizeFilt = str2double(get(handles.in1,'String'));  %  
RadFilt = str2double(get(handles.in2,'String'));  %  
Thr = str2double(get(handles.in3,'String'));  %  

N = str2double(get(handles.in4,'String'));  %  
smo = str2double(get(handles.in5,'String'));  %  
ThresMin = str2double(get(handles.in6,'String'));  %  
ThresMax = str2double(get(handles.in7,'String'));  %  

handles.mRNA1ChanFlat=handles.mRNA1Chan;

[LabmRNA1New,mRNA1ICenter,threshold] = segmentmRNANoFilt(handles.mRNA1ChanFlat,Thr,N,smo,ThresMin,ThresMax); % Removed filtering step

handles.savedata.mRNA1threshold = threshold;

handles.LabmRNA1New=sparse(LabmRNA1New);
handles.mRNA1ICenter=sparse(mRNA1ICenter);

handles.savedata.LabmRNA1New=sparse(LabmRNA1New);
handles.savedata.mRNA1ICenter=sparse(mRNA1ICenter);

mRNAdispBac=bwmorph(LabmRNA1New>0,'thicken',2);

boundary =  LabelObjectBoundary(mRNAdispBac,1);

imshowbig(bsxfun(@times,label2rgb(handles.mRNA1ChanFlat,'jet',[0,0,0]),uint8(~boundary))+bsxfun(@times,label2rgb(handles.mRNA1ChanFlat,'hsv',[0,0,0]),uint8(boundary)));

title(['The threshold is ' num2str(threshold) ', which gives ' num2str(max(LabmRNA1New(:))) ' objects'])

imshowbig(label2rgb(handles.mRNA1ChanFlat,'jet',[0,0,0]));




end

% Step f : Segment mRNA2

if handles.stepsorder(step) == 'f';
    
SaveParm(hObject, eventdata, handles);
 
SizeFilt = str2double(get(handles.in1,'String'));  %  
RadFilt = str2double(get(handles.in2,'String'));  %  
Thr = str2double(get(handles.in3,'String'));  %  

N = str2double(get(handles.in4,'String'));  %  
smo = str2double(get(handles.in5,'String'));  %  
ThresMin = str2double(get(handles.in6,'String'));  %  
ThresMax = str2double(get(handles.in7,'String'));  %  

handles.mRNA2ChanFlat=handles.mRNA2Chan;

[LabmRNA2New,mRNA2ICenter,threshold] = segmentmRNANoFilt(handles.mRNA1ChanFlat,Thr,N,smo,ThresMin,ThresMax);

handles.savedata.mRNA2threshold = threshold;

handles.LabmRNA2New=sparse(LabmRNA2New);
handles.mRNA2ICenter=sparse(mRNA2ICenter);


handles.savedata.LabmRNA2New=sparse(LabmRNA2New);
handles.savedata.mRNA2ICenter=sparse(mRNA2ICenter);
    

mRNAdispBac=bwmorph(LabmRNA2New>0,'thicken',2);

boundary =  LabelObjectBoundary(mRNAdispBac,1);

imshowbig(bsxfun(@times,label2rgb(handles.mRNA2ChanFlat,'jet',[0,0,0]),uint8(~boundary))+bsxfun(@times,label2rgb(handles.mRNA2ChanFlat,'hsv',[0,0,0]),uint8(boundary)));

title(['The threshold is ' num2str(threshold) ', which gives ' num2str(max(LabmRNA2New(:))) ' objects'])

imshowbig(label2rgb(handles.mRNA2ChanFlat,'jet',[0,0,0]));

    
end

% Step g : Find mRNA pairs

if handles.stepsorder(step) == 'g';

SaveParm(hObject, eventdata, handles);
    
n = str2double(get(handles.in1,'String'));  %  
    
[handles.mRNA1o, handles.mRNA2o, handles.mRNAp] =PairmRNAs(handles.mRNA1ICenter,handles.mRNA2ICenter,n);

handles.mRNA1o=sparse(handles.mRNA1o);

handles.mRNA2o=sparse(handles.mRNA2o);

handles.mRNAp=sparse(handles.mRNAp);

handles.savedata.mRNA1o=handles.mRNA1o;

handles.savedata.mRNA2o=handles.mRNA2o;

handles.savedata.mRNAp=handles.mRNAp;

imshowbig(label2rgb(full(handles.mRNA1o+2*handles.mRNA2o+3*handles.mRNAp),[1,0,0;0,1,0;1,1,0],[0,0,0]))

hold on

 [x,y] = ind2sub(size(handles.NucChan), find(handles.mRNA1o));
 scatter(y,x,50,'r');
 
 [x,y] = ind2sub(size(handles.NucChan), find(handles.mRNA2o));
 scatter(y,x,50,'g');

 [x,y] = ind2sub(size(handles.NucChan), find(handles.mRNAp));
 scatter(y,x,50,'y');
 
 hold off
 

 
 imshowbig(handles.mRNA1Chan)
 hold on
 [x,y] = ind2sub(size(handles.NucChan), find(handles.mRNA2o));
 scatter(y,x,50,'g');
 
 [x,y] = ind2sub(size(handles.NucChan), find(handles.mRNAp));
 scatter(y,x,50,'y');
 
 
 hold off
 
 
 imshowbig(handles.mRNA2Chan)
 hold on
 [x,y] = ind2sub(size(handles.NucChan), find(handles.mRNA1o));
 scatter(y,x,50,'r');
 
 
 [x,y] = ind2sub(size(handles.NucChan), find(handles.mRNAp));
 scatter(y,x,50,'y');
 
 
 hold off
 

 
 
end

% Step h : Disambiguate mRNA

if handles.stepsorder(step) == 'h';
    
    
SaveParm(hObject, eventdata, handles);
    
DilR = str2double(get(handles.in1,'String'));  %  
radarea = str2double(get(handles.in2,'String'));  %  
MaxClusterSize = str2double(get(handles.in3,'String'));  %  
fracsplit = str2double(get(handles.in4,'String'));  %  
genotype = str2double(get(handles.in5,'String'));  %  
SrchDepth = str2double(get(handles.in6,'String'));  %

try
    
[LmRNA,mRNA1oInd, mRNA2oInd, mRNApInd, AllmRNAinNucPre, AllmRNAinNucPost,NucmRNA1oInd,NucmRNA2oInd,NucmRNApInd,genotype]=AmbigiousmRNAsDiff(full(handles.mRNA1o).*handles.bwConfRegIncBnd,full(handles.mRNA2o).*handles.bwConfRegIncBnd,full(handles.mRNAp).*handles.bwConfRegIncBnd,double(handles.LabelNucSeg),handles.NucChan,DilR,radarea,MaxClusterSize,fracsplit,genotype,handles.ACON,SrchDepth);

DisplayDisAmb2Chan(LmRNA,mRNA1oInd, mRNA2oInd, mRNApInd, AllmRNAinNucPre,AllmRNAinNucPost,NucmRNA1oInd,NucmRNA2oInd,NucmRNApInd,handles);

[LmRNA,mRNA1oInd, mRNA2oInd, mRNApInd, AllmRNAinNucPost,NucmRNA1oInd,NucmRNA2oInd,NucmRNApInd] = ExcludingNucleiFromCellsExpressing(LmRNA,mRNA1oInd, mRNA2oInd, mRNApInd, AllmRNAinNucPost,handles.NucIndConRegProp,handles.LabelNucSeg);



RNANucAfter=zeros(1,max(handles.LabelNucSeg(:)))';
RNANucBefore=zeros(1,max(handles.LabelNucSeg(:)))';

for i=1:length(AllmRNAinNucPost);
    
    if ismember(i,mRNApInd)
    RNANucAfter(AllmRNAinNucPost(i))=RNANucAfter(AllmRNAinNucPost(i))+2;
    RNANucBefore(AllmRNAinNucPre(i))=RNANucBefore(AllmRNAinNucPre(i))+2;
    else
    RNANucAfter(AllmRNAinNucPost(i))=RNANucAfter(AllmRNAinNucPost(i))+1;    
    RNANucBefore(AllmRNAinNucPre(i))=RNANucBefore(AllmRNAinNucPre(i))+1;
    end
    
    
end

[X,N]=hist(nonzeros(RNANucAfter),[1:max(RNANucAfter)]);
[Y,NN]=hist(nonzeros(RNANucBefore),[1:max(RNANucBefore)]);


XL=length(X);
YL=length(Y);


if XL<=2*genotype
RNAfracAfterOK = sum(X(1:end));
RNAfracAfterNOOK =0;
else
    
RNAfracAfterOK = sum(X(1:2*genotype));
RNAfracAfterNOOK = sum(X(2*genotype+1:end))/sum(X(1:2*genotype));
end

if YL<=2*genotype
RNAfracBeforeOK = sum(Y(1:end));
RNAfracBeforeNOOK =0;    
else

RNAfracBeforeOK =  sum(Y(1:2*genotype));
RNAfracBeforeNOOK = sum(Y(2*genotype+1:end))/sum(Y(1:2*genotype));
end

DisplayDisAmb2ChanRed(LmRNA,mRNA1oInd, mRNA2oInd, mRNApInd, AllmRNAinNucPost,handles);


handles.LmRNA = LmRNA;
handles.mRNA1oInd = mRNA1oInd;
handles.mRNA2oInd = mRNA2oInd;
handles.mRNApInd =  mRNApInd;
handles.AllmRNAinNucPre = AllmRNAinNucPre;
handles.AllmRNAinNucPost = AllmRNAinNucPost;
handles.NucmRNA1oInd = NucmRNA1oInd;
handles.NucmRNA2oInd = NucmRNA2oInd;
handles.NucmRNApInd = NucmRNApInd;

handles.savedata.LmRNA = sparse(LmRNA);
handles.savedata.mRNA1oInd = mRNA1oInd;
handles.savedata.mRNA2oInd = mRNA2oInd;
handles.savedata.mRNApInd =  mRNApInd;
handles.savedata.AllmRNAinNucPre = AllmRNAinNucPre;
handles.savedata.AllmRNAinNucPost = AllmRNAinNucPost;
handles.savedata.NucmRNA1oInd = NucmRNA1oInd;
handles.savedata.NucmRNA2oInd = NucmRNA2oInd;
handles.savedata.NucmRNApInd = NucmRNApInd;
handles.savedata.genotype=genotype;


catch

[LmRNA1o, numRNA1o] = bwlabeln(full(handles.mRNA1o).*(handles.LabelNucSeg>0));         % labeled matrix with only mRNA1 in region
[LmRNA2o, numRNA2o] = bwlabeln(full(handles.mRNA2o).*(handles.LabelNucSeg>0));         % labeled matrix with only mRNA2 in region
[LmRNAp, numRNAp] = bwlabeln(full(handles.mRNAp).*(handles.LabelNucSeg>0));            % labeled matrix with only mRNA pairs in region


LmRNA2o(find(LmRNA2o))=LmRNA2o(find(LmRNA2o))+numRNA1o;                          % Incresaing labels of mRNA2 for big matrix
LmRNAp(find(LmRNAp))=LmRNAp(find(LmRNAp))+numRNA1o+numRNA2o;                     % Incresaing labels of mRNAp for big matrix
LmRNA=LmRNA1o+LmRNA2o+LmRNAp;                                                    % Label matrix that includes all mRNA types

AllmRNAinNucPre = [LmRNA(find(LmRNA)),handles.LabelNucSeg(find(LmRNA))];                      % Array specifying which mRNAs are assigned to which nuclei before 
AllmRNAinNucPre =sortrows(AllmRNAinNucPre,1);                                    % ambigioous mRNAs are reassigned
AllmRNAinNucPre=AllmRNAinNucPre(:,2);


mRNA1oInd=[1:numRNA1o];                                   %Big matrix lables for mRNA1
mRNA2oInd=[numRNA1o+1:numRNA1o+numRNA2o];                 %Big matrix lables for mRNA2
mRNApInd=[numRNA1o+numRNA2o+1:numRNAp+numRNA1o+numRNA2o]; %Big matrix lables for mRNAp
 
    
    
Nuc=[1:max(handles.LabelNucSeg(:))]';

NucmRNA1oInd=ismember(Nuc,unique(nonzeros(handles.LabelNucSeg(LmRNA1o>0))));
NucmRNA2oInd=ismember(Nuc,unique(nonzeros(handles.LabelNucSeg(LmRNA2o>0))));
NucmRNApInd=ismember(Nuc,unique(nonzeros(handles.LabelNucSeg(LmRNAp>0))));
    
    

end



handles.LmRNA = sparse(LmRNA);
handles.mRNA1oInd = mRNA1oInd;
handles.mRNA2oInd = mRNA2oInd;
handles.mRNApInd =  mRNApInd;

handles.AllmRNAinNucPre = AllmRNAinNucPre;


handles.NucmRNA1oInd = NucmRNA1oInd;
handles.NucmRNA2oInd = NucmRNA2oInd;
handles.NucmRNApInd = NucmRNApInd;




handles.savedata.LmRNA = sparse(LmRNA);
handles.savedata.mRNA1oInd = mRNA1oInd;
handles.savedata.mRNA2oInd = mRNA2oInd;
handles.savedata.mRNApInd =  mRNApInd;

handles.savedata.AllmRNAinNucPre = AllmRNAinNucPre;


handles.savedata.NucmRNA1oInd = NucmRNA1oInd;
handles.savedata.NucmRNA2oInd = NucmRNA2oInd;
handles.savedata.NucmRNApInd = NucmRNApInd;


end

% Step i : Display Nuclei

if handles.stepsorder(step) == 'i';

SaveParm(hObject, eventdata, handles);
Nucp=intersect(find(handles.NucmRNApInd),handles.NucIndConfReg);
Nuc1only=setdiff(intersect(find(handles.NucmRNA1oInd),handles.NucIndConfReg),Nucp);
Nuc2only=setdiff(intersect(find(handles.NucmRNA2oInd),handles.NucIndConfReg),Nucp);

NucOff = setdiff(1:max(handles.LabelNucSeg(:)),[Nuc1only;Nuc2only;Nucp]);

DISPLAY = label2rgbBackdrop((ismember(handles.LabelNucSeg,Nuc1only)+ismember(handles.LabelNucSeg,Nuc2only).*2 + ismember(handles.LabelNucSeg,Nucp).*3 + ismember(handles.LabelNucSeg,NucOff).*4),[1,0,0;0,1,0;1,1,0;0,0,1],[0,0,0],handles.NucChan);

imshowbig(DISPLAY);
handles.DISPLAY=DISPLAY;
handles.savedata.DISPLAY = DISPLAY;

end



% Step s : Save File

if handles.stepsorder(step) == 's';

        handles.DISPLAY(:,:,1) = imadjust(handles.DISPLAY(:,:,1)); % Adjust contrast for display image
        handles.DISPLAY(:,:,2) = imadjust(handles.DISPLAY(:,:,2));
        handles.DISPLAY(:,:,3) = imadjust(handles.DISPLAY(:,:,3));
        
    
        data=handles.savedata;                                     % Take data that was put into savedata

        data.filename=handles.filename;                            % Filename data
        data.datestamp=datestr(now, 'mmmm dd, yyyy HH:MM:SS AM');  % Date of processing 
        
        savename = get(handles.DataFile,'String'); 
        savename = ([savename '_' handles.emb '_EmbProcV3Data']); 
        save([handles.fin '\' savename], 'data');
        
        
        
        
        imwrite(handles.DISPLAY,[handles.fin '\' savename 'image.jpg'],'JPG','Quality',100)
        
         dir = {'Data has been saved'};
        set(handles.directions,'String',dir);
                      
        ExistArray = zeros(10,1);
        for i=1:10
            thename=[handles.fin '\' handles.fname '_' num2str(str2double(handles.emb) + i) '.mat'];
          ExistArray(i) = exist(thename,'file');
        end

          embnum = str2double(handles.emb) + find(ExistArray,1);
          embnum = num2str(embnum);
          handles.emb=embnum;
          set(handles.embin,'String',embnum);
          
          figs = sort(get(0, 'Children'),'descend'); %List of open figures 
          close(figs(2:end))                        %Close all but 
          
          handles = refreshcustom(hObject);

             
end

%catch
%end
%%%%%%AutoCycle %%%%%%%%

 if strcmp(get(handles.AutoCycleButton,'String'),'Auto Cycle On')

 save('trail.mat','hObject','eventdata'); %Way to clear all data from memory so autocycle does not overrun memory
 guidata(hObject,handles);
 clear all
 load('trail.mat');

nextstep_Callback(hObject, eventdata);
run_Callback(hObject, eventdata);
 
 else  
     
 end
 
    
guidata(hObject,handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%========================================================================%
% End of main processing part of code, miscellaneous processing bits below

% --- Function that saves the input parameters 
 
 function  [hObject,handles] = SaveParm(hObject, eventdata, handles)
 
 for i=1:8                                % Saving Inputs
 strut=get(eval(['handles.in' num2str(i)]),'String');
 trr=['handles.savedata.parms.' handles.stepsorder(handles.step) '.in' num2str(i) ' = '''  strut ''';'];
 eval(trr);
 end
 


% --- Executes on button press in AutoCycleButton.
function AutoCycleButton_Callback(hObject, eventdata)

handles=guidata(hObject);

if strcmp(get(handles.AutoCycleButton,'String'),'Auto Cycle Off')
 set(handles.AutoCycleButton,'String','Auto Cycle On')
else
   set(handles.AutoCycleButton,'String','Auto Cycle Off')
end

guidata(hObject,handles);




% --- Executes on button press in SaveDef.
function SaveDef_Callback(hObject, eventdata, handles)
   % record the values of the 6 input boxes for the step now showing
     p1 = (get(handles.in1,'String'));  
     p2 = (get(handles.in2,'String'));  
     p3 = (get(handles.in3,'String'));  
     p4 = (get(handles.in4,'String'));  
     p5 = (get(handles.in5,'String'));  
     p6 = (get(handles.in6,'String'));
     p7 = (get(handles.in7,'String'));
     p8 = (get(handles.in8,'String'));
     
     pars = {p1, p2, p3, p4, p5, p6, p7, p8}; % cell array of strings
  % Export parameters 
     stp_label = get(handles.stepnum,'String'); 
     savelabel = ['EmbProcV3Pars_',handles.stepsorder(handles.step)];  
     % labeled as nucdot_parsi.mat where "i" is the step number 
     save([handles.datafolder, '\', savelabel], 'pars');        % export values


     
%=========================================================================%
%                          Step Controls                                  %      
%=========================================================================%

% --- Executes on button press in nextstep.
function nextstep_Callback(hObject, eventdata)

    handles=guidata(hObject);
    handles.step = handles.step + 1; % forward 1 step
    handles.step = mod(handles.step-1,handles.totstepnum)+1;
    set(handles.stepnum,'String', handles.step); % change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
    setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels



% --- Executes on button press in back.
function back_Callback(hObject, eventdata)
    
    handles=guidata(hObject);
    handles.step = handles.step-1; % go back a step
    handles.step = mod(handles.step-1,handles.totstepnum)+1;
    set(handles.stepnum,'String',handles.step); % Change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles); % update GUI data with new handles
    setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels


% -------------------------------------------------------- %



% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %
%                   File managining scripts                               %  
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %
% This function sets up the new steps with the appropriate input labels and
% defalut label parameters


function setup(hObject,eventdata,handles)   %Setup loads appropriate lables for inputs and default input values

for i=1:8 %load label values for each step as specified in handles.steplabels
  set(eval(['handles.in' num2str(i) 'label']),'String',char(handles.steplabels.(handles.stepsorder(handles.step))(i)))
end

  %load default input values if saved for each step
  if exist( [handles.datafolder '\EmbProcV3Pars_' handles.stepsorder(handles.step) '.mat'],'file')>0   
     load([handles.datafolder '\EmbProcV3Pars_' handles.stepsorder(handles.step) '.mat'])
   
        for i=1:8
        set(eval(['handles.in' num2str(i)]),'String',pars{i})
        end
  else
       dir = {['Parameter file does not exist! Save one!!']};
        
  end


set(handles.directions,'String',handles.steplabels.dir.(handles.stepsorder(handles.step))); 


% set directions for each step as specified in handles.steplabels

guidata(hObject, handles); % update GUI data with new labels


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %



% Automatically return the program to step 1 if the image source directory,
% file name, or image number are changed.  

function froot_Callback(hObject, eventdata, handles)

    handles = refreshcustom(hObject);
    handles.step = 1;  % starting step is step 1 
    set(handles.stepnum,'String',handles.step); % change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
    setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels




function embin_Callback(hObject, eventdata, handles)

    
    handles = refreshcustom(hObject);    
    handles.step = 1;  % starting step is step 1 
    set(handles.stepnum,'String',handles.step); % change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
     setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels



function source_Callback(hObject, eventdata, handles)
    
    handles = refreshcustom(hObject);
    
    handles.step = 1;  % starting step is step 1 
    set(handles.stepnum,'String',handles.step); % change step label in GUI
    handles.output = hObject; % update handles object with new step number
    guidata(hObject, handles);  % update GUI data with new handles
    setup(hObject, eventdata, handles); % set up labels and default values for new step
    guidata(hObject, handles); % update GUI data with new labels
  

function handles = refreshcustom(hObject) % refresh such that only key info is kept in handles

handles = guidata(hObject);
handlesNew = guihandles(hObject);
handlesNew.step=handles.step;
handlesNew.steplabels=handles.steplabels;
handlesNew.stepsorder = handles.stepsorder;
handlesNew.datafolder = handles.datafolder;
handlesNew.totstepnum = handles.totstepnum;

handles=handlesNew;


% Open file browser to select source folder 
function SourceBrowse_Callback(hObject, eventdata, handles)
 sourcefile = uigetdir; % prompts user to select directory
  set(handles.source,'String',sourcefile);

% --- Outputs from this function are returned to the command line.
% Important!!!!
function varargout = EmbProcV3_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function embin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function source_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function froot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in1_Callback(hObject, eventdata, handles)
function in1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function in2_Callback(hObject, eventdata, handles)
function in2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function in3_Callback(hObject, eventdata, handles)

function in3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in4_Callback(hObject, eventdata, handles)
    
function in4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in5_Callback(hObject, eventdata, handles)
    
function in5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in6_Callback(hObject, eventdata, handles)
    
function in6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in7_Callback(hObject, eventdata, handles)
    
function in7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in8_Callback(hObject, eventdata, handles)
    
function in8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function DataFile_Callback(hObject, eventdata, handles)
    
% --- Executes during object creation, after setting all properties.
function DataFile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
