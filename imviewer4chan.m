%                               imviewer4chan.m
%
% Jacques Bothma                                    Last Modified: 09/13/11     
% Levine Lab, UC Berkeley                       
% Functionally complete                          
%
%
%% Attribution:
%  Feel free to use, modify and distribute this code provided that you
%  attribute Jacques Bothma and Alistair Boettiger for development.
%  This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
%  To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.
%
%  Important Notes:
%  This version written for Windows.
%
%  Overview:
%
%  This GUI is designed to take a stack of fluorescence confocal images of a Drosophila embryo
%  that has been stained for DNA and up to three different types of mRNA and project
%  this onto a 2D image that is processed by the EmbProc GUI to segment the
%  nuclei and mRNA. This image can be saved as a tiff or structure in a mat file.
%  I prefer saving the data as a structure in a mat file. It allows you to
%  add on arbitrary info easily and it saves space, the tiffs take up about
%  twice the space of the equivalent mat file. It also saves a lower res
%  jpg image and a lowe res jpg image of the first and last slices so one
%  can readily make sure the zstack included all releevtn regions of the
%  embryo.
%
%  We use a Zeiss LSM 700 confocal microscope and so the
%  associated scripts are designed to read the info in from the .lsm file format
%  saved by the Zeiss software Zen. Lsm files are essentially just big
%  tiff files with multiple subfiles where the microscope info is included in the
%  header of the first file. Extensive documentation exists and can be
%  obtained by constacting me (jpbothma through gmail) directly.
%
%  The scripts for reading from lsm were modified versions of code written
%  by Peter Hawley Li.
%
%  This GUI is a step based approach where all successive filtering and
%  segmentation steps are grouped into individual steps that allow for user
%  input. There is also a an autocycle mode in which the standard
%  parameters are applied blindly. Some parameter values are stored in
%  external files while most of the image data is stored in the hadles structure of
%  the gui figure.
% 
%  Load parameters for first run.
%


function varargout = imviewer4chan(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imviewer4chan_OpeningFcn, ...
                   'gui_OutputFcn',  @imviewer4chan_OutputFcn, ...
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





% --- Executes just before imviewer4chan is made visible.
function imviewer4chan_OpeningFcn(hObject, eventdata, handles, varargin)
handles.step = 1; 
handles.datafolder = 'C:\Users\bothma\Documents\MATLAB\Code\ImageProcessing\ParameterFiles';
setup(hObject,eventdata,handles); 
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = imviewer4chan_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Command that sets up the reading in of typed commands, myFunction is
% the function that ties the input to the output.
%set(handles.figure1,'KeyPressFcn',{@myFunction,hObject, eventdata})


 try %tries to load the name of the last processed image if it exists. this gets saved in step 1. 

     load([handles.datafolder, '\LastNamesimviewer4chan.mat']); 
     
     set(handles.foutput,'String',namestore.foutput);
     set(handles.fdir,'String',namestore.fdir);
     set(handles.froot,'String',namestore.froot);
     set(handles.outlabel,'String',namestore.outlabel);
 
 catch
 end



% --- Executes on button press in run.
function [hObject,handles] = run_Callback(hObject, eventdata, handles);
    
step = handles.step;


if step == 1;
     [hObject,handles] = step1(hObject, eventdata, handles);
end

if step == 2;
     [hObject,handles] = step2(hObject, eventdata, handles);
end

if step == 3;
     [hObject,handles] = step3(hObject, eventdata, handles);
end

handles.output = hObject;
guidata(hObject, handles); 

%=========================================================================%
% STEP 1: PROJECT SELECTED RAW DATA FILES

function [hObject,handles] = step1(hObject, eventdata, handles) 
    % Input information
    fin = get(handles.fdir,'String');  % folder containing raw image data 
    handles.name =  get(handles.froot,'String'); %  file root 
    handles.imnum = get(handles.imnumber,'String'); %  image number as string

    
    % Order to save channels in
    
    handles.chl =[get(handles.ch1Order,'String');...
                  get(handles.ch2Order,'String');...
                  get(handles.ch3Order,'String');...
                  get(handles.ch4Order,'String')];

    
    handles.chlnum = str2num(char(handles.chl));
              
              
    % RGB coordinates of the color to display different channels
   
    C(1,1) = str2double(get(handles.ch1R,'String'));
    C(1,2) = str2double(get(handles.ch1G,'String'));
    C(1,3) = str2double(get(handles.ch1B,'String'));
    C(2,1) = str2double(get(handles.ch2R,'String'));
    C(2,2) = str2double(get(handles.ch2G,'String'));
    C(2,3) = str2double(get(handles.ch2B,'String'));
    C(3,1) = str2double(get(handles.ch3R,'String'));
    C(3,2) = str2double(get(handles.ch3G,'String'));
    C(3,3) = str2double(get(handles.ch3B,'String'));
    C(4,1) = str2double(get(handles.ch4R,'String'));
    C(4,2) = str2double(get(handles.ch4G,'String'));
    C(4,3) = str2double(get(handles.ch4B,'String'));

    handles.C=C;

    
    
    
    % Choose which slices to project
    
    ch1subset = str2num(get(handles.in1,'String'));
    ch2subset = str2num(get(handles.in2,'String'));
    ch3subset = str2num(get(handles.in3,'String'));
    ch4subset = str2num(get(handles.in4,'String'));

    nucchan = str2num(get(handles.in5,'String')); % Specify which channel is the nuclear channel
    
    maxProj = str2num(cell2mat(get(handles.in6,'String'))); %  Which channels to max project
    
    LogFilterMaxProj = str2num(cell2mat(get(handles.in7,'String'))); % Which channels to filter
    
    FilterSize = str2num(cell2mat(get(handles.in8,'String')));

    FilterRadius = str2num(cell2mat(get(handles.in9,'String')));

    FracSat= str2num(cell2mat(get(handles.in10,'String')));
    
    ObjT = str2num(cell2mat(get(handles.in11,'String')));
    
    FilterParm=[FilterSize,FilterRadius,FracSat];
    
    subsets=[ch1subset;ch2subset;ch3subset;ch4subset];

    % This section parses the lsm file. It takes the info like where in the
    % lsm file that actual images are stored and saves that to a mat file.
    % The info in this mat file is then used by the loadlsm file to load
    % image data from the lsm file directly
    
    if exist([fin '\' handles.name '.mat'])~=2
        disp(['Parsing ' fin '\' handles.name '.lsm'])
        parselsm([fin '\' handles.name '.lsm'])
    else
    end
    
    %Load lsm file info
    load([fin '\' handles.name]);
    
    Datas.LSM_info.DataSummary.Fluorophores
    handles.VoxelSizeX = Datas.LSM_info.VoxelSizeX;
    handles.VoxelSizeY = Datas.LSM_info.VoxelSizeY;
    handles.VoxelSizeZ = Datas.LSM_info.VoxelSizeZ;
    handles.Fluorophores = Datas.LSM_info.DataSummary.Fluorophores;
    handles.LSM_info=Datas.LSM_info;
    
    I = fxnprojectLSM4Chan(fin,handles.name,handles.chlnum,subsets,maxProj,LogFilterMaxProj,FilterParm,nucchan,ObjT,str2num(handles.imnum)); % project image from raw data
    
    Io=I.Io;

    handles.I=I;
    handles.output = hObject;
    guidata(hObject, handles); 
    

    [a,b,c]=size(Io);
    
    Dis=zeros([a,b,3],class(Io));
    
    imshowbig(Io(:,:,1))
    imshowbig(Io(:,:,2))
    imshowbig(Io(:,:,3))
    
    chlnumfind = find(handles.chlnum);
    
    for i=1:length(chlnumfind)
    
   Dis = label2rgbBackdrop([],handles.C(chlnumfind(i),:),[0,0,0],Io(:,:,i)) +Dis;

    end
    imshowbig(Dis)
    
    
    [hObject,handles] = SaveParm(hObject, eventdata, handles); % Function below that saves input parameters to handles.
   
    

%%%%% Saves the name of the files currently being processed to reload the
%%%%% next time the gui is being run.

namestore=[];
namestore.foutput=get(handles.foutput,'String');
namestore.fdir=get(handles.fdir,'String');
namestore.froot = get(handles.froot,'String'); 
namestore.outlabel=get(handles.outlabel,'String'); 

save([handles.datafolder, '\LastNamesimviewer4chan.mat'],'namestore')

%%%%%%%%%%%%%


    
%%%%% End Step 1 %%%%%%%%%%%%%%%%%%
    
% STEP 2: CLEAN AND ORIENT IMAGE

    function [hObject,handles] = step2(hObject,eventdata,handles)

        % Image cleaning parameters
    emb_chn = str2double(get(handles.in1,'String'));
    ManualOrientation  = str2double(get(handles.in2,'String'));
     
    I = fxnclean4Chan(handles.I.Io,emb_chn,handles.I.bw,handles.I.NucBright,ManualOrientation); % clean and orient image

    [a,b,c]=size(I);
    
    Dis=zeros([a,b,3],class(I));
 
    chlnumfind = find(handles.chlnum);
    
    for i=1:length(chlnumfind)
    
   Dis = label2rgbBackdrop([],handles.C(chlnumfind(i),:),[0,0,0],I(:,:,i)) +Dis;

    end
    imshowbig(Dis)
    
    handles.If=I;

    [hObject,handles] = SaveParm(hObject, eventdata, handles); % Function below that saves input parameters to handles.
        
    guidata(hObject, handles); 
    handles.output = hObject;

%%%%% End Step 2 %%%%%%%%%%%%%%%%%%    

    
 % STEP 3: SAVE DATA
   function [hObject,handles] = step3(hObject,eventdata,handles)
     
       outname = get(handles.outlabel,'String');
      fout = get(handles.foutput,'String');
      
      SaveType= get(handles.in1,'String');
      strcmp(SaveType,'mat');

      
   [hObject,handles] = SaveParm(hObject, eventdata, handles); % Function below that saves input parameters to handles.
      
      
      if exist(fout)==7 %Checks to see if output folder exists and if it doesn't it creates it.
      else
          mkdir(fout);
      end
      
         handles.chlnum =str2num(char([get(handles.ch1Order,'String');...
                  get(handles.ch2Order,'String');...
                  get(handles.ch3Order,'String');...
                  get(handles.ch4Order,'String')]));
              
      

       temp = zeros(size(handles.If),class(handles.If));
      
       
         chlnumfind = find(handles.chlnum);
       
       for i=1:length(chlnumfind)
           
       temp(:,:,handles.chlnum(chlnumfind(i)))=handles.If(:,:,i);
       
       end

       [xx,yy,zz] = size(temp);
       
       [x,y,z] = size(handles.I.IoF);
       
    IoFDis = zeros([x,y,3],class(handles.I.IoF));
    IoLDis = zeros([x,y,3],class(handles.I.IoL));
    Display= zeros([xx,yy,3],class(temp));
   
    for i=1:length(chlnumfind)
    
    IoFDis = label2rgbBackdrop([],handles.C(chlnumfind(i),:),[0,0,0],handles.I.IoF(:,:,i)) +IoFDis;
    IoLDis = label2rgbBackdrop([],handles.C(chlnumfind(i),:),[0,0,0],handles.I.IoL(:,:,i)) +IoLDis;
    Display = label2rgbBackdrop([],handles.C(chlnumfind(i),:),[0,0,0],handles.If(:,:,i)) +Display;
   
    end

    
    
    
      fluorophores='';
      
      for jj=1:length(chlnumfind)
            fluorophores=[fluorophores, 'Ch' num2str(jj) ' : ' handles.Fluorophores{find(handles.chlnum==jj)} ';'];
      end
      
    
      if strcmp(SaveType,'tif')
   


      figurecomment=['Date: ' datestr(now, 'mmmm dd, yyyy HH:MM:SS AM') '; Filename :' handles.name ...
        '; VoxelsizeX (m),VoxelSizeY (m):' num2str(handles.VoxelSizeX) ',' num2str(handles.VoxelSizeY) ';'];
      
      imwrite(temp,[fout,'\', outname,'_'...
               handles.imnum,'.tif'],'tif','Description',figurecomment);
      
      elseif strcmp(SaveType,'mat')
      
      Data.Image=temp;
      Data.Date = datestr(now, 'mmmm dd, yyyy HH:MM:SS AM');
      Data.Filename = handles.name;
      Data.fluorophores = fluorophores;
      Data.VoxelsizeX = handles.VoxelSizeX;
      Data.VoxelsizeY = handles.VoxelSizeY;
      Data.VoxelsizeZ = handles.VoxelSizeZ;
      Data.NumberSlices = handles.LSM_info.DimensionZ;
      Data.DataSummary=handles.LSM_info.DataSummary;
      Data.SavedParameters=handles.savedata;
      
      save([fout,'\', outname,'_' handles.imnum],'Data')
      
      else
          
          disp('Error, Need to save as either tif or mat')
      
      end
           

           
      imwrite(IoFDis,[fout,'\', outname,'_'...
               handles.imnum,'_First','.jpg'],'jpg');

        imwrite(IoLDis,[fout,'\', outname,'_'...
               handles.imnum,'_Last','.jpg'],'jpg');         

        
        imwrite(Display,[fout,'\', outname,'_'...
               handles.imnum,'_Preview','.jpg'],'jpg');    
           

               
    guidata(hObject, handles); 
    handles.output = hObject;

           
          imnum = str2double(get(handles.imnumber,'String'));
          imnum = imnum + 1; % advance to next image

          handles.imnum = imnum;
          set(handles.imnumber,'String',imnum);
          

%%%%% End Step 3 %%%%%%%%%%%%%%%%%%               
           
%=========================================================================%

function setup(hObject,eventdata,handles);
if handles.step == 1; 
    
        load([handles.datafolder '\imviewer_1.mat']);
    
        set(handles.in1label,'String','Start/End Ch1');
        set(handles.in1,'String', pars{1});
        set(handles.in2label,'String','Start/End Ch2');
        set(handles.in2,'String', pars{2});
        set(handles.in3label,'String','Start/End Ch3');
        set(handles.in3,'String',pars{3});
        set(handles.in4label,'String','Start/End Ch4');
        set(handles.in4,'String',pars{4});

        set(handles.in5label,'String','Nuclei channel'); 
        set(handles.in5,'String',pars{5});

        set(handles.in6label,'String','Chan to MaxP'); 
        set(handles.in6,'String',pars{6});
        
        set(handles.in7label,'String','Chan to Filt+MaxP'); 
        set(handles.in7,'String',pars{7});

        set(handles.in8label,'String','Log Filter Size'); 
        set(handles.in8,'String',pars{8});

        set(handles.in9label,'String','Log Filter Width'); 
        set(handles.in9,'String',pars{9});
        
        set(handles.in10label,'String','Frac Max'); 
        set(handles.in10,'String',pars{10});
        
             
        set(handles.in11label,'String','Nuc Threshold'); 
        set(handles.in11,'String',pars{11});
        
end

if handles.step == 2;
     load([handles.datafolder '\imviewer_2.mat']);
     
        set(handles.in1label,'String','Object Channel'); 
        set(handles.in1,'String', pars{1} );
        
        set(handles.in2label,'String','Manual Orientation'); 
        set(handles.in2,'String', pars{2} );
        
        for i=3:11
              set(handles.(['in' num2str(i) 'label']),'String',''); 
              set(handles.(['in' num2str(i)]),'String',''); 
        end
            
end


if handles.step == 3;
     load([handles.datafolder '\imviewer_3.mat']);
     
        set(handles.in1label,'String','Save Type'); 
        set(handles.in1,'String', pars{1});


        for i=2:11
              set(handles.(['in' num2str(i) 'label']),'String',''); 
              set(handles.(['in' num2str(i)]),'String',''); 
        end
            
end

load([handles.datafolder '\imviewercolor.mat'])

set(handles.ch1R,'String', parschan{1});
set(handles.ch1G,'String', parschan{2});        
set(handles.ch1B,'String', parschan{3});

set(handles.ch2R,'String', parschan{4});
set(handles.ch2G,'String', parschan{5});
set(handles.ch2B,'String', parschan{6});

set(handles.ch3R,'String', parschan{7});
set(handles.ch3G,'String', parschan{8});
set(handles.ch3B,'String', parschan{9});

set(handles.ch4R,'String', parschan{10});
set(handles.ch4G,'String', parschan{11});
set(handles.ch4B,'String', parschan{12});



Chcell=cell(4,1);
for i=1:4
Chcell{i} = parschan{12+i};
end


set(handles.ch1Order,'String', Chcell{1});
set(handles.ch2Order,'String', Chcell{2});
set(handles.ch3Order,'String', Chcell{3});
set(handles.ch4Order,'String', Chcell{4});



% --- Executes on button press in Next.
function [hObject,handles] = Next_Callback(hObject, eventdata, handles);
    
    
    handles.step = handles.step + 1;
    handles.step = mod(handles.step-1,3)+1;
    
    
    set(handles.stepnum,'String',handles.step);    
    
    handles.output = hObject; 
   
    guidata(hObject, handles);  
    
    setup(hObject, eventdata, handles); 
    
    guidata(hObject, handles);
    

 
% --- Executes on button press in back.
function [hObject,handles] = back_Callback(hObject, eventdata, handles);
    
handles.step = handles.step - 1;
   handles.step = mod(handles.step-1,3)+1;
   
 set(handles.stepnum,'String',handles.step);
    handles.output = hObject; 
    guidata(hObject, handles);  
    setup(hObject, eventdata, handles); 
    guidata(hObject, handles); 
 
    

% ---- Update colors
function [hObject,handles] = cupdate_Callback(hObject, eventdata, handles)
    
    
C(1,1) = str2double(get(handles.ch1R,'String'));
C(1,2) = str2double(get(handles.ch1G,'String'));
C(1,3) = str2double(get(handles.ch1B,'String'));
C(2,1) = str2double(get(handles.ch2R,'String'));
C(2,2) = str2double(get(handles.ch2G,'String'));
C(2,3) = str2double(get(handles.ch2B,'String'));
C(3,1) = str2double(get(handles.ch3R,'String'));
C(3,2) = str2double(get(handles.ch3G,'String'));
C(3,3) = str2double(get(handles.ch3B,'String'));
C(4,1) = str2double(get(handles.ch4R,'String'));
C(4,2) = str2double(get(handles.ch4G,'String'));
C(4,3) = str2double(get(handles.ch4B,'String'));

handles.C=C;

if isfield(handles, 'I')
I = handles.I;


    [a,b,c]=size(I)
    Dis=zeros([a,b,3],class(I));
    
    for i=1:length(handles.chlnum)
    Dis=label2rgbBackdrop([],handles.C(i,:),[0,0,0],I(:,:,i))+Dis;
    end
    imshowbig(Dis)
    
else
I = handles.Io;


    [a,b,c]=size(I)
    Dis=zeros([a,b,3],class(I));
    
    for i=1:length(handles.chlnum)
    Dis=label2rgbBackdrop([],handles.C(i,:),[0,0,0],I(:,:,i))+Dis;
    end
    imshowbig(Dis)
    
    
end




% --- Autocycle
function cycle_Callback(hObject, eventdata, handles)
    
  h = 1; 
   while h==1;
         
 
          [hObject,handles] = step1(hObject, eventdata, handles);
          [hObject,handles] = Next_Callback(hObject, eventdata, handles);
          [hObject,handles] = step2(hObject, eventdata, handles);
          [hObject,handles] = Next_Callback(hObject, eventdata, handles);
          [hObject,handles] = step3(hObject, eventdata, handles);

          
          imnum = str2double(get(handles.imnumber,'String'));
          handles.imnum = imnum;
          set(handles.imnumber,'String',imnum); 
          
          % restart at step 1
            handles.step = 1;  
            set(handles.stepnum,'String',handles.step);
            setup(hObject, eventdata, handles); 
            handles.output = hObject; 
            guidata(hObject, handles);
          
          disp(['cycle: ',handles.imnum]);

          figs = sort(get(0, 'Children'),'descend'); %List of open figures 
    
          close(figs(2:end))

  end


% --- update step default parameters 

function [hObject,handles] = sdefaults_Callback(hObject, eventdata, handles);
    
    
    
stp_label = get(handles.stepnum,'String');
savelabel = ['imviewer_',stp_label];

pars = {}; 


for i=1:15
    pars{i} = get(handles.(['in' num2str(i)]),'String');
end

save([handles.datafolder, '\' savelabel],'pars'); 

parschan = {};

 parschan{1} = get(handles.ch1R,'String');
 parschan{2} = get(handles.ch1G,'String');
 parschan{3} = get(handles.ch1B,'String');
 parschan{4} = get(handles.ch2R,'String');
 parschan{5} = get(handles.ch2G,'String');
 parschan{6} = get(handles.ch2B,'String');
 parschan{7} = get(handles.ch3R,'String');
 parschan{8} = get(handles.ch3G,'String');

 parschan{9} = get(handles.ch3B,'String');
 parschan{10} = get(handles.ch4R,'String');
 parschan{11} = get(handles.ch4G,'String');
 parschan{12} = get(handles.ch4B,'String');
 
 parschan{13} = get(handles.ch1Order,'String');
 parschan{14} = get(handles.ch2Order,'String');
 parschan{15} = get(handles.ch3Order,'String');
 parschan{16} = get(handles.ch4Order,'String');
  
save([handles.datafolder, '\imviewercolor.mat'],'parschan')




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%%%% Save parameter script %%%%%%%%%%

 function  [hObject,handles] = SaveParm(hObject, eventdata, handles)
 
 for i=1:15                                % Saving Inputs
 par=get(eval(['handles.in' num2str(i)]),'String');
 lab = get(eval(['handles.in' num2str(i) 'label']),'String');
 
 handles.savedata.(['step' num2str(handles.step)]).(['in' num2str(i)]) = par;
 handles.savedata.(['step' num2str(handles.step)]).(['in' num2str(i) 'label']) = lab;
 
 end
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








function imnumber_Callback(hObject, eventdata, handles)
    
handles.step =  1;
 set(handles.stepnum,'String',handles.step);
    handles.output = hObject; 
    guidata(hObject, handles);  
    setup(hObject, eventdata, handles); 
    guidata(hObject, handles);
  
function froot_Callback(hObject, eventdata, handles);
    
    froot = get(handles.froot,'String');
    guidata(hObject, handles);
        
        
        
function froot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end
    
function imnumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function fdir_Callback(hObject, eventdata, handles)
    
function fdir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end




function in1_Callback(hObject, eventdata, handles)


function in1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function in2_Callback(hObject, eventdata, handles)

function in2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end


function in3_Callback(hObject, eventdata, handles)

function in3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function in4_Callback(hObject, eventdata, handles)

function in4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end


function redname_Callback(hObject, eventdata, handles)

function redname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function greenname_Callback(hObject, eventdata, handles)

function greenname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function bluename_Callback(hObject, eventdata, handles)

function bluename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function ch1R_Callback(hObject, eventdata, handles)

function ch1R_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch1G_Callback(hObject, eventdata, handles)

function ch1G_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch1B_Callback(hObject, eventdata, handles)

function ch1B_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch3R_Callback(hObject, eventdata, handles)

function ch3R_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch3G_Callback(hObject, eventdata, handles)

function ch3G_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch3B_Callback(hObject, eventdata, handles)

function ch3B_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ch2R_Callback(hObject, eventdata, handles)

function ch2R_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch2G_Callback(hObject, eventdata, handles)

function ch2G_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ch2B_Callback(hObject, eventdata, handles)

function ch2B_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in5_Callback(hObject, eventdata, handles)

function in5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function foutput_Callback(hObject, eventdata, handles)

function foutput_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outlabel_Callback(hObject, eventdata, handles)
    


% --- Executes during object creation, after setting all properties.
function outlabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function othername_Callback(hObject, eventdata, handles)
% hObject    handle to othername (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of othername as text
%        str2double(get(hObject,'String')) returns contents of othername as a double


% --- Executes during object creation, after setting all properties.
function othername_CreateFcn(hObject, eventdata, handles)
% hObject    handle to othername (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch4R_Callback(hObject, eventdata, handles)
% hObject    handle to ch4R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch4R as text
%        str2double(get(hObject,'String')) returns contents of ch4R as a double


% --- Executes during object creation, after setting all properties.
function ch4R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch4R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch4G_Callback(hObject, eventdata, handles)
% hObject    handle to ch4G (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch4G as text
%        str2double(get(hObject,'String')) returns contents of ch4G as a double


% --- Executes during object creation, after setting all properties.
function ch4G_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch4G (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch4B_Callback(hObject, eventdata, handles)
% hObject    handle to ch4B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch4B as text
%        str2double(get(hObject,'String')) returns contents of ch4B as a double


% --- Executes during object creation, after setting all properties.
function ch4B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch4B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch1Order_Callback(hObject, eventdata, handles)
% hObject    handle to ch1Order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch1Order as text
%        str2double(get(hObject,'String')) returns contents of ch1Order as a double


% --- Executes during object creation, after setting all properties.
function ch1Order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1Order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch2order_Callback(hObject, eventdata, handles)
% hObject    handle to ch2order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch2order as text
%        str2double(get(hObject,'String')) returns contents of ch2order as a double


% --- Executes during object creation, after setting all properties.
function ch2order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch2order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch3Order_Callback(hObject, eventdata, handles)
% hObject    handle to ch3Order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch3Order as text
%        str2double(get(hObject,'String')) returns contents of ch3Order as a double


% --- Executes during object creation, after setting all properties.
function ch3Order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3Order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch4Order_Callback(hObject, eventdata, handles)
% hObject    handle to ch4Order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch4Order as text
%        str2double(get(hObject,'String')) returns contents of ch4Order as a double


% --- Executes during object creation, after setting all properties.
function ch4Order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch4Order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch2Order_Callback(hObject, eventdata, handles)
% hObject    handle to ch2Order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch2Order as text
%        str2double(get(hObject,'String')) returns contents of ch2Order as a double


% --- Executes during object creation, after setting all properties.
function ch2Order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch2Order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in10_Callback(hObject, eventdata, handles)
% hObject    handle to in10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in10 as text
%        str2double(get(hObject,'String')) returns contents of in10 as a double


% --- Executes during object creation, after setting all properties.
function in10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in9_Callback(hObject, eventdata, handles)
% hObject    handle to in9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in9 as text
%        str2double(get(hObject,'String')) returns contents of in9 as a double


% --- Executes during object creation, after setting all properties.
function in9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in8_Callback(hObject, eventdata, handles)
% hObject    handle to in8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in8 as text
%        str2double(get(hObject,'String')) returns contents of in8 as a double


% --- Executes during object creation, after setting all properties.
function in8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in7_Callback(hObject, eventdata, handles)
% hObject    handle to in7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in7 as text
%        str2double(get(hObject,'String')) returns contents of in7 as a double


% --- Executes during object creation, after setting all properties.
function in7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in6_Callback(hObject, eventdata, handles)
% hObject    handle to in6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in6 as text
%        str2double(get(hObject,'String')) returns contents of in6 as a double


% --- Executes during object creation, after setting all properties.
function in6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in15_Callback(hObject, eventdata, handles)
% hObject    handle to in15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in15 as text
%        str2double(get(hObject,'String')) returns contents of in15 as a double


% --- Executes during object creation, after setting all properties.
function in15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in14_Callback(hObject, eventdata, handles)
% hObject    handle to in14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in14 as text
%        str2double(get(hObject,'String')) returns contents of in14 as a double


% --- Executes during object creation, after setting all properties.
function in14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in13_Callback(hObject, eventdata, handles)
% hObject    handle to in13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in13 as text
%        str2double(get(hObject,'String')) returns contents of in13 as a double


% --- Executes during object creation, after setting all properties.
function in13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in12_Callback(hObject, eventdata, handles)
% hObject    handle to in12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in12 as text
%        str2double(get(hObject,'String')) returns contents of in12 as a double


% --- Executes during object creation, after setting all properties.
function in12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function in11_Callback(hObject, eventdata, handles)
% hObject    handle to in11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of in11 as text
%        str2double(get(hObject,'String')) returns contents of in11 as a double


% --- Executes during object creation, after setting all properties.
function in11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
