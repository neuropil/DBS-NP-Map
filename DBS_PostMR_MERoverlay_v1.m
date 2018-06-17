function varargout = DBS_PostMR_MERoverlay_v1(varargin)
% DBS_POSTMR_MEROVERLAY_V1 MATLAB code for DBS_PostMR_MERoverlay_v1.fig
%      DBS_POSTMR_MEROVERLAY_V1, by itself, creates a new DBS_POSTMR_MEROVERLAY_V1 or raises the existing
%      singleton*.
%
%      H = DBS_POSTMR_MEROVERLAY_V1 returns the handle to a new DBS_POSTMR_MEROVERLAY_V1 or the handle to
%      the existing singleton*.
%
%      DBS_POSTMR_MEROVERLAY_V1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DBS_POSTMR_MEROVERLAY_V1.M with the given input arguments.
%
%      DBS_POSTMR_MEROVERLAY_V1('Property','Value',...) creates a new DBS_POSTMR_MEROVERLAY_V1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DBS_PostMR_MERoverlay_v1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DBS_PostMR_MERoverlay_v1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DBS_PostMR_MERoverlay_v1

% Last Modified by GUIDE v2.5 16-Jun-2018 09:04:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DBS_PostMR_MERoverlay_v1_OpeningFcn, ...
    'gui_OutputFcn',  @DBS_PostMR_MERoverlay_v1_OutputFcn, ...
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

% TO DO
% ADD - INPUTs to RUN and TEST

% --- Executes just before DBS_PostMR_MERoverlay_v1 is made visible.
function DBS_PostMR_MERoverlay_v1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DBS_PostMR_MERoverlay_v1 (see VARARGIN)

grid off
set(handles.showMR,'XTickLabel',[])
set(handles.showMR,'YTickLabel',[])
set(handles.showMR,'ZTickLabel',[])
set(handles.showMR,'XTick',[])
set(handles.showMR,'YTick',[])
set(handles.showMR,'ZTick',[])

handles.colorB.Enable = 'off';
handles.sizeB.Enable = 'off';

% Starting off

% Step 1
handles.step1.ForegroundColor = [0.5 0.5 0.5];

% Step 2
handles.step2.ForegroundColor = [0.5 0.5 0.5];

% Step 3
handles.step3.ForegroundColor = [0.5 0.5 0.5];

% Step 4
handles.step4.ForegroundColor = [0.5 0.5 0.5];

% Step 5
handles.step5.ForegroundColor = [0.5 0.5 0.5];

% Step 6
handles.step6.ForegroundColor = [0.5 0.5 0.5];

% Step 7
handles.step7.ForegroundColor = [0.5 0.5 0.5];

% Step 8
handles.step8.ForegroundColor = [0.5 0.5 0.5];
handles.mni_button.Enable = 'off';
handles.native_button.Enable = 'off';

handles.definput = {'s','f','0','1','1.3'};
handles.mri2use = 'native';


% handles.run.Enable = 'off';

% Choose default command line output for DBS_PostMR_MERoverlay_v1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DBS_PostMR_MERoverlay_v1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DBS_PostMR_MERoverlay_v1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.step2.ForegroundColor,  [0 0 0]) && isequal(handles.step6.ForegroundColor ,[0 0 0])
    
    ele_n_nii = [handles.EN_path , handles.EN_Mask];
    ele_m_nii = [handles.EM_path , handles.EM_Mask];
    
    mri_n_nii = [handles.PN_mripath , handles.PN_MRI];
    mri_m_nii = [handles.PM_mripath , handles.PM_MRI];
    
    mni_nii = [handles.mniPath , handles.mniMRI];
    neuroDAT = [handles.ephysPath , handles.ephysCSV]; 
    
    % ele_nii  , solidType , actELEDIM , numPts
    
    if strcmp(handles.mri2use,'mni')
        [ handles.output_args ] = ExtractDBSLeadPoly_02( ele_m_nii  , str2double(handles.definput{4}) , str2double(handles.definput{5}) , 80);
        handles.mriMatDat = load_nii(mni_nii);
    elseif strcmp(handles.mri2use,'native')
        [ handles.output_args ] = ExtractDBSLeadPoly_02( ele_n_nii  , str2double(handles.definput{4}) , str2double(handles.definput{5}) , 80);
        handles.mriMatDat = load_nii(mri_n_nii);
    end
    
    [handles.brainIm] = Process_MRI_01_GUI(handles.mriMatDat);
    [handles.Xsl , handles.Ysl , handles.Zsl ] = size(handles.brainIm);
    
    if strcmp(handles.mri2use,'mni')
    
    handles.xSliceN = round(handles.Xsl/2) + 15;
    handles.ySliceN = round(handles.Ysl/2) - 15;
    handles.zSliceN = round(handles.Zsl/2) - 15;
    
    else
            handles.xSliceN = round(handles.Xsl/2);
    handles.ySliceN = round(handles.Ysl/2);
    handles.zSliceN = round(handles.Zsl/2);
    end
    handles.neuroDat = readtable(neuroDAT);
    handles.color = 1;
    handles.size = 0;
    
    [handles.x,handles.y,handles.z,handles.sz,handles.col] =...
        getNeuroGUIparams(handles.neuroDat, handles.output_args,handles.color,handles.size);
    
%     set(handles.xAxVal,'Min',1);
%     set(handles.xAxVal,'Max',handles.Xsl);
%     set(handles.xAxVal,'Value',handles.xSliceN);
%     set(handles.xAxVal, 'SliderStep', [1/handles.Xsl, 10/handles.Xsl]);
%     set(handles.xVl,'String',num2str(handles.xSliceN));
%     
%     set(handles.yAxVal,'Min',1);
%     set(handles.yAxVal,'Max',handles.Ysl);
%     set(handles.yAxVal,'Value',handles.ySliceN);
%     set(handles.yAxVal, 'SliderStep', [1/handles.Ysl, 10/handles.Ysl]);
%     set(handles.yVl,'String',num2str(handles.ySliceN));
%     
%     set(handles.zAxVal,'Min',1);
%     set(handles.zAxVal,'Max',handles.Zsl);
%     set(handles.zAxVal,'Value',handles.zSliceN);
%     set(handles.zAxVal, 'SliderStep', [1/handles.Zsl, 10/handles.Zsl]);
%     set(handles.zVl,'String',num2str(handles.zSliceN));
    
    makePlot(handles)
    
    if isequal(handles.step1.ForegroundColor,  [0 0 0])
        
        handles.mni_button.Enable = 'on';
        handles.native_button.Enable = 'on';
        
    end
    
    handles.actAopts.Enable = 'on';
    handles.colorB.Enable = 'on';
    handles.sizeB.Enable = 'on';
end

guidata(hObject, handles);




function handles = makePlot(handles)
cla(handles.showMR)
axes(handles.showMR)
hold on
Cs = slice(handles.brainIm,handles.xSliceN,[],[]);
CsquzData = squeeze(double(handles.brainIm(:,handles.xSliceN,:)));
hold on
set(Cs, 'alphadata', CsquzData, 'facealpha','interp');alim([0 0.5]);

Ss = slice(handles.brainIm,[],handles.ySliceN,[]);
SsquzData = squeeze(double(handles.brainIm(handles.ySliceN,:,:)));
set(Ss, 'alphadata', SsquzData, 'facealpha','interp');alim([0 0.5]);

As = slice(handles.brainIm,[],[],handles.zSliceN);
AsquzData = squeeze(double(handles.brainIm(:,:,handles.zSliceN)));
set(As, 'alphadata', AsquzData, 'facealpha','interp');alim([0 0.5]);

shading('interp')
colormap('bone')

hold on
[~] = Plot3D_ElecBound_01_GUI( handles.output_args );
hold on
scatter3(handles.x.C , handles.y.C , handles.z.C, handles.sz.C, handles.col.C, 'filled');
scatter3(handles.x.A , handles.y.A , handles.z.A, handles.sz.A, handles.col.A, 'filled');
scatter3(handles.x.P , handles.y.P , handles.z.P, handles.sz.P, handles.col.P, 'filled');
scatter3(handles.x.M , handles.y.M , handles.z.M, handles.sz.M, handles.col.M, 'filled');
scatter3(handles.x.L , handles.y.L , handles.z.L, handles.sz.L, handles.col.L, 'filled');

grid off
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'ZTickLabel',[])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])

xlim([0 handles.Xsl])
ylim([0 handles.Ysl])
zlim([0 handles.Zsl])

set(gca,'View', [-116.0000  22.0000])
set(gca,'Color','none')

rotate3d on


% --- Executes on button press in colorB.
function colorB_Callback(hObject, eventdata, handles)
% hObject    handle to colorB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of colorB
get(hObject)
if handles.colorB.Value == 1
    handles.color = 1;
else
    handles.color = 0;
end
[handles.x,handles.y,handles.z,handles.sz,handles.col] =...
    getNeuroGUIparams(handles.neuroDat, handles.output_args,handles.color,handles.size);
makePlot(handles)
guidata(hObject, handles);

% --- Executes on button press in sizeB.
function sizeB_Callback(hObject, eventdata, handles)
% hObject    handle to sizeB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sizeB
get(hObject)
if handles.sizeB.Value == 1
    handles.size = 1;
else
    handles.size = 0;
end
[handles.x,handles.y,handles.z,handles.sz,handles.col] =...
    getNeuroGUIparams(handles.neuroDat, handles.output_args,handles.color,handles.size);
makePlot(handles)
guidata(hObject, handles);


% --- Executes on selection change in ImageSet.
function ImageSet_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImageSet contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImageSet


% --- Executes during object creation, after setting all properties.
function ImageSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in mniLoad.
function mniLoad_Callback(hObject, eventdata, handles)
% hObject    handle to mniLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[mni_in, mni_path] = uigetfile(...
    {'*.nii;*.gz',...
     'NeuroImage (*.nii,*.gz)'},'Select MNI');
handles.mniMRI = mni_in;
handles.mniPath = mni_path;
handles.step1.ForegroundColor = [0 0 0];

handles.mniText.String = mni_in;

guidata(hObject, handles);

% --- Executes on button press in mri_pmni_Load.
function mri_pmni_Load_Callback(hObject, eventdata, handles)
% hObject    handle to mri_pmni_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[pMmri , mri_path] = uigetfile({'*.nii;*.gz',...
     'NeuroImage (*.nii,*.gz)'},'Select MRI');
handles.PM_MRI = pMmri;
handles.PM_mripath = mri_path;
handles.step2.ForegroundColor = [0 0 0];

handles.pmniText.String = pMmri;

% mr_nii = 'Case_260_post_brainextraction.nii';

guidata(hObject, handles);


% --- Executes on button press in patMRInat.
function patMRInat_Callback(hObject, eventdata, handles)
% hObject    handle to patMRInat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[pNmri , mri_path] = uigetfile({'*.nii;*.gz',...
     'NeuroImage (*.nii,*.gz)'},'Select MRI');
handles.PN_MRI = pNmri;
handles.PN_mripath = mri_path;
handles.step3.ForegroundColor = [0 0 0];

handles.pNatText.String = pNmri;

% mr_nii = 'Case_260_post_brainextraction.nii';

guidata(hObject, handles);


% --- Executes on button press in ele_mni_Load.
function ele_mni_Load_Callback(hObject, eventdata, handles)
% hObject    handle to ele_mni_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[EM_in , ele_path] = uigetfile({'*.nii;*.gz',...
     'ImageMask (*.nii,*.gz)'},'Select Electrode Mask');
handles.EM_Mask = EM_in;
handles.EM_path = ele_path;
handles.step4.ForegroundColor = [0 0 0];

handles.elemniText.String = EM_in;

% ele_nii = 'Case_260_post_R_electrode_isolation.nii.gz';

guidata(hObject, handles);


% --- Executes on button press in ele_nat_Load.
function ele_nat_Load_Callback(hObject, eventdata, handles)
% hObject    handle to ele_nat_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[EN_in , ele_path] = uigetfile({'*.nii;*.gz',...
     'ImageMask (*.nii,*.gz)'},'Select Electrode Mask');
handles.EN_Mask = EN_in;
handles.EN_path = ele_path;
handles.step5.ForegroundColor = [0 0 0];

handles.eleNatText.String = EN_in;

% ele_nii = 'Case_260_post_R_electrode_isolation.nii.gz';

guidata(hObject, handles);

% --- Executes on button press in ephysLoad.
function ephysLoad_Callback(hObject, eventdata, handles)
% hObject    handle to ephysLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ephys_in , ephys_path] = uigetfile('*.csv','Select Ephys CSV');
handles.ephysCSV = ephys_in;
handles.ephysPath = ephys_path;
handles.step6.ForegroundColor = [0 0 0];

% ele_nii = 'Case_260_post_R_electrode_isolation.nii.gz';

guidata(hObject, handles);

% --- Executes on button press in options.
function options_Callback(hObject, eventdata, handles)
% hObject    handle to options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



prompt = {'Data Bubble Color: s = solid (d) , g = gradient',...
    'Data Bubble Size: f = fixed (d) , r = relative' ,...
    'Display borders: 1 or 0 (d)',...
    'Electrode mask type: 1-3',...
    'Electrode diameter: 1.3 (d)'
    };
title = 'OPTIONS (d) = defaults';
dims = [1 50];
definput = {'s','f','0','1','1.3'};

% FIGURE OUT A WAY to UPDATE

answer = inputdlg(prompt,title,dims,definput);

handles.step7.ForegroundColor = [0 0 0];


guidata(hObject, handles);


% --- Executes when selected object is changed in mniORnat.
function mniORnat_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in mniORnat 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(hObject.String,'MNI')
    handles.mri2use = 'mni';
    run_Callback(hObject, eventdata, handles)
    
elseif strcmp(hObject.String,'Native')
    handles.mri2use = 'native';
    run_Callback(hObject, eventdata, handles)
    
    
end

handles.step8.ForegroundColor = [0 0 0];


guidata(hObject, handles);


% --------------------------------------------------------------------
function defaultT_Callback(hObject, eventdata, handles)
% hObject    handle to defaultT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


figLOC = which('DBS_PostMR_MERoverlay_v1.m');

[defaultLoc , ~ , ~ ] = fileparts(figLOC);

handles.mniMRI = 'MNI152lin_T1_1mm_brain.nii.gz';
handles.mniPath = [defaultLoc , '\'];
handles.step1.ForegroundColor = [0 0 0];
handles.mniText.String = 'MNI152lin_T1_1mm_brain.nii.gz';

handles.PM_MRI = 'c260_MNIbrain.nii.gz';
handles.PM_mripath = [defaultLoc , '\'];
handles.step2.ForegroundColor = [0 0 0];
handles.pmniText.String = 'c260_MNIbrain.nii.gz';

handles.PN_MRI = 'c260_NATbrain.nii';
handles.PN_mripath = [defaultLoc , '\'];
handles.step3.ForegroundColor = [0 0 0];
handles.pNatText.String = 'c260_NATbrain.nii';

handles.EM_Mask = 'c260_eleMNI.nii.gz';
handles.EM_path = [defaultLoc , '\'];
handles.step4.ForegroundColor = [0 0 0];
handles.elemniText.String = 'c260_eleMNI.nii.gz';

handles.EN_Mask = 'c260_NATele.nii.gz';
handles.EN_path = [defaultLoc , '\'];
handles.step5.ForegroundColor = [0 0 0];
handles.eleNatText.String = 'c260_NATele.nii.gz';

handles.ephysCSV = 'neurodata.csv';
handles.ephysPath = [defaultLoc , '\'];
handles.step6.ForegroundColor = [0 0 0];

handles.actAopts.Enable = 'on';

guidata(hObject, handles);


% --------------------------------------------------------------------
function setDatapath_Callback(hObject, eventdata, handles)
% hObject    handle to setDatapath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataDIR = uigetdir;

checkisDir = 1;
while checkisDir
    if dataDIR == 0
        dataDIR = uigetdir;
    else
        checkisDir = 0;
    end
end

addpath(genpath(dataDIR));
handles.DataDIR = [dataDIR , '\'];

guidata(hObject, handles);
