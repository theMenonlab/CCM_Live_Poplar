function varargout = leanDataEditor(varargin)
% LEANDATAEDITOR MATLAB code for leanDataEditor.fig
%      LEANDATAEDITOR, by itself, creates a new LEANDATAEDITOR or raises the existing
%      singleton*.
%
%      H = LEANDATAEDITOR returns the handle to a new LEANDATAEDITOR or the handle to
%      the existing singleton*.
%
%      LEANDATAEDITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEANDATAEDITOR.M with the given input arguments.
%
%      LEANDATAEDITOR('Property','Value',...) creates a new LEANDATAEDITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before leanDataEditor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to leanDataEditor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help leanDataEditor

% Last Modified by GUIDE v2.5 07-Jul-2017 11:35:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @leanDataEditor_OpeningFcn, ...
                   'gui_OutputFcn',  @leanDataEditor_OutputFcn, ...
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

%--------------------------------------------------------------------------
function updateImages(handles)
    axes(handles.Axes_imgRaw);
    imagesc(handles.frames(handles.fInd).imgRaw);                           colormap gray; axis off equal xy;
    
    axes(handles.Axes_imgPro);
    imagesc(handles.frames(handles.fInd).imgPro);                           colormap gray; axis off equal xy;
    
    axes(handles.Axes_objRaw);
    imagesc(handles.frames(handles.fInd).objRaw(:,:,handles.zInd));         colormap gray; axis off equal;
    
    axes(handles.Axes_objPro);
    imagesc(handles.frames(handles.fInd).objPro(:,:,handles.zInd));         colormap gray; axis off equal xy;
    
    axes(handles.Axes_objRef);
    imagesc(handles.frames(handles.fInd).objRef);                           colormap gray; axis off equal xy;
    
    figOpt = handles.figOpt;
    % Additional figure
    if figOpt == 1
        % Pixel mean array
        axes(handles.Axes_oneMore); %hold on;
        pixelMeanArr = extractfield(handles.frames,'pixelMean'); 
        plot(pixelMeanArr,'--b'); xlabel('Frame #'); ylabel('Pixel Mean (a.u)'); axis tight;
        hold on;
        plot(handles.fInd,pixelMeanArr(handles.fInd),'ro'); 
        hold off;
        %set(gca,'Position',[0 0 100 100]);
        %hold off;   
    elseif figOpt == 2
        % Reconstructed Cannula Image
        axes(handles.Axes_oneMore);
        imagesc(handles.frames(handles.fInd).imgRec); colormap gray; axis off equal xy;
    elseif figOpt == 3
        % Residue
        axes(handles.Axes_oneMore);
        imagesc(handles.frames(handles.fInd).imgRes); colormap gray; axis off equal xy;
    elseif figOpt == 4
        % Residue array
        axes(handles.Axes_oneMore);
        residueMeanArr = extractfield(handles.frames,'residueMean');
        plot(residueMeanArr,'--b'); xlabel('Frame #'); ylabel('Residue Mean (a.u)'); axis tight;
        hold on;
        plot(handles.fInd,residueMeanArr(handles.fInd),'ro'); 
        hold off;
    elseif figOpt == 5
        % SVD decomposition
        axes(handles.Axes_oneMore);
        cal = handles.calibration;
        lambda = handles.frames(handles.fInd).ctrlFactor;
        y = double(handles.frames(handles.fInd).imgPro(cal.imgFOV));
        
        if handles.frames(handles.fInd).method == 0
            svdPlot = (cal.s(:,1).*(cal.U'*y))./(cal.s.^2 + lambda^2);
        elseif handles.frames(handles.fInd).method == 2
            svdPlot = (cal.U'*y)./cal.s;
            k = round(handles.frames(handles.fInd).ctrlFactor*length(handles.calibration.s)/100);
            svdPlot(k+1:end) = 0;
        end
        
        semilogy(svdPlot); xlabel('SVD #'); ylabel('Intensity'); axis tight;
    elseif figOpt == 6
        axes(handles.Axes_oneMore);
        imagesc(handles.calibration.DC); colormap gray; axis off equal xy;
    end
    
    % Position update
    set(handles.Edit_xPos,'String',handles.frames(handles.fInd).xPos);
    set(handles.Edit_yPos,'String',handles.frames(handles.fInd).yPos);
    set(handles.Edit_zPos,'String',handles.frames(handles.fInd).zPos);
    set(handles.Edit_time,'String',handles.frames(handles.fInd).timeStamp);
    set(handles.Edit_pixelMean,'String',handles.frames(handles.fInd).pixelMean);
    

function handles = processReconImageFunc(handles)
    
    f   = handles.fInd;
    frames = handles.frames;
    calibration.zSize = handles.zIndMax;
    processObject;
    
    if handles.denoise
        %Denoise reconstruction image
        [thr,sorh,keepapp] = ddencmp('den','wv',handles.frames(f).objPro);
        frames(f).objPro = wdencmp('gbl',handles.frames(f).objPro,'sym4',2,thr*handles.denoiseThr,'h',1);
    end
    
    handles.frames = frames;
    
function handles = reconstructImageFunc(handles)
    
    f   = handles.fInd;
    frames = handles.frames;
    calibration = handles.calibration;
    
    processImage;
    reconstructObject;
    
    % Compute reconstructed image and residue
    %for f = 1:handles.frameMax
        %frames(f).objRaw=zeros(size(calibration.objFOV));
        yR = calibration.A*frames(f).objRaw(calibration.objFOV);
        frames(f).imgRec = zeros(size(calibration.imgFOV));
        frames(f).imgRec(calibration.imgFOV==1) = yR;
        %frames(f).imgPro=zeros(size(frames(f).imgRec));
        frames(f).imgRes = frames(f).imgRec - frames(f).imgPro;
        frames(f).residueMean = mean(frames(f).imgRes(calibration.imgFOV));
    %end 
    
    
    %handles.psfSize = psfSize;
    handles.frames = frames;
    
function handles = loadToHandles(handles, filename)
    
    vars = whos('-file',filename);
    % Check if it's cal file
    if ismember('A', {vars.name})
        errordlg('This is a calibration file. Load data file.');
        return;
    elseif ~ismember('calID', {vars.name})
        errordlg('This is not a proper data file.');
        return;
    end
    
    % Load only if it's a data file
    load(filename);
    
    if ~ismember('frames', {vars.name})
        saveObjPro = 0;
        saveImgPro = 1;
        for f = 1:size(imgRawFrames,2)
            frames(f).imgRaw        = imgRawFrames{f};
            frames(f).imgPro        = imgProFrames{f};
            frames(f).objRaw        = objRawFrames{f};
            frames(f).exposure      = Exp;
            frames(f).timeStamp     = timeStamp(f);
            frames(f).sfX           = sfX;
            frames(f).sfY           = sfY;
            frames(f).ctrlFactor    = ctrlFactor;
            frames(f).gain          = gain;
            frames(f).threshold     = threshold;
            frames(f).brightness    = brightness;
            frames(f).convOpt       = convOpt;
            frames(f).exCom         = exCom;
            frames(f).method        = method;
            calPath = saveCal;
            
            % Masking feature was not available at first
            if exist('maskOpt','var')
                frames(f).maskOpt       = maskOpt;
                frames(f).maskSize      = maskSize;
                frames(f).maskShiftX    = maskShiftX;
                frames(f).maskShiftY    = maskShiftY;
            else
                frames(f).maskOpt       = 0;
                frames(f).maskSize      = 100;
                frames(f).maskShiftX    = 0;
                frames(f).maskShiftY    = 0;
            end
            
            % objRef used to have different names
            if exist('objRefFrames','var')
                frames(f).objRef        = objRefFrames{f};
            elseif exist('refRawFrames','var')
                frames(f).objRef        = refRawFrames{f};
            else
                frames(f).objRef        = 0;
            end
            
            % Older save files didn't save xyz position
            if exist('xPosFrames','var')
                frames(f).xPos          = xPosFrames(f);
                frames(f).yPos          = yPosFrames(f);
                frames(f).zPos          = zPosFrames(f);
            else
                frames(f).xPos          = 0;
                frames(f).yPos          = 0;
                frames(f).zPos          = 0;
            end
            
        end
        
    end
    
    handles.thisFile = filename;
    set(handles.figure1, 'Name', filename);
    
    handles.frames      = frames;
    handles.calID       = calID;
    handles.calPath     = calPath;
    handles.saveObjPro  = saveObjPro;
    if ~exist('saveImgPro','var')
        saveImgPro = 1;
    end
    handles.saveImgPro  = saveImgPro;
    
    % Misc
    handles.fInd       = 1;
    handles.zInd       = 1;
    handles.frameMax   = size(frames,2);
    handles.psfSize    = size(frames(1).objRaw,2);
    handles.zIndMax    = size(frames(1).objRaw,3);
    handles.msPerFrame = 0;
    %handles.msPerFrame = round(seconds(datetime(frames(handles.frameMax).timeStamp) ... 
    %                         - datetime(frames(1).timeStamp))/handles.frameMax*1000);
    
    % Post-processing
    if exist('denoiseThr','var')==0
        handles.denoiseThr = 0;
    else
        handles.denoiseThr = denoiseThr;
    end
        
    for f = 1:handles.frameMax
        % Recreate objPro
        if ~saveObjPro
            handles.fInd = f;
            handles = processReconImageFunc(handles);
        end
        
        % calculate pixelMean
        imgPro = handles.frames(f).imgPro;
        y = double(handles.frames(f).imgPro(handles.frames(f).imgPro~=0));
        handles.frames(f).pixelMean = mean(y);
    end
    
    % Reset fInd
    handles.fInd = 1;
    
    % Set control values to the saved parameters
    setFrameControlValues(handles);
    updateImages(handles);
        
function setFrameControlValues(handles)
    fInd = handles.fInd;
    frames = handles.frames;
    
    % Set control values to the saved parameters
    set(handles.Edit_shiftX,'String',frames(fInd).sfX);
    set(handles.Edit_shiftY,'String',frames(fInd).sfY);
    set(handles.Edit_ctrlFactor,'String',frames(fInd).ctrlFactor);
    set(handles.Popup_method,'Value',frames(fInd).method+1);
    
    set(handles.Checkbox_convolution,'Value',frames(fInd).convOpt);
    set(handles.Edit_gain,'String',frames(fInd).gain);
    set(handles.Slider_gain,'Value',frames(fInd).gain);
    set(handles.Edit_threshold,'String',frames(fInd).threshold);
    set(handles.Slider_threshold,'Value',frames(fInd).threshold);
    set(handles.Edit_brightness,'String',frames(fInd).brightness);
    set(handles.Slider_brightness,'Value',frames(fInd).brightness);
    
    set(handles.Edit_frame,'String',handles.fInd);
    set(handles.Slider_frame,'Value',handles.fInd);
    set(handles.Slider_frame,'Min',1);
    set(handles.Slider_frame,'Max',handles.frameMax);
    set(handles.Slider_frame,'SliderStep', [1/handles.frameMax , 10/handles.frameMax]);
    
    set(handles.Edit_zInd,'String',handles.zInd);
    set(handles.Slider_zInd,'Value',handles.zInd);
    set(handles.Slider_zInd,'Min',1);
    set(handles.Slider_zInd,'Max',handles.zIndMax);
    set(handles.Slider_zInd,'SliderStep', [1/handles.zIndMax , 10/handles.zIndMax]);
    
    set(handles.Checkbox_maskOpt,'Value',frames(fInd).maskOpt);
    set(handles.Edit_maskSize,'String',frames(fInd).maskSize);
    set(handles.Edit_maskShiftX,'String',frames(fInd).maskShiftX);
    set(handles.Edit_maskShiftY,'String',frames(fInd).maskShiftY);
    
    set(handles.Edit_msPerFrame,'String',handles.msPerFrame);
    set(handles.Edit_denoiseThr,'String',handles.denoiseThr);
    
    exp = frames(fInd).exposure;
    for i = 1:size(exp,2)
        expStr{i} = num2str(exp(i));
    end
    set(handles.Listbox_exposure,'String',expStr);
    

function saveHandleToFile(handles,filename)

    frames  = handles.frames;
    if exist('handles.calibration')
        if handles.calID ~= handles.newCalID
            calID   = handles.newCalID;
            calPath = handles.newCalPath;
        else
            calID   = handles.calID;
            calPath = handles.calPath;
        end
    else
        calID   = handles.calID;
        calPath = handles.calPath;
    end
    saveObjPro = handles.saveObjPro;
    saveData;

% --- Executes just before leanDataEditor is made visible. --------------------
function leanDataEditor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to leanDataEditor (see VARARGIN)

% Choose default command line output for leanDataEditor
    handles.output  = hObject;
    handles.denoise = 0;
    
    addpath([pwd]);
    
    handles.saveFig = figure('Position', [100, 100, 800, 200]);
    handles.figOpt  = 1;
    set(handles.saveFig,'visible','off');
    
    demofile = 'demo_v5.mat';
    handles = loadToHandles(handles,demofile);
    updateImages(handles);
    
    set(handles.Edit_shiftX,'Enable','off');
    set(handles.Edit_shiftY,'Enable','off');
    set(handles.Popup_method,'Enable','off');
    set(handles.Edit_ctrlFactor,'Enable','off');
    set(handles.Button_reconstructFrames,'Enable','off');
    
    set(handles.Edit_xPos,'Enable','off');
    set(handles.Edit_yPos,'Enable','off');
    set(handles.Edit_zPos,'Enable','off');
    set(handles.Edit_time,'Enable','off');
    
    addpath('C:\Cannula Microscopy\IVCCM_BPRB Software V1\Matlab Library\');

guidata(hObject, handles);

% --- Outputs from this function are returned to the command line. --------
function varargout = leanDataEditor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %rmpath('D:\Box Sync\Lab PC\CCM Software v6\Matlab Library\');
    handles.output = 1;

% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------
function Slider_threshold_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    
function Slider_threshold_Callback(hObject, eventdata, handles)
    
    handles.frames(handles.fInd).threshold = get(hObject,'Value');
    set(handles.Edit_threshold,'String',handles.frames(handles.fInd).threshold);
    handles = processReconImageFunc(handles);
    updateImages(handles);
    
guidata(hObject, handles);

function Edit_threshold_Callback(hObject, eventdata, handles)
    
    handles.frames(handles.fInd).threshold = str2double(get(hObject,'String'));
    if handles.frames(handles.fInd).threshold > get(handles.Slider_threshold,'Max')
        handles.frames(handles.fInd).threshold = get(handles.Slider_threshold,'Max');
    elseif handles.frames(handles.fInd).threshold < get(handles.Slider_threshold,'Min')
        handles.frames(handles.fInd).threshold = get(handles.Slider_threshold,'Min');
    end
    set(handles.Edit_threshold,'String',handles.frames(handles.fInd).threshold);
    set(handles.Slider_threshold,'Value',handles.frames(handles.fInd).threshold);
    handles = processReconImageFunc(handles);
    updateImages(handles);
    
guidata(hObject, handles);    

function Edit_threshold_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
%--------------------------------------------------------------------------
function Slider_brightness_Callback(hObject, eventdata, handles)

    handles.frames(handles.fInd).brightness = get(hObject,'Value');
    set(handles.Edit_brightness,'String',handles.frames(handles.fInd).brightness);
    handles = processReconImageFunc(handles);    
    updateImages(handles);

guidata(hObject, handles);

function Slider_brightness_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function Edit_brightness_Callback(hObject, eventdata, handles)
    
    handles.frames(handles.fInd).brightness = str2double(get(hObject,'String'));
    if handles.frames(handles.fInd).brightness > get(handles.Slider_brightness,'Max')
        handles.frames(handles.fInd).brightness = get(handles.Slider_brightness,'Max');
    elseif handles.frames(handles.fInd).brightness < get(handles.Slider_brightness,'Min')
        handles.brightness = get(handles.Slider_brightness,'Min');
    end
    set(handles.Edit_brightness,'String',handles.frames(handles.fInd).brightness);
    set(handles.Slider_brightness,'Value',handles.frames(handles.fInd).brightness);
    handles = processReconImageFunc(handles);
    updateImages(handles);

guidata(hObject, handles);

function Edit_brightness_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%--------------------------------------------------------------------------
function Slider_gain_Callback(hObject, eventdata, handles)

    handles.frames(handles.fInd).gain = get(hObject,'Value');
    set(handles.Edit_gain,'String',handles.frames(handles.fInd).gain);
    handles = processReconImageFunc(handles);
    updateImages(handles);

guidata(hObject, handles);

function Slider_gain_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function Edit_gain_Callback(hObject, eventdata, handles)

    handles.frames(handles.fInd).gain = str2double(get(hObject,'String'));
    if handles.frames(handles.fInd).gain > get(handles.Slider_gain,'Max')
        handles.frames(handles.fInd).gain = get(handles.Slider_gain,'Max');
    elseif handles.frames(handles.fInd).gain < get(handles.Slider_gain,'Min')
        handles.frames(handles.fInd).gain = get(handles.Slider_gain,'Min');
    end
    set(handles.Edit_gain,'String',handles.frames(handles.fInd).gain);
    set(handles.Slider_gain,'Value',handles.frames(handles.fInd).gain);
    handles = processReconImageFunc(handles);
    updateImages(handles);
    
guidata(hObject, handles);

function Edit_gain_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
function Slider_frame_Callback(hObject, eventdata, handles)
    
    handles.fInd = round(get(hObject,'Value'));
    if handles.fInd > handles.frameMax
        handles.fInd = handles.frameMax;
    end
    
    set(handles.Slider_frame,'Value',handles.fInd);
    set(handles.Edit_frame,'String',handles.fInd);
    
    setFrameControlValues(handles)
    handles = processReconImageFunc(handles);
    updateImages(handles);
    
guidata(hObject, handles);

function Slider_frame_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Edit_frame_Callback(hObject, eventdata, handles)
    
    handles.fInd = round(str2double(get(hObject,'String')));
    if handles.fInd > get(handles.Slider_frame,'Max')
        handles.fInd = get(handles.Slider_frame,'Max');
    elseif handles.fInd < get(handles.Slider_frame,'Min')
        handles.fInd = get(handles.Slider_frame,'Min');
    end
    set(handles.Slider_frame,'Value',handles.fInd);
    set(handles.Edit_frame,'String',handles.fInd);
    setFrameControlValues(handles)
    handles = processReconImageFunc(handles);
    updateImages(handles);
    
guidata(hObject, handles);

function Edit_frame_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
function Button_reconstructFrames_Callback(hObject, eventdata, handles)
    
    fIndRef = handles.fInd;
    for f = 1:handles.frameMax
        handles.fInd = f;
        set(handles.Slider_frame,'Value',handles.fInd);
        set(handles.Edit_frame,'String',handles.fInd);
        
        handles.frames(handles.fInd).method     = handles.frames(fIndRef).method;
        handles.frames(handles.fInd).sfX        = handles.frames(fIndRef).sfX;
        handles.frames(handles.fInd).sfY        = handles.frames(fIndRef).sfY;
        handles.frames(handles.fInd).ctrlFactor = handles.frames(fIndRef).ctrlFactor;
        
        handles = reconstructImageFunc(handles);
        handles = processReconImageFunc(handles);
        
        updateImages(handles);
        drawnow;
        guidata(hObject,handles);
        
    end
guidata(hObject,handles);


%--------------------------------------------------------------------------
function Button_processFrames_Callback(hObject, eventdata, handles)

    fIndRef = handles.fInd;
    for f = 1:handles.frameMax
        handles.fInd = f;
        set(handles.Slider_frame,'Value',handles.fInd);
        set(handles.Edit_frame,'String',handles.fInd);
        
        handles.frames(handles.fInd).gain       = handles.frames(fIndRef).gain;
        handles.frames(handles.fInd).threshold  = handles.frames(fIndRef).threshold;
        handles.frames(handles.fInd).brightness = handles.frames(fIndRef).brightness;
        handles.frames(handles.fInd).maskOpt    = handles.frames(fIndRef).maskOpt;
        handles.frames(handles.fInd).maskSize   = handles.frames(fIndRef).maskSize;
        handles.frames(handles.fInd).maskShiftX = handles.frames(fIndRef).maskShiftX;
        handles.frames(handles.fInd).maskShiftY = handles.frames(fIndRef).maskShiftY;
        handles.frames(handles.fInd).convOpt    = handles.frames(fIndRef).convOpt;
        
        handles = processReconImageFunc(handles);
        updateImages(handles);
        
        drawnow;
        guidata(hObject,handles);
    end
guidata(hObject,handles);

%--------------------------------------------------------------------------
function Button_play_Callback(hObject, eventdata, handles)
    for f = 1:handles.frameMax
        handles.fInd = f;
        set(handles.Slider_frame,'Value',handles.fInd);
        set(handles.Edit_frame,'String',handles.fInd);
        handles = processReconImageFunc(handles);
        updateImages(handles);
        drawnow;
        guidata(hObject,handles);
        pause(handles.msPerFrame/1000);
    end
guidata(hObject,handles);

%--------------------------------------------------------------------------
function Button_subtract_Callback(hObject, eventdata, handles)
    
    minPixelMean = 999999999;
    avgFrame = zeros(size(handles.imgRawFrames{1}));
    for f = 1:handles.frameMax
        avgFrame = avgFrame + handles.imgRawFrames{f};
        
        imgPro = handles.imgProFrames{handles.frame};
        y = double(imgPro(imgPro~=0));
        pixelMean = mean(y);
        
        if minPixelMean > pixelMean
            minFrame = handles.imgRawFrames{f};
        end
    end
    avgFrame = avgFrame / handles.frameMax;
    
    for f = 1:handles.frameMax
        handles.frame = f;
        
        %handles.imgRawFrames{f} =  handles.imgRawFrames{f} - avgFrame;
        handles.imgRawFrames{f} =  handles.imgRawFrames{f} - minFrame;
        
        set(handles.Slider_frame,'Value',handles.frame);
        set(handles.Edit_frame,'String',handles.frame);
        handles = reconstructImageFunc(handles);
        handles = processReconImageFunc(handles);
        updateImages(handles);
        drawnow;
        guidata(hObject,handles);
    end
guidata(hObject,handles);

%--------------------------------------------------------------------------
function Edit_msPerFrame_Callback(hObject, eventdata, handles)
    handles.msPerFrame = round(str2double(get(hObject,'String')));
guidata(hObject,handles);

function Edit_msPerFrame_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
function Edit_shiftX_Callback(hObject, eventdata, handles)
    handles.frames(handles.fInd).sfX = str2double(get(hObject,'String'));
    handles = reconstructImageFunc(handles);
    handles = processReconImageFunc(handles);
    updateImages(handles);
guidata(hObject,handles);

function Edit_shiftX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
function Edit_shiftY_Callback(hObject, eventdata, handles)
    handles.frames(handles.fInd).sfY = str2double(get(hObject,'String'));
    handles = reconstructImageFunc(handles);
    handles = processReconImageFunc(handles);
    updateImages(handles);
guidata(hObject,handles);
    
function Edit_shiftY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
function Edit_ctrlFactor_Callback(hObject, eventdata, handles)
    handles.frames(handles.fInd).ctrlFactor = str2double(get(hObject,'String'));
    handles = reconstructImageFunc(handles);
    handles = processReconImageFunc(handles);
    updateImages(handles);
guidata(hObject,handles);

function Edit_ctrlFactor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function Menu_file_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Menu_saveData_Callback(hObject, eventdata, handles)
    
    set(handles.figure1, 'pointer', 'watch'); drawnow;
    saveHandleToFile(handles,handles.thisFile);    
    set(handles.figure1, 'pointer', 'arrow');

% --------------------------------------------------------------------
function Menu_saveAs_Callback(hObject, eventdata, handles)
    
    [name,path] = uiputfile({'*.mat'},'Save Computational Microscope data file as');
    set(handles.figure1, 'pointer', 'watch'); drawnow;
    if name
        cd(path);
        savefile = [path,name];
        saveHandleToFile(handles,savefile);
    end
    set(handles.figure1, 'pointer', 'arrow');

% --------------------------------------------------------------------
function Menu_exportFrame_Callback(hObject, eventdata, handles)
    
    [name,path] = uiputfile({'*.eps';'*.png';'*.pdf'},'Export the frame to');
    set(handles.figure1, 'pointer', 'watch'); drawnow;
    if name
        cd(path);
        savefile = [path,name];
        saveas(handles.figure1, savefile);
    end
    set(handles.figure1, 'pointer', 'arrow');
    
guidata(hObject, handles);

% --------------------------------------------------------------------
function Menu_exportVideo_Callback(hObject, eventdata, handles)
    
    [name,path] = uiputfile('*.mp4','Export video frames to');
    set(handles.figure1, 'pointer', 'watch'); drawnow;
    if name
        cd(path);
        savefile = [path,name];
        writerObj = VideoWriter(savefile,'MPEG-4');
        writerObj.FrameRate = 5;
        open(writerObj);

        for f = 1:handles.frameMax
            handles.fInd = f;
            set(handles.Slider_frame,'Value',handles.fInd);
            set(handles.Edit_frame,'String',handles.fInd);
            handles = processReconImageFunc(handles);
            updateImages(handles);
            drawnow;
            
            %set(0,'CurrentFigure',handles.saveFig);
            %subplot(1,2,1); imagesc(handles.imgProFrames{handles.frame}); title('Cannula image');           colormap gray; axis off equal xy;
            %subplot(1,2,2); imagesc(handles.objProFrames{handles.frame}); title('Reconstruction image');    colormap gray; axis off equal xy;
            %subplot(1,3,3); imagesc(handles.objRefFrames{handles.frame}); title('Reference image');         colormap gray; axis off equal;
            %frame = getframe(handles.saveFig);
            
            frame = getframe(handles.figure1);
            writeVideo(writerObj,frame);
        end
        close(writerObj);
    end
    set(handles.figure1, 'pointer', 'arrow');
guidata(hObject, handles);

% --------------------------------------------------------------------
function Menu_loadData_Callback(hObject, eventdata, handles)
    
    [name,path] = uigetfile('*.mat','Load Computational Microscope data file');
    set(handles.figure1, 'pointer', 'watch'); drawnow;
    if name
        cd(path);
        loadfile = [path,name];
        handles = loadToHandles(handles, loadfile);        
    end
    set(handles.figure1, 'pointer', 'arrow');
    
guidata(hObject,handles);

% --------------------------------------------------------------------
function Menu_loadCal_Callback(hObject, eventdata, handles)
    
    [name,path] = uigetfile('*.mat','Load Computational Microscope calibration file');
    
    set(handles.figure1, 'pointer', 'watch'); drawnow;
    drawnow;
    if name
        calfile = [path,name];
        
        % Check if it's cal file
        vars = whos('-file',calfile);
        if ~ismember('A', {vars.name})
            errordlg('Not a proper calibration file.');
            set(handles.figure1, 'pointer', 'arrow');
            return;
        end
        
        load(calfile);
        
        handles.calibration.A           = A;
        handles.calibration.U           = U;
        handles.calibration.s           = s;
        handles.calibration.V           = V;
        handles.calibration.DC          = DC;
        handles.calibration.calID       = calID;
        handles.calibration.imgFOV      = imgFOV;
        handles.calibration.objFOV      = objFOV;
        handles.calibration.objFOVInd   = find(objFOV==1);
        %handles.calibration.objFOVInd   = objFOVInd;
        if exist('zSize','var')
            handles.calibration.zSize       = zSize;
        else
            handles.calibration.zSize       = 1;
        end
        handles.calibration.psfSize     = psfSize;
        handles.calibration.normArray   = normArray;
        
        handles.newCalID   = calID;
        handles.newCalPath = calfile;
        handles.zIndMax = handles.calibration.zSize;
        
        setFrameControlValues(handles);
        set(handles.Edit_shiftX,'Enable','on');
        set(handles.Edit_shiftY,'Enable','on');
        set(handles.Popup_method,'Enable','on');
        set(handles.Edit_ctrlFactor,'Enable','on');
        set(handles.Button_reconstructFrames,'Enable','on');
        
        % Compute reconstructed image and residue
        for f = 1:handles.frameMax
            handles.frames(f).objRaw=zeros(size(handles.calibration.objFOV));
            yR = handles.calibration.A*handles.frames(f).objRaw(handles.calibration.objFOV);
            handles.frames(f).imgRec = zeros(size(handles.calibration.imgFOV));
            handles.frames(f).imgRec(handles.calibration.imgFOV==1) = yR;
            handles.frames(f).imgPro=zeros(size(handles.frames(f).imgRec));
            handles.frames(f).imgRes = handles.frames(f).imgRec - handles.frames(f).imgPro;
            handles.frames(f).residueMean = mean(handles.frames(f).imgRes(handles.calibration.imgFOV));
        end 
        
    end
    set(handles.figure1, 'pointer', 'arrow');
    
guidata(hObject,handles);

% --------------------------------------------------------------------
function Popup_method_Callback(hObject, eventdata, handles)
    
    str = get(hObject, 'String');
    val = get(hObject,'Value');
    
    switch str{val}
        case 'Tikhonov'
            selectedMethod = 0;
        case 'DBS' % User selects membrane.
            selectedMethod = 1;
        case 'TSVD'
            selectedMethod = 2;
        otherwise
            selectedMethod = 0;
    end
    handles.frames(handles.fInd).method = selectedMethod;
    handles = reconstructImageFunc(handles);
    handles = processReconImageFunc(handles);
    updateImages(handles);

guidata(hObject,handles)

function Popup_method_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function Checkbox_convolution_Callback(hObject, eventdata, handles)
    handles.frames(handles.fInd).convOpt = get(hObject,'Value');
    handles = processReconImageFunc(handles);
    updateImages(handles);
guidata(hObject,handles)

% --------------------------------------------------------------------
function Menu_quitGUI_Callback(hObject, eventdata, handles)
    close(handles.figure1);

% --------------------------------------------------------------------
function Checkbox_denoise_Callback(hObject, eventdata, handles)
    handles.denoise = get(hObject,'Value');
    
    %if handles.denoise
    %    handles.denoiseThr = 1;
    %    set(handles.Edit_denoiseThr,'String',handles.denoiseThr);
    %end
    
    handles = processReconImageFunc(handles);
    updateImages(handles);
    
guidata(hObject,handles)

% --------------------------------------------------------------------
function Edit_denoiseThr_Callback(hObject, eventdata, handles)
    handles.denoiseThr = str2double(get(hObject,'String'));
    handles = processReconImageFunc(handles);
    updateImages(handles);
guidata(hObject, handles);

function Edit_denoiseThr_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
%--------------------------------------------------------------------------
function Edit_pixelMean_Callback(hObject, eventdata, handles)

function Edit_pixelMean_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
function Edit_zInd_Callback(hObject, eventdata, handles)
    
    handles.zInd = round(str2double(get(hObject,'String')));
    if handles.zInd > get(handles.Slider_zInd,'Max')
        handles.zInd = get(handles.Slider_zInd,'Max');
    elseif handles.zInd < get(handles.Slider_zInd,'Min')
        handles.zInd = get(handles.Slider_zInd,'Min');
    end
    set(handles.Slider_zInd,'Value',handles.zInd);
    set(handles.Edit_zInd,'String',handles.zInd);
    handles = processReconImageFunc(handles);
    updateImages(handles);

function Edit_zInd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
function Slider_zInd_Callback(hObject, eventdata, handles)
    handles.zInd = round(get(hObject,'Value'));
    if handles.zInd > handles.zIndMax
        handles.zInd = handles.zIndMax;
    end
    set(handles.Slider_zInd,'Value',handles.zInd);
    set(handles.Edit_zInd,'String',handles.zInd);
    handles = processReconImageFunc(handles);
    updateImages(handles);
guidata(hObject, handles);

function Slider_zInd_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%--------------------------------------------------------------------------
function Edit_maskShiftX_Callback(hObject, eventdata, handles)
    handles.frames(handles.fInd).maskShiftX = str2double(get(hObject,'String'));
    set(handles.Edit_maskShiftX,'String',handles.frames(handles.fInd).maskShiftX);
    handles = processReconImageFunc(handles);
    updateImages(handles);
guidata(hObject, handles);

function Edit_maskShiftX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Edit_maskShiftY_Callback(hObject, eventdata, handles)
    handles.frames(handles.fInd).maskShiftY = str2double(get(hObject,'String'));
    set(handles.Edit_maskShiftY,'String',handles.frames(handles.fInd).maskShiftY);
    handles = processReconImageFunc(handles);
    updateImages(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Edit_maskShiftY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Checkbox_maskOpt.
function Checkbox_maskOpt_Callback(hObject, eventdata, handles)
    handles.frames(handles.fInd).maskOpt = get(hObject,'Value');
    handles = processReconImageFunc(handles);
    updateImages(handles);
guidata(hObject, handles);

function Edit_maskSize_Callback(hObject, eventdata, handles)
    handles.frames(handles.fInd).maskSize = str2double(get(hObject,'String'));
    set(handles.Edit_maskSize,'String',handles.frames(handles.fInd).maskSize);
    handles = processReconImageFunc(handles);
    updateImages(handles);
guidata(hObject, handles);

function Edit_maskSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Edit_xPos_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_xPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function Edit_xPos_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Edit_yPos_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_yPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function Edit_yPos_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Edit_zPos_Callback(hObject, eventdata, handles)

function Edit_zPos_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Edit_time_Callback(hObject, eventdata, handles)

function Edit_time_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Listbox_exposure.
function Listbox_exposure_Callback(hObject, eventdata, handles)

function Listbox_exposure_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Popup_extraFig.
function Popup_extraFig_Callback(hObject, eventdata, handles)
    val = get(hObject,'Value');
    
    if ~isfield(handles,'calibration') && val > 1
        val = 1;
        errordlg('Load calibration first.');
        set(hObject,'Value',val);
        return;
    end
    handles.figOpt = val;
    updateImages(handles);
guidata(hObject,handles);
    
function Popup_extraFig_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end