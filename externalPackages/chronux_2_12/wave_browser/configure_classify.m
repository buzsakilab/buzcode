function varargout = configure_classify(varargin)
% CONFIGURE_CLASSIFY M-file for configure_classify.fig
%      CONFIGURE_CLASSIFY, by itself, creates a new CONFIGURE_CLASSIFY or raises the existing
%      singleton*.
%
%      H = CONFIGURE_CLASSIFY returns the handle to a new CONFIGURE_CLASSIFY or the handle to
%      the existing singleton*.
%
%      CONFIGURE_CLASSIFY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONFIGURE_CLASSIFY.M with the given input arguments.
%
%      CONFIGURE_CLASSIFY('Property','Value',...) creates a new CONFIGURE_CLASSIFY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before configure_classify_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to configure_classify_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help configure_classify

% Last Modified by GUIDE v2.5 01-May-2006 23:22:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @configure_classify_OpeningFcn, ...
                   'gui_OutputFcn',  @configure_classify_OutputFcn, ...
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


% --- Executes just before configure_classify is made visible.
function configure_classify_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to configure_classify (see VARARGIN)


% Choose default command line output for configure_classify
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes configure_classify wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% configure_classify(0,7500,400,500,44100,[0.01 0.001],[3 5],[0 20000])

handles.lowerfreq = varargin{1};  % Lower frequency for zooming
handles.upperfreq = varargin{2}; % Upper frequency for zooming

handles.classified_height = varargin{3} ; % the height of the image in the classified axes
handles.classified_width = varargin{4};  % the width of the image in the classified axes

handles.Fs = varargin{5}; % Frequency of audio sampling per second
handles.movingwin = varargin{6}; % Size of the moving window in seconds; the first number is the window size and the second is the step size
handles.tapers = varargin{7}; % Tapers for smoothing
handles.fpass = varargin{8}; % Range of frequency sampling 

handles.fixed = varargin{9}; % Fixed scaling of the classified axes

set(handles.ZoomLowerFreq,'String',num2str(handles.lowerfreq));
set(handles.ZoomUpperFreq,'String',num2str(handles.upperfreq));
set(handles.ClassifiedWidth,'String',num2str(handles.classified_width));
set(handles.ClassifiedHeight,'String',num2str(handles.classified_height));
set(handles.Frequency,'String',num2str(handles.Fs));
set(handles.WinSize,'String',num2str(handles.movingwin(1) * 1000));
set(handles.StepSize,'String',num2str(handles.movingwin(2) * 1000));
set(handles.TW,'String',num2str(handles.tapers(1)));
set(handles.MinFreq,'String',num2str(handles.fpass(1)));
set(handles.MaxFreq,'String',num2str(handles.fpass(2)));

set(handles.FixedCheckbox,'Value',handles.fixed);

% set(handles.ZoomLowerFreq,'Enable','off');
% set(handles.ZoomUpperFreq,'Enable','off');
% set(handles.ClassifiedWidth,'Enable','off');
% set(handles.ClassifiedHeight,'Enable','off');
% set(handles.Frequency,'Enable','off');
% set(handles.WinSize,'Enable','off');
% set(handles.StepSize,'Enable','off');
% set(handles.TW,'Enable','off');
% set(handles.MinFreq,'Enable','off');
% set(handles.MaxFreq,'Enable','off');

uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = configure_classify_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

close;

function WinSize_Callback(hObject, eventdata, handles)
% hObject    handle to WinSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WinSize as text
%        str2double(get(hObject,'String')) returns contents of WinSize as a double


% --- Executes during object creation, after setting all properties.
function WinSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WinSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in OKButton.
function OKButton_Callback(hObject, eventdata, handles)
% hObject    handle to OKButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lowerfreq = str2num(get(handles.ZoomLowerFreq,'String'));  % Lower frequency for zooming
upperfreq = str2num(get(handles.ZoomUpperFreq,'String')); % Upper frequency for zooming

classified_height = str2num(get(handles.ClassifiedHeight,'String')); % the height of the image in the classified axes
classified_width = str2num(get(handles.ClassifiedWidth,'String'));  % the width of the image in the classified axes

Fs = str2num(get(handles.Frequency,'String')); % Frequency of audio sampling per second

winsizeS = str2num(get(handles.WinSize,'String')) / 1000;
stepS = str2num(get(handles.StepSize,'String')) / 1000;
movingwin = [winsizeS stepS]; % Size of the moving window in seconds; the first number is the window size and the second is the step size

tw = str2num(get(handles.TW,'String'));

fpasslower = str2num(get(handles.MinFreq,'String')); % Range of frequency sampling 
fpassupper = str2num(get(handles.MaxFreq,'String'));
fpass = [fpasslower fpassupper];

ierror = 1; % Indicates no errors encountered

if isempty(classified_height) || (classified_height < 1)
    ierror = 0;
end
    
if isempty(tw) || tw < 0
    ierror = 0;
end

if isempty(lowerfreq) || lowerfreq < 0
    ierror = 0;
end

if isempty(fpasslower) || fpasslower < 0
    ierror = 0;
end

if isempty(fpassupper) || fpassupper < fpasslower
    ierror = 0;
end

if isempty(upperfreq) || lowerfreq > upperfreq 
    ierror = 0;
end

if isempty(winsizeS) || winsizeS < 0
    ierror = 0;
end

if isempty(stepS) || stepS < 0
    ierror = 0
end

if isempty(tw) || tw < 0
    ierror = 0
else
    tapers = [tw,floor(2*tw-1)]; % Tapers for smoothing
end

fixed = get(handles.FixedCheckbox,'Value');

if ierror == 0
    ;
else
    handles.output = {lowerfreq,upperfreq,classified_height,classified_width,Fs,movingwin,tapers,fpass,fixed};
    guidata(hObject,handles);
    uiresume(handles.figure1);
end


%uiresume;
%close;

% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = 0;
guidata(hObject,handles);
uiresume(handles.figure1);

function StepSize_Callback(hObject, eventdata, handles)
% hObject    handle to StepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StepSize as text
%        str2double(get(hObject,'String')) returns contents of StepSize as a double


% --- Executes during object creation, after setting all properties.
function StepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function TW_Callback(hObject, eventdata, handles)
% hObject    handle to TW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TW as text
%        str2double(get(hObject,'String')) returns contents of TW as a double


% --- Executes during object creation, after setting all properties.
function TW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function MinFreq_Callback(hObject, eventdata, handles)
% hObject    handle to MinFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinFreq as text
%        str2double(get(hObject,'String')) returns contents of MinFreq as a double


% --- Executes during object creation, after setting all properties.
function MinFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function MaxFreq_Callback(hObject, eventdata, handles)
% hObject    handle to MaxFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxFreq as text
%        str2double(get(hObject,'String')) returns contents of MaxFreq as a double


% --- Executes during object creation, after setting all properties.
function MaxFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ZoomLowerFreq_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomLowerFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZoomLowerFreq as text
%        str2double(get(hObject,'String')) returns contents of ZoomLowerFreq as a double


% --- Executes during object creation, after setting all properties.
function ZoomLowerFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZoomLowerFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ZoomUpperFreq_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomUpperFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZoomUpperFreq as text
%        str2double(get(hObject,'String')) returns contents of ZoomUpperFreq as a double


% --- Executes during object creation, after setting all properties.
function ZoomUpperFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZoomUpperFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ClassifiedWidth_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifiedWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ClassifiedWidth as text
%        str2double(get(hObject,'String')) returns contents of ClassifiedWidth as a double


% --- Executes during object creation, after setting all properties.
function ClassifiedWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClassifiedWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ClassifiedHeight_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifiedHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ClassifiedHeight as text
%        str2double(get(hObject,'String')) returns contents of ClassifiedHeight as a double


% --- Executes during object creation, after setting all properties.
function ClassifiedHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClassifiedHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function Frequency_Callback(hObject, eventdata, handles)
% hObject    handle to Frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Frequency as text
%        str2double(get(hObject,'String')) returns contents of Frequency as a double


% --- Executes during object creation, after setting all properties.
function Frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in FixedCheckbox.
function FixedCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to FixedCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FixedCheckbox

handles.fixed = get(hObject,'Value');
guidata(gcbo,handles);
