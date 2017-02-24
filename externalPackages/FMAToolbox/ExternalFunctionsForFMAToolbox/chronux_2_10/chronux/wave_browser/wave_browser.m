function varargout = wave_browser(varargin)
% WAVE_BROWSER M-file for wave_browser.fig
%      WAVE_BROWSER, by itself, creates a new WAVE_BROWSER or raises the existing
%      singleton*.
%
%      H = WAVE_BROWSER returns the handle to a new WAVE_BROWSER or the handle to
%      the existing singleton*.
%
%      WAVE_BROWSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVE_BROWSER.M with the given input arguments.
%
%      WAVE_BROWSER('Property','Value',...) creates a new WAVE_BROWSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wave_browser_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wave_browser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help wave_browser

% Last Modified by GUIDE v2.5 29-May-2007 16:30:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wave_browser_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_browser_OutputFcn, ...
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


% --- Executes just before wave_browser is made visible.
function wave_browser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wave_browser (see VARARGIN)

handles.firsttime = 0; % indicates the firsttime that segment has been precomputed
handles.precomputed_spec = 0; % indicates that the spectra has not been precomputed
handles.longfile = 0;  % indicates whether the file is a long file
handles.maxwavsize = 10 * 44100; % I will have to explore what number works best here
handles.maxspec_t = 30; % duration of the max size of a spectra

handles.Fs = 44100; % default size to start with

handles.segments = []; % holds regular segments in the current chunk
handles.allsegments = []; % holds segments across the maximum wave size
handles.loadedsegment = 0; % indicates no segments have been loaded

handles.lastmarkerstart = 1; % largest segment

handles.segmentmode = 0; % by default start off with segmenting turned off
handles.dontcutsegments = 0; % by default do not adapt to segments

handles.automethod = 'threshold'; % use threshold or ratiof method

handles.indexthresh = 10; % for ration method the threshold which to cut the curve off
handles.lower_range = [10 10000];  % the numerator in the ratio
handles.upper_range = [15000 20000]; % the denomitor in the ratio

handles.nsmooth = 0; % moving average parameter for the thresholds curves

positionP = get(handles.OptionsUiPanel,'Position');
positionF = get(gcf,'Position');

positionF(3) = positionF(3) - positionP(3);

% set(gcf,'Position',positionF); % untested

% Choose default command line output for wave_browser
handles.output = hObject;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wave_browser wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wave_browser_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

function Frequency_Callback(hObject, eventdata, handles)
handles.Fs = eval(get(hObject,'String'));
guidata(gcbo,handles); 

function Frequency_CreateFcn(hObject, eventdata, handles)
set(hObject,'String', '44100');
handles.Fs = 44100;
guidata(gcbo,handles);
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in LoadFile.
function LoadFile_Callback(hObject, eventdata, handles)

[fname pname]=uigetfile({'*.wav';'*.*'},'Load Time Series');

if fname == 0
    return
end



set(handles.FileNameString, 'String',fname);
handles.filename = [pname fname];
[path] = cell2mat(regexp( handles.filename, '^.*\', 'match' ));
[extension] = cell2mat(regexp( handles.filename, '\w+$', 'match' ));
set(handles.Path,'String',path);
% set(handles.Extensions,'String',extension);

handles.segments = [];
handles.allsegments = [];

handles = loadfile(hObject, eventdata, handles);
guidata(hObject,handles);

function [wavesize channels] = wavsizeget(filename);
% Provides usable information about a file
    wavesize = 0;
    [filestatus info] = wavfinfo(filename);
    info = regexp(info,'[0-9]+','match');
    channels = str2num(info{2});
    wavesize = str2num(info{1});
        
function handles = loadfile(hObject, eventdata, handles, varargin)
% Function for loading a file or using the optional varargin to load a
% specified position and size in the file

% contents=get(handles.endian,'String');
% precision=contents{get(handles.endian,'Value')};

[datasize channels] = wavsizeget(handles.filename);
handles.wavsize = datasize; % total number of samples in the files

try
    handles.maxwavsize = round(handles.Fs * str2num(get(handles.MaximumWavSize,'String')));
catch
    handles.maxwavsize = 20;
    set(handles.MaximumWavSize,'String',num2str(handles.maxwavsize));
    handles.maxwavsize = handles.maxwavsize * handles.Fs;
end

if isempty(varargin)
    handles.markerstart = 1;
    if datasize > handles.maxwavsize
        handles.markerend = handles.maxwavsize;
        handles.longfile = 1; % indicates that the file is long and will be loaded in chunks
    else
        handles.markerend = datasize;
    end
else % passed in optional parameter
    handles.markervec = varargin{1};
    handles.markerstart = handles.markervec(1);
    handles.markerend = handles.markervec(2);
end

if handles.markerstart <= 1 % make sure the range is possible
    handles.markerstart = 1;   
    set(handles.PreviousChunk,'Enable','off');
else
    set(handles.PreviousChunk,'Enable','on');
end

if handles.markerend >= handles.wavsize
    handles.markerend = handles.wavsize;   
    set(handles.NextChunk,'Enable','off');
else
    set(handles.NextChunk,'Enable','on');
end


if handles.maxwavsize < handles.wavsize
   total_chunk = ceil(handles.wavsize / handles.maxwavsize);
   i = 1;
   while (i < total_chunk) &&  (handles.markerend >= (handles.maxwavsize * i))
       i = i + 1;
   end
   
   current_chunk = i-1;  
   
   if handles.markerend == handles.wavsize
       current_chunk = total_chunk;
   end
   
   set(handles.ChunkText,'String',['Chunk ' num2str(current_chunk) '/' num2str(total_chunk)]);
end

try
    handles.maxseglength = round(handles.Fs * str2num(get(handles.MaxSegLength,'String')));
catch
    handles.maxseglength = handles.Fs * 1;
end

set(handles.RealDuration,'String',num2str(handles.wavsize/handles.Fs,'%.1f'));

if handles.segmentmode % only if in segment mode make sure segments are not cut
    if (handles.markerstart - handles.lastmarkerstart > 0) % only do this in terms of forward movement
        if handles.dontcutsegments % this code is added so segments are not cut off when segmenting
            if not(isempty(handles.segments)) % at least one segment has been defined previously
                maxsegend = handles.segments(1).end; % find last defined segment in previous view
                for i = 2:length(handles.segments)
                    if handles.segments(i).end > maxsegend
                        maxsegend = handles.segments(i).end;
                    end
                end
        
                maxsegend = round(maxsegend * handles.Fs);
        
                if (handles.lastmarkerend - maxsegend) < (handles.lastmarkerend - handles.maxseglength)
                    handles.markerstart = (handles.lastmarkerstart + maxsegend) + 1; % defined segment is closer to the end
                else
                    handles.markerstart = handles.lastmarkerend - handles.maxseglength;
                end
            else
                handles.markerstart = handles.lastmarkerend - handles.maxseglength;
            end
    
            handles.markerend = handles.markerstart + handles.maxwavsize - 1;
            if handles.markerend > handles.wavsize 
                handles.markerend = handles.wavsize;
            end
        end
    end
end

hw=waitbar(0,'Loading ...'); waitbar(0.5,hw); drawnow;

[handles.markerstart handles.markerend]/handles.Fs
[handles.ts,handles.Fs] = wavread(handles.filename, [handles.markerstart handles.markerend]);
channel = str2double(get(handles.channel,'String'));
handles.ts = handles.ts(:,channel);

count = length(handles.ts);
handles.ts = handles.ts/std(handles.ts); % variance normalisation
set(handles.Frequency,'String', num2str(handles.Fs));
set( handles.Duration, 'String', count/handles.Fs );
Tim=eval(get(handles.DisplayWindow,'String'));
display_frac = 1;%max(1,Tim*handles.Fs/count);

set( handles.slider1, 'Value', 0 );
set( handles.SegmentButton, 'Enable', 'on' );

if handles.longfile
    set(handles.SeekButton,'Enable', 'on');
end

set(handles.LoadNext, 'Enable', 'on' );
set(handles.PlayAll, 'Enable', 'on' );
set(handles.PlayWindow, 'Enable', 'on' );
set(handles.Plot, 'Enable', 'on' );
set(handles.PlotAllButton, 'Enable', 'on');
set(handles.Precompute, 'Enable','on');
set(handles.Jump,'Enable','on');
set(handles.JumpBack,'Enable','on');

handles.segments = []; % remove the current segments
handles.segments = filtersegments(handles,handles.allsegments);

set(handles.Precompute,'Value',1); % Set into precompute mode

Precompute_Callback(handles.Precompute, eventdata, handles);
handles = guidata(gcbo);

% set(handles.Precompute,'Value',1);

handles.precomputed_spec = 1;
handles.dontcutsegments = 1; % make sure segments are not cut off
handles.lastmarkerstart = handles.markerstart;
handles.lastmarkerend = handles.markerend;

close(hw);
% guidata(gcbo,handles); 


% Plot_Callback(hObject, eventdata, handles);

% --- Executes on selection change in endian.
function endian_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function endian_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function edit2_Callback(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function FileNameString_Callback(hObject, eventdata, handles)

function FileNameString_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function WinSize_Callback(hObject, eventdata, handles)

function WinSize_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)

function slider1_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function StepSize_Callback(hObject, eventdata, handles)

function StepSize_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function TW_Callback(hObject, eventdata, handles)

function TW_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function DisplayWindow_Callback(hObject, eventdata, handles)

function DisplayWindow_CreateFcn(hObject, eventdata, handles)
set(hObject, 'String', '4');
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%function axes2ButtonDownCallback(hObject, eventdata, handles)
%h=handles.axesS; P=get(h,'CurrentPoint');
%fprintf( 'worked %f %f!\n', P(1),P(2));

function indexinS = getindexpre_c(t,timestart,timeend)
% A function for indexing correctly in time to a spectra stored in memory
tlen = length(t);

i=1;

while (i <= tlen) && (t(i) < timestart)
    i = i + 1;
end

firstindex = i;

while (i <= tlen) && (t(i) < timeend)
    i = i + 1;
end

secondindex = i;

indexinS = [firstindex secondindex];
% --- Executes on button press in Plot.
function Plot_Callback(hObject, eventdata, handles)

hw=waitbar(0.5,'Spectrogram calculation');drawnow

params.Fs=handles.Fs;

window=eval(get(handles.WinSize,'String'));
winstep=eval(get(handles.StepSize,'String'));
movingwin=[window winstep]*0.001;

fmin=eval(get(handles.MinFreq,'String'));
fmax=eval(get(handles.MaxFreq,'String'));
params.fpass=[fmin fmax];

p=eval(get(handles.TW,'String'));
params.tapers=[p floor(2*p-1)];

params.pad=1;

Tslider=get(handles.slider1,'Value');
Tim=eval(get(handles.DisplayWindow,'String'));
NT=min(round(Tim*handles.Fs),length(handles.ts));
handles.Tmin=1+floor(Tslider*length(handles.ts));
handles.Tmax=min(handles.Tmin+NT,length(handles.ts));

if handles.Tmax < length(handles.ts)
    set( handles.Jump, 'Enable', 'on' );
else 
    set( handles.Jump, 'Enable', 'off' );
end
if handles.Tmin > 1
    set( handles.JumpBack, 'Enable', 'on' );
else
    set( handles.JumpBack, 'Enable', 'off' );    
end

data=handles.ts(handles.Tmin:handles.Tmax);data=data(:);

handles.upper_range =  eval(get(handles.RatioLower,'String'));
handles.lower_range  = eval(get(handles.RatioUpper,'String'));
handles.indexthresh =  eval(get(handles.RatioThresh,'String'));
handles.nsmooth = eval(get(handles.SmoothFactor,'String'));

% determine spectrum type

contents=get(handles.SpectrumType,'String');
stype=contents{get(handles.SpectrumType,'Value')};

axes(handles.axesW); plot(((handles.markerstart - 1)/handles.Fs) + [handles.Tmin:handles.Tmax]/handles.Fs,handles.ts(handles.Tmin:handles.Tmax)); axis tight;

switch stype
    case 'Original'
        
     if not(handles.precomputed_spec) || handles.firsttime
        [S,t,f]=mtspecgramc(diff(data),movingwin,params);
        timeax=(handles.Tmin/handles.Fs)+t;
     else
         indexinS = getindexpre_c(handles.t,(handles.Tmin-1)/handles.Fs,(handles.Tmax-1)/handles.Fs);
%          indexinS = round(([handles.Tmin-1, handles.Tmax-1]/handles.Fs)/movingwin(2))+1;
         
         if indexinS(1) < 1
             indexinS(1) = 1;
         end
         
         SLen = length(handles.S(:,1));
         
         if indexinS(2) > SLen
             indexinS(2) = SLen;
         end
         
         f = handles.f;
         t = handles.t(indexinS(1):indexinS(2));
         S = handles.S(indexinS(1):indexinS(2),:);
         
         timeax=t;
     end
    
    cmap='default'; 
    
    th=eval(get(handles.AmpThresh,'String'));
    
    % This sets up the automatic segmenting algorithm
    if strcmp(handles.automethod,'threshold')
        [Stot boxcurve] = compute_threshold_free(S,th,handles.nsmooth);
        axes(handles.axesP); 
        semilogy(timeax,Stot); 
        axis tight;
    elseif strcmp(handles.automethod,'ratiof')
        [ratiof boxcurve] = compute_index(S,handles.lower_range,handles.upper_range,fmin,fmax,handles.indexthresh,handles.nsmooth);
        axes(handles.axesP); 
        semilogy(timeax,ratiof); 
        axis tight;
    end
    
     
    hold on; semilogy(timeax,boxcurve,'r'); hold off;
    axes(handles.axesS);
    imagesc(timeax,f,log(S)'); axis xy; colormap(cmap);
    %imagesc(t,f,log(S)'); axis xy; colormap(cmap);
    
%    set(h,'ButtonDownFcn',axes2ButtonDownCallback);
    
    case 'Time Derivative'
    
    if not(handles.precomputed_spec) || handles.firsttime   
        [S,t,f]=mtdspecgramc(diff(data),movingwin,0,params);S = S';
        timeax=handles.Tmin/handles.Fs+t;
    else
        indexinS = getindexpre_c(handles.t,(handles.Tmin-1)/handles.Fs,(handles.Tmax-1)/handles.Fs);
%          indexinS = round(([handles.Tmin-1, handles.Tmax-1]/handles.Fs)/movingwin(2))+1;
         
         if indexinS(1) < 1
             indexinS(1) = 1;
         end
         
         SLen = length(handles.S(1,:));
         
         if indexinS(2) > SLen
             indexinS(2) = SLen;
         end
         
         f = handles.f;
         t = handles.t(indexinS(1):indexinS(2));
         S = handles.S(:,indexinS(1):indexinS(2));
         timeax = t;
    end
    
    cmap='gray';
    th=eval(get(handles.TDerThresh,'String'));
    
    if strcmp(handles.automethod,'threshold')
        [Stot boxcurve] = compute_threshold_free(abs(S'),th.handles.nsmooth);
        axes(handles.axesP); 
        semilogy(timeax,Stot); 
        axis tight;
    elseif strcmp(handles.automethod,'ratiof')
        [ratiof boxcurve] = compute_index(abs(S)',handles.lower_range,handles.upper_range,fmin,fmax,handles.indexthresh,handles.nsmooth);
        axes(handles.axesP); 
        semilogy(timeax,ratiof); 
        axis tight;
    end
    
    hold on; semilogy(timeax,boxcurve,'r'); hold off;
    axes(handles.axesS);
    imagesc(timeax,f,S); axis xy; colormap(cmap);
    cmin=0.02*min(min(S)); cmax=0.02*max(max(S)); caxis([cmin cmax]);
    
    case 'Frequency Derivative'

    if not(handles.precomputed_spec) || handles.firsttime
        [S,t,f]=mtdspecgramc(diff(data),movingwin,pi/2,params);S=S';
        timeax=handles.Tmin/handles.Fs+t;
    else
        indexinS = getindexpre_c(handles.t,(handles.Tmin-1)/handles.Fs,(handles.Tmax-1)/handles.Fs);
         
         if indexinS(1) < 1
             indexinS(1) = 1;
         end
         
         SLen = length(handles.S(1,:));
         
         if indexinS(2) > SLen
             indexinS(2) = SLen;
         end
         
         f = handles.f;
         t = handles.t(indexinS(1):indexinS(2));
         S = handles.S(:,indexinS(1):indexinS(2));
         timeax = t;
    end
    
    cmap='gray';
    th=eval(get(handles.TDerThresh,'String'));
    
    if strcmp(handles.automethod,'threshold')
        [Stot boxcurve] = compute_threshold_free(abs(S'),th,handles.nsmooth);
        axes(handles.axesP); 
        semilogy(timeax,Stot); 
        axis tight;
    elseif strcmp(handles.automethod,'ratiof')
        [ratiof boxcurve] = compute_index(abs(S)',handles.lower_range,handles.upper_range,fmin,fmax,handles.indexthresh,handles.nsmooth);
        axes(handles.axesP); 
        semilogy(timeax,ratiof); 
        axis tight;
    end
    
    hold on; semilogy(timeax,boxcurve,'r'); hold off;
    axes(handles.axesS);
    imagesc(timeax,f,S); axis xy; colormap(cmap);
    cmin=0.02*min(min(S)); cmax=0.02*max(max(S)); caxis([cmin cmax]);

end;

if handles.firsttime % first time precomputing the spectra
    handles.S = S;
    handles.t = t;
    handles.f = f;
    handles.precomputed_spec = 1;
    handles.firstime = 0;
end

% S = log(S)';
% Smax = max(max(S));
% Smin = min(min(S));
% Ssmall = uint8(round(((S - Smin)/(Smax-Smin))*255));
% 
%  save('uint8_test.mat','Ssmall','-mat');
%  save('full_rest.mat','S','-mat');

handles.times=timeax(:);
handles.transition=[diff(boxcurve(:)); 0];

set( handles.axesS, 'XTick', [] );
set( handles.axesP, 'XTick', [] );
 
if exist('handles.datacursor')
    delete( handles.datacursor );
    delete( handles.segmentLineP );
    delete( handles.segmentLineS );
    delete( handles.segmentLineW );
end 

handles.datacursor=datacursormode(handles.figure1);  
axes(handles.axesP);
handles.segmentLineP = line('Visible','off');
axes(handles.axesS);
handles.segmentLineS = line('Visible','off');
axes(handles.axesW);
handles.segmentLineW = line('Visible','off');

if get( handles.SegmentButton, 'Value' )
    set(handles.datacursor,'Enable','on','DisplayStyle','datatip','SnapToDataVertex','off','UpdateFcn',@datacursorfunc);
end

guidata(gcbo,handles); 
close(hw);
handles = draw_segments(handles);

function  [Stot boxcurve] = compute_threshold_free(S,th,n)
% Computes the threshold based on a floating percentage of the maximum
% summed intensity
    Stot=sum(S,2); 
    boxcurve=Stot; 
    smax=max(Stot); 
    
    Stot = smooth_curve(Stot',n); % for removing extremes
    
    boxcurve(find(Stot<th*smax))= smax*th; 
    boxcurve(find(Stot>th*smax))= smax;
    
function [ratiof boxcurve] = compute_index(S,lower_range,upper_range,lowerfreq,upperfreq,indexthresh,n)
    % This algorithm is based on the method described in Aylin's
    % dissertation.
    
    S = S';
    nfreqs = length(S(:,1)); 
    freqspern = (upperfreq - lowerfreq) / nfreqs;
    
    indexinlower = fliplr(nfreqs - round((lower_range - lowerfreq)/freqspern));
    indexinupper = fliplr(nfreqs - round((upper_range - lowerfreq)/freqspern)) + 1;
    
    nrangelower = indexinlower(2)-indexinlower(1);
    nrangeupper = indexinupper(2)-indexinupper(1);
    
    ratiof = ( sum(S(indexinupper(1) : indexinupper(2),:)) / nrangeupper )... 
    ./ ( sum(S( indexinlower(1) : indexinlower(2),:)) / nrangelower );
     
    
    ratiof = smooth_curve(ratiof,n); % for smoothing the curve
    
    maxrf = max(ratiof);
    
    boxcurve = ratiof;
    
    boxcurve(find(ratiof<indexthresh))= indexthresh; 
    boxcurve(find(ratiof>=indexthresh))= maxrf;    
    
function smoothedcurve = smooth_curve(curve2smooth,n);
% Computes the moving average of the curve where n is an integer
% for example n = 1 averages the current point with the point before and afterwards

    m = length(curve2smooth);
    if m > 0
        curve2smooth = [repmat(curve2smooth(1),1,n) curve2smooth repmat(curve2smooth(m),1,n)];
        smoothedcurve = zeros(m,1);
    
        for i = 1:m
            smoothedcurve(i) = sum(curve2smooth(i:i + 2 * n)) / (2 * n + 1);
        end
    else % just to save computation time
        smoothed_curve = curve2smooth;
    end
    
function MinFreq_Callback(hObject, eventdata, handles)

function MinFreq_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function MaxFreq_Callback(hObject, eventdata, handles)
function MaxFreq_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in SpectrumType.
function SpectrumType_Callback(hObject, eventdata, handles)
function SpectrumType_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function AmpThresh_Callback(hObject, eventdata, handles)
function AmpThresh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function TDerThresh_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function TDerThresh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in PlayAll.
function PlayAll_Callback(hObject, eventdata, handles)
wavplay(handles.ts,handles.Fs);


% --- Executes on button press in PlayWindow.
function PlayWindow_Callback(hObject, eventdata, handles)
wavplay(handles.ts(handles.Tmin:handles.Tmax),handles.Fs,'async');

%h=handles.axesP; axes(h); semilogy(timeax,Stot); axis tight; 

function txt = datacursorfunc(empt,event_obj)
pos = get(event_obj,'Position');
handles=guidata(get(event_obj,'Target'));

  set(handles.segmentLineP,'Xdata',[pos(1) pos(1)],'Ydata',[0.00000000000001*pos(2) 1000000000000*pos(2)],'Visible','on' );
  set(handles.segmentLineS,'Xdata',[pos(1) pos(1)],'Ydata',[0.00000000000001*pos(2) 1000000000000*pos(2)],'Visible','on' );
  set(handles.segmentLineW,'Xdata',[pos(1) pos(1)],'Ydata',[-100000000000*pos(2) 1000000000000*pos(2)],'Visible','on' );

  if handles.start_stop_enable == 1
      set( handles.SegStartButton, 'Enable', 'on' );
  else
      set( handles.SegEndButton, 'Enable', 'on' );
  end

txt = {[num2str(pos(1))]};
guidata(gcbo,handles); 

function handles = draw_segments( handles )
n = 1;
while n <= length( handles.segments )
    handles.segments(n).lines=[];
    handles.segments(n) = draw_all_x( handles, handles.segments(n) );
    n = n + 1;
end
guidata(gcbo,handles); 

% --- Executes on button press in SegmentButton.
function SegmentButton_Callback(hObject, eventdata, handles)
toggled = get( handles.SegmentButton, 'Value' );

if toggled
    handles.segments = [];
    handles.segmentmode = 1;
    set( handles.SegmentButton, 'String', 'Segment On' );
    set( handles.SegmentButton, 'Enable', 'off' );
    if not(exist([handles.filename '.seg.txt']));
        set( handles.LoadSegments, 'Enable', 'off' );
    else
        set( handles.LoadSegments, 'Enable', 'on' );
    end
    set( handles.AutoSegmentFile, 'Enable','on');
    set( handles.AutoSegButton, 'Enable', 'on' );  
    set( handles.SegmentLengthEdit, 'Enable', 'on' );  
    set( handles.SegmentLengthText, 'Enable', 'on ' );  
    set( handles.SaveSegments, 'Enable', 'on' );  
    set( handles.DeleteSegment, 'Enable', 'on' );
    set( handles.DeleteAllButton, 'Enable', 'on' );
    set( handles.SegCancel, 'Enable', 'on' );  
    set( handles.PlotSegments, 'Enable', 'on' );
    set( handles.LoadFile, 'Enable', 'off' );
    set( handles.LoadNext, 'Enable', 'off' );
    handles.start_stop_enable = 1;
    set(handles.datacursor,'Enable','on','DisplayStyle','datatip','SnapToDataVertex','off','UpdateFcn',@datacursorfunc);
    fprintf( 'Segment mode on!\n' );
else
    handles.segmentmode = 0;
    set( handles.SegmentButton, 'String', 'Segment Off' );
    set( handles.AutoSegButton, 'Enable', 'off' );  
    set( handles.AutoSegmentFile, 'Enable','off');
    set( handles.SegmentLengthEdit, 'Enable', 'off' );  
    set( handles.SegmentLengthText, 'Enable', 'off' );  
    set( handles.LoadSegments, 'Enable', 'off' );
    set( handles.SaveSegments, 'Enable', 'off' );  
    set( handles.SegStartButton, 'Enable', 'off' );
    set( handles.SegEndButton, 'Enable', 'off' );
    set( handles.DeleteSegment, 'Enable', 'off' );
    set( handles.DeleteAllButton, 'Enable', 'off' );
    set( handles.SegCancel, 'Enable', 'off' );  
    set( handles.PlotSegments, 'Enable', 'off' );
    set( handles.LoadFile, 'Enable', 'on' );
    set( handles.LoadNext, 'Enable', 'on' );
    set(handles.datacursor,'Enable','off')
    fprintf( 'Segment mode off!\n' );
end
guidata(gcbo,handles); 


% --- Executes on button press in SegStartButton.
function SegStartButton_Callback(hObject, eventdata, handles)
set( handles.LoadSegments, 'Enable', 'off' );
set( handles.SegStartButton, 'Enable', 'off' );
handles.start_stop_enable = 0;
xy=get(handles.segmentLineP,'Xdata');
handles.segment.start=xy(1);
handles.segment.lines=[];
axes(handles.axesP);
set(handles.segmentLineP,'LineWidth',3);
handles.segment.lines(1) = handles.segmentLineP;
handles.segmentLineP = line('Visible','off');
axes(handles.axesS);
set(handles.segmentLineS,'LineWidth',3);
handles.segment.lines(2) = handles.segmentLineS;
handles.segmentLineS = line('Visible','off');
axes(handles.axesW);
set(handles.segmentLineW,'LineWidth',3);
handles.segment.lines(3) = handles.segmentLineW;
handles.segmentLineW = line('Visible','off');

guidata(gcbo,handles); 

% --- Executes on button press in SegEndButton.
function SegEndButton_Callback(hObject, eventdata, handles)
set( handles.SegEndButton, 'Enable', 'off' );
handles.start_stop_enable = 1;
xy=get(handles.segmentLineP,'Xdata');
handles.segment.end=xy(1);
handles.segment=draw_all_x( handles, handles.segment );
handles.segments = [handles.segments handles.segment];
guidata(gcbo,handles); 

function out=draw_all_x( handles, segment )
segment=draw_x( handles.axesP, segment );
segment=draw_x( handles.axesS, segment );
segment=draw_x( handles.axesW, segment );
out=segment;

function out=draw_x( theaxes, segment )
axes(theaxes);
ylim = get(theaxes,'YLim');
segment.lines = [segment.lines line('Xdata',[segment.start segment.start],'Ydata',ylim,'LineWidth',3)];
segment.lines = [segment.lines line('Xdata',[segment.end segment.end],'Ydata',ylim,'LineWidth',3)];
segment.lines = [segment.lines line('Xdata',[segment.start segment.end],'Ydata',ylim,'LineWidth',3)];
segment.lines = [segment.lines line('Xdata',[segment.start segment.end],'Ydata',[ylim(2) ylim(1)],'LineWidth',3)];
out=segment;

% --- Executes on button press in JumpBack.
function JumpBack_Callback(hObject, eventdata, handles)
Jump_shared(hObject, eventdata, handles, -1 )

% --- Executes on button press in Jump.
function Jump_Callback(hObject, eventdata, handles)
Jump_shared(hObject, eventdata, handles, 1 )

function Jump_shared(hObject, eventdata, handles, jump_dir )
Tim=eval(get(handles.DisplayWindow,'String'));
tDuration = str2num(get(handles.Duration,'String'));
maxTslider = (tDuration - Tim)/tDuration;
NT=min(round(Tim*handles.Fs),length(handles.ts));
Tslider=get(handles.slider1,'Value');
Tslider = Tslider + jump_dir * Tim * handles.Fs / length(handles.ts);

if Tim > tDuration
    set(handles.DisplayWindow,'String',num2str(tDuration));
    Tslider = 0;
end


if jump_dir == 1 % jumping forward
    if Tslider > maxTslider
        Tslider = maxTslider;
    end
end

if jump_dir == -1 % jumping backwards
    if Tslider < 0
        Tslider = 0
    end
end

    
% if Tslider > 1
%     Tslider = ( length(handles.ts) - NT ) / length(handles.ts);
% end
% if Tslider < 0
%     Tslider = 0
% end
set(handles.slider1,'Value',Tslider);
guidata(gcbo,handles); 
Plot_Callback(hObject, eventdata, handles)


function LoadNext_Callback(hObject, eventdata, handles)

% Get filename, extension.  Look for next file with same extension, no seg
% file associated

exclude_name = [handles.filename, get(handles.ExcludeExt,'String')];
if not(exist(exclude_name))
     fid=fopen( exclude_name, 'w' );
     fclose( fid);
end

[path] = cell2mat(regexp( handles.filename, '^.*\', 'match' ));
[extension] = cell2mat(regexp( handles.filename, '\w+$', 'match' ));
dirlist = dir( [path '*' extension] );
ndir = length(dirlist);
n = 1;
while n <= ndir
    file = dirlist(n).name;
    if not(exist([path file get(handles.ExcludeExt,'String')]))
        break;
    end
       n = n + 1;
end
if n <= ndir
    set( handles.FileNameString, 'String',file);
    handles.filename = [path file];
    guidata(gcbo,handles); 
    handles = loadfile(hObject, eventdata, handles);
else
    error('No more files found matching desired pattern');
end

% --- Executes on button press in Precompute.
function Precompute_Callback(hObject, eventdata, handles)
% handles = guidata(gcbo);
toggled = get( hObject, 'Value' );
if toggled
    
    % Disable spectra configuration parameters
    
%     set(handles.DisplayWindow, 'Enable', 'off');
    set(handles.WinSize, 'Enable', 'off');
    set(handles.StepSize, 'Enable', 'off');
    set(handles.TW, 'Enable', 'off');
    set(handles.MinFreq, 'Enable', 'off');
    set(handles.MaxFreq, 'Enable', 'off');
    set(handles.SpectrumType, 'Enable', 'off');
%     set(handles.AmpThresh, 'Enable', 'off');
%     set(handles.TDerThresh, 'Enable', 'off');
    set(handles.LoadNext, 'Enable','off');
%    set(handles.LoadFile, 'Enable','off');
    
    valueTslider = get(handles.slider1,'Value');
    set(handles.slider1,'Value',0);
    strDuration = get(handles.Duration,'String');
    strWindow = get(handles.DisplayWindow,'String');
    
    handles.firsttime = 1; % indicates that the spectra need to be calculated
    
    if str2num(strDuration) > handles.maxspec_t
        strDuration = num2str(handles.maxspec_t);
    end
    
    set(handles.DisplayWindow,'String',strDuration);
    
    Plot_Callback(handles.Plot, eventdata, handles);
    
    handles = guidata(hObject);
    
    handles.firsttime = 0;
    handles.precomputed_spec = 1;
    set(handles.DisplayWindow,'String',strWindow);
    set(handles.slider1,'Value',valueTslider);
    Plot_Callback(handles.Plot, eventdata, handles);
    handles = guidata(hObject);
    handles.precomputed_spec = 1;
else
  handles.precomputed_spec = 0;
  
  
  % Enable spectra configuration parameters
  
  handles.S = []; % release memory
  handles.t = [];
  handles.f = [];
  
  set(handles.WinSize, 'Enable', 'on');
    set(handles.StepSize, 'Enable', 'on');
    set(handles.TW, 'Enable', 'on');
    set(handles.MinFreq, 'Enable', 'on');
    set(handles.MaxFreq, 'Enable', 'on');
    set(handles.SpectrumType, 'Enable', 'on');
%     set(handles.AmpThresh, 'Enable', 'on');
    set(handles.TDerThresh, 'Enable', 'on');
    set(handles.LoadNext, 'Enable','on');
    set(handles.LoadFile, 'Enable','on');
  
end

guidata(hObject,handles);

function Precompute_CreateFcn(hObject, eventdata, handles)
   
function Path_Callback(hObject, eventdata, handles)
path=get(hObject,'String')


function Path_CreateFcn(hObject, eventdata, handles)
set(hObject,'String',pwd);
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Extensions_Callback(hObject, eventdata, handles)

function Extensions_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','wav');
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Duration_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function LoadSegments_Callback(hObject, eventdata, handles)

handles = load_segment(handles, [handles.filename '.seg.txt'] );
set( handles.LoadSegments, 'Enable', 'off' );
handles = draw_segments( handles );
guidata(gcbo,handles);

function handles=load_segment(handles,filename)
fid=fopen( filename, 'r' );
segments = [];
scanned=fscanf( fid, '%g %g',[2 inf] );
n = 1;

while n <= size(scanned, 2)
    segment.start = scanned(1,n);
    segment.end = scanned(2,n);
    segment.lines = [];
    segments = [ segments segment ];
    n = n + 1;
end

handles.allsegments = segments; % all segments holds all segments for the file
handles.segments = filtersegments(handles,handles.allsegments); % get segments for the current chunk
handles.loadedsegment = 1; % indicates segments have been filtered

guidata(gcf,handles); 


function filteredsegments = filtersegments(handles,segments)
% Returns segments which are in the current defined view. Returns segments
% which are not cut off.

realstart = handles.markerstart / handles.Fs;
realend = handles.markerend / handles.Fs;

filteredsegments = [];

for i = 1:length(segments) % no garuantee segments are in the same order
    if (segments(i).start >= realstart) && (segments(i).end <= realend)
        filteredsegments = [filteredsegments segments(i)];
    end
end
    
for i=1:length(filteredsegments)
    filteredsegments(i).start =  filteredsegments(i).start - realstart;
    filteredsegments(i).end =  filteredsegments(i).end - realstart;
end

function ExcludeExt_Callback(hObject, eventdata, handles)

% Hints: get(hObject,'String') returns contents of ExcludeExt as text
%        str2double(get(hObject,'String')) returns contents of ExcludeExt as a double

function ExcludeExt_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function DeleteSegment_Callback(hObject, eventdata, handles)


   pos=get(handles.segmentLineP,'Xdata');    
    n = 1;
    while n <= length( handles.segments )
        if pos(1) >= handles.segments(n).start && pos(1) <= handles.segments(n).end
            handles=delete_segment( handles, n );
        else
           n = n + 1;
        end 
    end
    drawnow;
guidata(gcbo,handles); 

function handles=delete_segment( handles, n )
            nl = 1;
            while nl <= length( handles.segments(n).lines )
                set( handles.segments(n).lines(nl), 'Visible', 'off'); 
                nl = nl + 1;
            end
            handles.segments(n) = [];
            fprintf('deleted!\n');

function SaveSegments_Callback(hObject, eventdata, handles)

% For the currently defined segments append to the segment list

handles = savesegments2mem(handles);

segment_file = fopen( [handles.filename '.seg.txt'], 'wt' );
n = 1;
while n <= size(handles.allsegments, 2)
    fprintf( segment_file, '%f %f\n', handles.allsegments(n).start, handles.allsegments(n).end );
    n = n + 1;
end
fclose(segment_file);
set( handles.SegmentButton, 'Enable', 'on' );
guidata(gcbo,handles); 

function handles=savesegments2mem(handles)
% Updates the handles allsegments in memory

% first remove in all segments all segments which are in the current chunk

oldsegments = [];

realstart = handles.markerstart / handles.Fs; % readjust time
realend = handles.markerend / handles.Fs; % readjust time

for i = 1:length(handles.allsegments)
    if not((handles.allsegments(i).start) >= realstart && (handles.allsegments(i).end <= realend))
        oldsegments = [oldsegments handles.allsegments(i)];
    end
end

% now put in the new segments

newsegments = [];

for i = 1:length(handles.segments)
    segment = handles.segments(i);
    segment.start = segment.start + realstart;
    segment.end = segment.end + realstart;
    newsegments = [newsegments segment];
end

handles.allsegments = [oldsegments newsegments];

function SegCancel_Callback(hObject, eventdata, handles)
set( handles.SegmentButton, 'Enable', 'on' );
guidata(gcbo,handles); 


function PlotSegments_Callback(hObject, eventdata, handles)

% Load Segments in directory

[path] = cell2mat(regexp( handles.filename, '^.*\', 'match' ));
[extension] = cell2mat(regexp( handles.filename, '\w+$', 'match' ));

path=get(handles.Path,'String');
extension=get(handles.Extensions,'String');
dirlist = dir( [path '\*' extension '.seg.txt'] );
ndir = length(dirlist);
n = 1;
all_segments = [];
while n <= ndir
    file = dirlist(n).name;
    segments = load_segment([path '\' file]);
    all_segments = [all_segments segments];
    n = n + 1;
end

% Plot info
if length(all_segments) > 2
    
    figure();
    axes();
    nbin= max(length([all_segments.end])/5,10);
    syllable_lengths=[all_segments.end]-[all_segments.start];
    hi=hist( syllable_lengths ,nbin);
    tl=min( syllable_lengths );
    th=max( syllable_lengths );
    times=tl:((th-tl)/(nbin-1)):th;
    plot(times,hi);
    xlabel('Segment Length (s)');
    ylabel('N');
    title(['All segments in ' path]);
else
    error('too few segments to plot');
end
guidata(gcbo,handles); 
  


function AutoSegButton_Callback(hObject, eventdata, handles)
    
    n = 1;
    segments = [];
    segment.start = 0;
    segment.end = 0;
    segment.lines = [];
    minlen = eval(get( handles.SegmentLengthEdit, 'String' ));
    while n < length( handles.times )
        
       if ( handles.transition(n) > 0 )
           segment.start = handles.times(n);
       end
       if ( handles.transition(n) < 0 )
           segment.end = handles.times(n);
       end
       if (segment.start > 0) && (segment.end) > 0 && (segment.end - segment.start) > minlen
           segments = [ segments segment ];
           segment.start = 0;
            segment.end = 0;
       end
       n = n + 1; 
    end
    
    handles.segments = [handles.segments segments]; 
    handles = draw_segments( handles );
    guidata(gcbo,handles);
    
function DeleteAllButton_Callback(hObject, eventdata, handles)

while length( handles.segments )
    handles = delete_segment( handles, 1 );
end
guidata(gcf,handles); 


% --- Executes on button press in PlotAllButton.
function PlotAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.slider1,'Value',0);

strDuration = get(handles.Duration,'String');

if str2num(strDuration) > handles.maxspec_t
    strDuration = num2str(handles.maxspec_t);
end


set(handles.DisplayWindow,'String',strDuration);
Plot_Callback(hObject, eventdata, handles);


% --- Executes on button press in PreviousChunk.
function PreviousChunk_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousChunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% [handles.markerstart handles.markerend]
handles = savesegments2mem(handles);
handles = loadfile(hObject, eventdata, handles,[handles.markerstart-handles.maxwavsize-1,handles.markerstart-1]);
% [handles.markerstart handles.markerend]
;
guidata(gcf,handles);

% --- Executes on button press in NextChunk.
function NextChunk_Callback(hObject, eventdata, handles)
% hObject    handle to NextChunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = savesegments2mem(handles);
handles = loadfile(hObject, eventdata, handles, [handles.markerend+1,handles.markerend+1+handles.maxwavsize]);
guidata(gcf,handles);

% --- Executes on button press in AutoSegmentFile.
function AutoSegmentFile_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSegmentFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

while (handles.markerend < handles.wavsize)

PlotAllButton_Callback(hObject, eventdata, handles);
handles = guidata(gcbo);
AutoSegButton_Callback(hObject, eventdata, handles);
handles = guidata(gcbo);
NextChunk_Callback(hObject, eventdata, handles);
handles = guidata(gcbo);

end


PlotAllButton_Callback(hObject, eventdata, handles);
handles = guidata(gcbo);
AutoSegButton_Callback(hObject, eventdata, handles);
handles = guidata(gcbo);
guidata(gcbo,handles);


function MaxSegLength_Callback(hObject, eventdata, handles)
% hObject    handle to MaxSegLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxSegLength as text
%        str2double(get(hObject,'String')) returns contents of MaxSegLength as a double


% --- Executes during object creation, after setting all properties.
function MaxSegLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxSegLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function MaximumWavSize_Callback(hObject, eventdata, handles)
% hObject    handle to MaximumWavSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaximumWavSize as text
%        str2double(get(hObject,'String')) returns contents of MaximumWavSize as a double


% --- Executes during object creation, after setting all properties.
function MaximumWavSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaximumWavSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in SeekButton.
function SeekButton_Callback(hObject, eventdata, handles)
% hObject    handle to SeekButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Seek anywhere in a long file

handles = savesegments2mem(handles);

try 
    timetoseek = str2num(get(handles.SeektoEdit,'String'));
catch
    timetoseek = 0;
    set(handles.SeektoEdit,'String','0');
end

if timetoseek < 0
    timetoseek = 0;
    set(handles.SeektoEdit,'String','0');
end

timetoseek = round(timetoseek * handles.Fs);

if timetoseek >= handles.wavsize
    timetoseek = timetoseek - handles.maxwavsize;
end

timetoseek = timetoseek  + 1;
timetoseekend = timetoseek + handles.maxwavsize;

if timetoseekend > handles.wavsize
    timetoseekend = handles.wavsize;
end

oldstate = handles.dontcutsegments;
handles.dontcutsegments = 0;
handles = loadfile(hObject,eventdata,handles,[timetoseek timetoseekend]);
handles.dontcutsegments = oldstate;

guidata(gcbo,handles);

function SeektoEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SeektoEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SeektoEdit as text
%        str2double(get(hObject,'String')) returns contents of SeektoEdit as a double


% --- Executes during object creation, after setting all properties.
function SeektoEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SeektoEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function RealDuration_Callback(hObject, eventdata, handles)
% hObject    handle to RealDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RealDuration as text
%        str2double(get(hObject,'String')) returns contents of RealDuration as a double


% --- Executes during object creation, after setting all properties.
function RealDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RealDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on selection change in AutoMethodPopupMenu.
function AutoMethodPopupMenu_Callback(hObject, eventdata, handles)
% hObject    handle to AutoMethodPopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns AutoMethodPopupMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AutoMethodPopupMenu
contents = get(hObject,'String');
method = contents{get(hObject,'Value')}

if strcmp(method,'Summed intensity')
    handles.automethod = 'threshold';
    set(handles.AmpThresh,'Visible','on');
    set(handles.RatioThresh,'Visible','off');
elseif strcmp(method,'Ratio')
    handles.automethod = 'ratiof';
    set(handles.AmpThresh,'Visible','off');
    set(handles.RatioThresh,'Visible','on');
end

guidata(gcbo,handles);

% --- Executes during object creation, after setting all properties.
function AutoMethodPopupMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AutoMethodPopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function RatioThresh_Callback(hObject, eventdata, handles)
% hObject    handle to RatioThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RatioThresh as text
%        str2double(get(hObject,'String')) returns contents of RatioThresh as a double


% --- Executes during object creation, after setting all properties.
function RatioThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RatioThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function RatioLower_Callback(hObject, eventdata, handles)
% hObject    handle to RatioLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RatioLower as text
%        str2double(get(hObject,'String')) returns contents of RatioLower as a double


% --- Executes during object creation, after setting all properties.
function RatioLower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RatioLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function RatioUpper_Callback(hObject, eventdata, handles)
% hObject    handle to RatioUpper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RatioUpper as text
%        str2double(get(hObject,'String')) returns contents of RatioUpper as a double


% --- Executes during object creation, after setting all properties.
function RatioUpper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RatioUpper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in OptionsDisplay.
function OptionsDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to OptionsDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OptionsDisplay

% positionP = get(handles.OptionsUiPanel,'Position');
% positionF = get(gcf,'Position');

state = get(hObject,'Value');

if state
    set(handles.OptionsUiPanel,'Visible','on');
%     positionF(3) = positionF(3) + positionP(3);
else
    set(handles.OptionsUiPanel,'Visible','off');
%     positionF(3) = positionF(3) - positionP(3);
end

%  set(gcf,'Position',positionF); % untested

guidata(gcbo,handles)



function Duration_Callback(hObject, eventdata, handles)
% hObject    handle to Duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Duration as text
%        str2double(get(hObject,'String')) returns contents of Duration as a double





function SmoothFactor_Callback(hObject, eventdata, handles)
% hObject    handle to SmoothFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SmoothFactor as text
%        str2double(get(hObject,'String')) returns contents of SmoothFactor as a double


% --- Executes during object creation, after setting all properties.
function SmoothFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SmoothFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function channel_Callback(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channel as text
%        str2double(get(hObject,'String')) returns contents of channel as a double



% --- Executes during object creation, after setting all properties.
function channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


