function varargout = classify_spectra(varargin)
% CLASSIFY_SPECTRA M-file for classify_spectra.fig
%      CLASSIFY_SPECTRA, by itself, creates a new CLASSIFY_SPECTRA or raises the existing
%      singleton*.
%
%      H = CLASSIFY_SPECTRA returns the handle to a new CLASSIFY_SPECTRA or
%      the handle to
%      the existing singleton*.
%
%      CLASSIFY_SPECTRA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLASSIFY_SPECTRA.M with the given input arguments.
%
%      CLASSIFY_SPECTRA('Property','Value',...) creates a new CLASSIFY_SPECTRA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before classify_spectra_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to classify_spectra_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help classify_spectra

% Last Modified by GUIDE v2.5 26-Jun-2006 23:19:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @classify_spectra_OpeningFcn, ...
                   'gui_OutputFcn',  @classify_spectra_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before classify_spectra is made visible.
function classify_spectra_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to classify_spectra (see VARARGIN)

% Choose default command line output for classify_spectra
handles.output = hObject;

% Set defaults
handles.recompute = logical(1); % Whether to recompute a spectra
handles.cwavfile = ''; % The current wave file

handles.directory = pwd;

handles.Fs = 44100; % Frequency of audio sampling per second
handles.movingwin=[0.01 0.002]; % Size of the moving window in seconds; the first number is the window size and the second is the step size
handles.tapers=[3 5]; 
handles.pad=1; 
handles.fpass=[0 20000]; % Range of frequency sampling 

handles.nsegments = 0; % total number of segments
handles.NextIndex = 1; % the index for segments
handles.maxseglength = 0; % set in seconds

% ClassifyAxes handles

handles.classified_height = 560 ; % the height of the image in the classified axes
handles.classified_width = 450;  % the width of the image in the classified axes

% 

handles.plotmode = 'spectra'; % The main spectra plot mode can also be
                              % see plot modes
                              
handles.plotmodes = {'spectra' 'waveform' 'spectra_dt' 'spectra_df' };

handles.plotmodevalue = 1;

% set up a density measurement which will allow scaling

classaxpos = get(handles.ClassifiedAxes,'Position');
handles.classified_height_density =  handles.classified_height / classaxpos(4);
handles.classified_width_density =   handles.classified_width / classaxpos(3);

handles.ispecheight = 100; % fixed height of the iconized spectogram
handles.ismaxwidth = SmallAxes_width(handles);
handles.specpad = 0.02; % pad in image sizes

handles.xspacer = 5; % fixed space between the images in the horizontal direction
handles.yspacer = 10; % fixed space between the row of the images

handles.xpsacer_density = handles.xspacer / classaxpos(3);
handles.ypsacer_density = handles.yspacer / classaxpos(4);

handles.image_list = {}; % holds an array of the specicons
handles.positions = [];  % holds the position of images on infinitely long canvas
handles.images_dim = []; % holds the size of the images

handles.mapindex = []; % holds the position number for the segment

handles.nimages = 0; % total number of images or spectra
handles.number_rows = floor(handles.classified_height/(handles.yspacer + handles.ispecheight)); % the number of rows allowed on a page 

handles.cnrows = 0; % current number of rows

handles.startpage = 1; % an index for the first image on the page
handles.endpage = 1; % an index for the last image on the page 

handles.startx = 1; % for the classified axes holds the start position for the image
handles.endx = handles.classified_width; % for the classified axes holds the end position for the image

handles.mode = 'browse';  % A string representing the current major mode which is either 
                          % 'browse','classify','class-view','class-members'

handles.submode = 'select'; %A string representing the minor mode which is either
                            % 'select','select-class','remove-class', 'compare', 
                            % 'typify'
 
handles.quickmode = logical(0); % Quick mode allows quick classification with minimum
                                % work for the user. By default this is set
                                % off

handles.lastsegment = 1; % The last segment classified 
                                
                                
handles.sortclass = 'length'; % Tells how classes are to be ordered in the ClassifiedAxes
                                % 'original' is the order the classes were created or loaded from
                                % 'popularity' sort the classes with most popular first
                                % 'length' sort the classes by longest
                                %  class first

set(handles.SortPopupMenu,'Value',3);
                               
handles.lastclass = 0; 

handles.lowerfreq = 0;  % Lower frequency for zooming
handles.upperfreq = 7500; % Upper frequency for zooming
handles.rezoom = logical(1);

handles.IconListf = {};

handles.baseclassname = 'mockingbird'; % This string should be set by the user
                                       % used as the base class name                                       
handles.nclasses = 0; % total number of syllable classes
handles.classes = []; % structure for holding class information
handles.current_class = 0; % used by compare to go through classes

handles.configschanged = logical(0); % indicates whether the configs for spectra has changed

handles.precomputed = logical(0);

handles.configfile = 'class_spec.conf';

handles.originalsize = [0 0 170 44]; % original position of the form
handles.originalaxessize = [80.5 6.375 86.167 34.438];
handles.blank = logical(1); % indicate that the ClassifiedAxes is blank
handles.prevsize = handles.originalsize;

handles.fixed = logical(1); % Whether to use fixed scaling when resizing the form
% initially set to true so not to call the repositioning algorithm when the
% form is blank
handles.nfeatures = 10;
handles.ncepestral = 10; % Number of cepestral coefficients to include

set(gcf, 'ResizeFcn', {@ResizeFcn});
%classify_spectra('ResizeFcn',gcbo,[],guidata(gcbo))
% Update handles structure
guidata(hObject, handles);
 
% UIWAIT makes classify_spectra wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = classify_spectra_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function DirectoryEditBox_Callback(hObject, eventdata, handles)
handles.directory = get(hObject,'String');

function DirectoryEditBox_CreateFcn(hObject, eventdata, handles)
set(hObject,'string',pwd);
guidata(gcbo,handles); 
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Functions for loading segments  %    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segments=LoadSegmentsFromDirectory(directory)
segfilelist = dir( [directory '\*' '.seg.txt'] );
%segfilelist;
nsegfile = length(segfilelist);
segments = [];
n = 1;
while n <= nsegfile
    segfilename = segfilelist(n).name;
    fid=fopen( segfilename, 'rt' );
    if fid ~= -1 % if there is a seg file
        scanned=fscanf( fid, '%g %g',[2 inf] );
        fclose(fid);
        %fprintf( 'File %d of %d has %d segments: %s\n', n, nsegfile, size(scanned,2),segfilename );
        wavfile = segfilename(1:(length(segfilename)-8)); % can this be made more general
        i = 1;
        while i <= size(scanned, 2) % Load the start and stop of segments
            segment.wavfile = wavfile;
            segment.class = ''; % Loaded segments start out unclassified
            segment.features = [];
            segment.start = scanned(1,i);
            segment.end = scanned(2,i);
            segment.specfilename = [segment.wavfile '.' num2str(segment.start) '-' num2str(segment.end) '.spec'];
            segments = [ segments segment];
            i = i + 1;
        end 
    end
    n = n + 1;
end

function handles=Load_InitialSegments(hObject,handles)
  %handles.directory = pwd;
  handles.segments = LoadSegmentsFromDirectory(handles.directory);
  handles.nsegments = length(handles.segments);
  guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Functions for handling syllable classes  %    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = cindex2imagelist(handles)
% Takes a cindex and generates a a list of images
handles.image_list = {};
%load('-mat','specicons');
for i = 1:length(handles.mapindex)
    handles.image_list{i} = handles.classes( handles.mapindex(i) ).iconS; 
end

function mapindex = sortindexbypop(handles)
% Takes a cindex and sorts by the class popularity
nindexes = length(handles.mapindex);
mapindex = [];
for i = 1:nindexes
    mapindex(i,1:2) = [handles.mapindex(i) handles.classes(handles.mapindex(i)).nmembers];
end
mapindex = sortrows(mapindex, 2);
mapindex = flipud(mapindex(:,1));

function mapindex = sortindexbylength(handles)
% Takes a cindex and sorts by the class popularity
nindexes = length(handles.mapindex);
mapindex = [];
for i = 1:nindexes
    mapindex(i,1:2) = [handles.mapindex(i) handles.classes(handles.mapindex(i)).length];
end
mapindex = sortrows(mapindex, 2);
mapindex = flipud(mapindex(:,1));

function class_string = newclassname(handles)
% Generates a new name for the class string using the baseclassname
% variable. Classes are numbered sequentially from the class with the largest number.

nbaseclass = length(handles.baseclassname);
classnum = 0;
for i = 1:handles.nclasses % make sure the largest class number is gotten
    classname = handles.classes(i).name;
    curr_classnum = str2num(classname(nbaseclass + 1:length(classname)));
    if curr_classnum > classnum
        classnum = curr_classnum;
    end
end
class_string = strcat(handles.baseclassname,num2str(classnum + 1));
;

function cindex = returnclassindex(handles,classname)
% Return the index to the class
    cindex = 0;
    i = 1;
    while (i <= handles.nclasses) && not(strcmp(classname,handles.classes(i).name))
        i = i + 1;
    end
    cindex = i;
;

function handles = add_new_class(handles,segment)
% Segment is the class that will be used to typify the class
handles.nclasses = length(handles.classes);
class.name = newclassname(handles);
class.nmembers = 1; % the number of segments which are members of this class
class.specfilename = segment.specfilename; % specfilename will be used as a unique identifier
class.index = handles.NextIndex;
%load('-mat',class.specfilename);
class.iconS = handles.IconList{handles.NextIndex}; % this is the icon which typifies the class
class.length = segment.end - segment.start; % used to hold the lengt of the length

handles.classes = [handles.classes class];
handles.nclasses = handles.nclasses + 1;
%guidata(gcbo,handles);
;

function handles=ConfigureClassSegment(handles)
% Handles the gui configuration of the class information when navigating
segment = handles.segments(handles.NextIndex);
if strcmp(segment.class,'') % Unclassified segment
    set(handles.ClassifyButton,'String','Classify');
    set(handles.ClassifyButton,'Enable','on');
else % Classified segment
    set(handles.ClassifyButton,'String','Declassify');
    set(handles.ClassifyButton,'Enable','on');
end
%guidata(gcbo,handles);

function handles = blankaxes(handles)
    set(handles.NextRowButton,'Enable','off');
    set(handles.PreviousRowButton,'Enable','off');
    axes(handles.ClassifiedAxes);
    handles.hiclass = image(uint8(zeros(handles.classified_height,handles.classified_width)));
    set(handles.ClassifiedAxes,'Xtick',[]);
    set(handles.ClassifiedAxes,'Ytick',[]);
    handles.blank = logical(1);

% --- Executes on button press in ClassifyButton.
function ClassifyButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
status = get(handles.ClassifyButton,'String');
set(handles.TypifyClassButton,'Visible', 'off');
set(handles.RemoveClassButton, 'Visible', 'off');
set(handles.NextClassButton,'Visible','off');
% set(handles.CompareToggleButton,'Visible','off');
set(handles.RenameClassButton,'Visible','off');
handles.startx = 1; % if you have been scrolling reset to defaults
handles.endx = handles.classified_width;
set(handles.SliderClassified,'Value',0);

if strcmp(status,'Classify')
    set_status(handles,'');
    set(handles.ModePopupMenu,'Enable','off');
    set(handles.NextSpectraButton,'Enable','off');
    set(handles.PreviousSpectraButton,'Enable','off');
    set(handles.ClassifyButton,'Enable','off');
    set(handles.CompareToggleButton,'Enable','off');
    set(handles.QuickModeButton,'Enable','off');
    set(handles.AutoClassifyButton,'Enable','off');
    set(handles.NewClassButton,'Visible','on');
    if handles.nclasses > 0 % Make sure there is at least one class
        if strcmp(handles.mode,'comparison') % if you are in comparison mode
            set(handles.NewClassButton,'Visible','off');
            handles.segments(handles.NextIndex).class = handles.classes(handles.lastclass).name;
            set(handles.ModePopupMenu,'Enable','on');
            set(handles.ClassifyButton,'String','Declassify');
            set(handles.ClassifyButton,'Enable','on');
            set(handles.CompareToggleButton,'value',0);
            set(handles.CompareToggleButton,'Enable','on');
            handles = configureclassview(handles,'select-class');
            set_status(handles, ['Viewing all ', num2str(handles.nclasses),' classes']);
            set(handles.RemoveClassButton,'Enable','on');
            set(handles.RemoveClassButton,'Visible','on');
            set(handles.QuickModeButton,'Enable','on');
            set(handles.AutoClassifyButton,'Enable','on');
            setnavigationbuttons(handles);
        else % regular mode
            set(handles.SortText,'Visible','on');
            set(handles.SortPopupMenu,'Visible','on');
            handles = configureclassview(handles,'select-class');
            set_status(handles, ['Select a class']);
        end
    else % draw a blank image
        handles = blankaxes(handles);
        handles.mode = 'class-view';
        handles.submode = 'select-class';
        set(handles.hiclass,'ButtonDownFcn',{@DummyClassifyAxesClickCallBack,handles});
    end
    % configure the remaining gui
    handles = SetModePopupMenu(handles,'class view');
elseif strcmp(status,'Declassify')
    cindex = returnclassindex(handles,handles.segments(handles.NextIndex).class);
    if handles.classes(cindex).nmembers == 1 % Only one member left of that class
        if length(handles.classes) == 1 % Only one class remaining
            axes(handles.ClassifiedAxes);
            handles.classes = [];
            handles.nclasses = 0;
            handles = blankaxes(handles);
        else % remove the class
            handles.classes = [handles.classes(1:cindex-1) handles.classes(cindex + 1:length(handles.classes))];
            handles.nclasses = handles.nclasses - 1;
        end     
    else 
        if strcmp(handles.classes(cindex).specfilename,handles.segments(handles.NextIndex).specfilename)
            i = 1; %Test if the class you are removing the type class
            segs = [ handles.segments(1:(handles.NextIndex - 1)) handles.segments((handles.NextIndex + 1) : handles.nclasses)];
            while (i <= length(segs)) && strcmp(segs(i).class,handles.classes(cindex).name)
                i = i + 1;
            end
            handles.classes(cindex).specfilename = handles.segments(i).specfilename;
            handles.classes(cindex).length = handles.segments(i).end - handles.segments(i).start;
            handles.iconS = handles.IconList{i};
        end   
        handles.classes(cindex).nmembers = handles.classes(cindex).nmembers - 1;
    end
    handles.segments(handles.NextIndex).class = ''; % Remove class information
    if handles.nclasses >= 1 % Redraw axes
        if strcmp(handles.mode,'class-view')
            handles = configureclassview(handles,'select');
            set_status(handles,['Viewing all ' num2str(handles.nclasses)  ' classes']);
        elseif strcmp(handles.mode,'class-members')
            handles = configureclassmembers(handles,handles.classes(cindex).name);
            set_status(handles,['Viewing ' num2str(handles.classes(cindex).nmembers) ' members of ' num2str(handles.classes(cindex).name)]);
        elseif strcmp(handles.mode,'browse')
            ; % do nothing
        end
    end
    set(handles.ClassifyButton,'String','Classify');
    
    setnavigationbuttons(handles);
end
guidata(gcbo,handles);


% --- Executes on button press in NewClassButton.
function NewClassButton_Callback(hObject, eventdata, handles)
% hObject    handle to NewClassButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.startx = 1; % if you have been scrolling reset to defaults
handles.endx = handles.classified_width;
set(handles.SliderClassified,'Value',0);

handles = add_new_class(handles,handles.segments(handles.NextIndex));
handles.segments(handles.NextIndex).class = handles.classes(handles.nclasses).name;

set(handles.SortText,'Visible','off');
set(handles.SortPopupMenu,'Visible','off');
set(handles.ClassifyButton,'Enable','on');
set(handles.ClassifyButton,'String','Declassify');
set(handles.NewClassButton,'Visible','off');
set(handles.ModePopupMenu,'Enable','on');
set(handles.NextSpectraButton,'Enable','on');
set(handles.CompareToggleButton,'Enable','on');
set(handles.QuickModeButton,'Enable','on');
set(handles.AutoClassifyButton,'Enable','on');
set(handles.PreviousSpectraButton,'Enable','on');
handles = configureclassview(handles,'xxx');
setnavigationbuttons(handles);
set_status(handles,'');
handles = SetModePopupMenu(handles,'class view');
guidata(gcbo,handles);
;

% --- Executes on button press in RemoveClassButton.
function RemoveClassButton_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveClassButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set(handles.TypifyClassButton,'Enable','on');
set(handles.RemoveClassButton,'Enable','off');
handles.mode = 'class-view';
handles.submode = 'remove-class';
set_status(handles,'Select a class to remove');
guidata(gcbo,handles);


% --- Executes on button press in TypifyClassButton.
function TypifyClassButton_Callback(hObject, eventdata, handles)
% hObject    handle to TypifyClassButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.TypifyClassButton,'Enable','off');
%set(handles.RemoveClassButton,'Enable','on');
handles.mode = 'class-members';
handles.submode = 'typify';
set_status(handles,'Select an icon to change the type');
guidata(gcbo,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Functions for computing the spectra    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segment=precompute_spectra(handles, segment)    
% This function will do the precomputing of the spectragram
spec_ipad = round(handles.specpad * handles.Fs);
spec_istart = round(segment.start * handles.Fs);
spec_iend = round(segment.end * handles.Fs);

% This is to catch an over and under run errors in the wav file because of the padding
try 
    [data] = wavread(segment.wavfile, [spec_istart - spec_ipad, spec_iend + spec_ipad]);
catch
    errmsg = lasterr;
    if strfind(errmsg, 'Sample limits out of range')
        if (segment.start - handles.specpad) < 0 % Make sure the starting point is not negative
            [data] = wavread(segment.wavfile, [1 spec_iend + spec_ipad]);
        else % over run of the buffer
            [data] = wavread(segment.wavfile);
            [data] = data((spec_istart - spec_ipad):length(data));
        end
    end
end

[Sfull tfull f] = compute_spectra(data,handles.tapers,handles.Fs,handles.fpass,handles.movingwin); % precompute the portion

wavlength = length(data);
Ssize = size(Sfull);

Slength = Ssize(2);
RatioWS = Slength / (wavlength / handles.Fs); % this allows us to index by time through spec file

Sstart = round(RatioWS * segment.start);
Send = round(RatioWS * segment.end);
Spad = round(RatioWS * handles.specpad);

Spre = Sfull(:,1:Spad);
S = Sfull(:,Spad+1:Spad + Send-Sstart);
Spost = Sfull(:,Spad + (Send-Sstart)+1:Slength);
t=[segment.start, segment.end];

iconS = iconify_spec(S,handles.ispecheight);

save(segment.specfilename,'S','t','f','Spre','Spost','RatioWS','tfull','iconS','-mat');

fprintf('Saving %s file\n',segment.specfilename);

function handles = precompute_AllSpectra(handles)
% This function precomputes all the spectra in a directory
hw = waitbar(0,'Precomputing spectra. . .');
if handles.nsegments >= 1
    for i = 1:handles.nsegments
        precompute_spectra(handles,handles.segments(i));
        waitbar(i/handles.nsegments);
    end
end
close(hw);
handles.IconList = get_SpecIcons(handles);
;
function [S t f]=compute_spectra(data,tapers,Fs,fpass,movingwin)
    data = data / std(data); % normalize the variance of the spectra
    params.tapers=tapers; params.Fs=Fs; params.fpass=fpass;
    [S t f] = mtspecgramc( diff(data), movingwin, params );
    S = log(S)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for plotting the spectragram   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles=get_and_plot(handles, segment)
load('-mat',segment.specfilename);
axes(handles.ToClassifyAxes);
% RatioWS

cmap = jet(256);

if strcmp(handles.plotmode,'spectra')
    SFull = cat(2,Spre,S,Spost);
elseif strcmp(handles.plotmode,'spectra_dt') || strcmp(handles.plotmode,'spectra_df') || strcmp(handles.plotmode,'waveform')
    wav_ipad = round(handles.specpad * handles.Fs);
    wav_istart = round(segment.start * handles.Fs);
    wav_iend = round(segment.end * handles.Fs);

    % This is to catch an over and under run errors in the wav file because of the padding
    try 
        [data] = wavread(segment.wavfile, [wav_istart - wav_ipad, wav_iend + wav_ipad]);
    catch
        errmsg = lasterr;
        if strfind(errmsg, 'Sample limits out of range')
            if (segment.start - handles.specpad) < 0 % Make sure the starting point is not negative
                [data] = wavread(segment.wavfile, [1 wav_iend + wav_ipad]);
            else % over run of the buffer
                [data] = wavread(segment.wavfile);
                [data] = data((wav_istart - wav_ipad):length(data));
            end
        end
    end
    
    data = data / std(data);
    
    params.Fs=handles.Fs;
    params.tapers = handles.tapers;
    params.fpass=handles.fpass;
    params.pad = 1;
    
    if strcmp(handles.plotmode,'spectra_dt')
         cmap =  gray(256);
        [SFull t f]= mtdspecgramc(diff(data),handles.movingwin,0,params); SFull=SFull';
    elseif strcmp(handles.plotmode,'spectra_df')
          cmap =  gray(256);
        [SFull t f]= mtdspecgramc(diff(data),handles.movingwin,pi/2,params); SFull=SFull';
    end
end

if strcmp(handles.plotmode,'spectra') || strcmp(handles.plotmode,'spectra_dt') || strcmp(handles.plotmode,'spectra_df')

    cmap(1,:) = [1, 1, 1];
    colormap(cmap);
    
    SFmin = min(min(SFull));
    SFmax = max(max(SFull));
    SFull = uint8(1 + round(255 * (SFull-SFmin) / (SFmax-SFmin)));

    hi = image(tfull + segment.start - handles.specpad,f,SFull);
    set(hi,'ButtonDownFcn',{@PlotModeCallBack});
    
    axis xy;
    
    hline1 = line([segment.start segment.start],[f(1) max(f)],'Color',[0 0 0],'LineWidth',3);
    hline2 = line([segment.end segment.end],[f(1),max(f)],'Color',[0 0 0],'LineWidth',3);
    
else

% xlim([tfull(1) tfull(length(tfull))] + segment.start - handles.specpad);

hp = plot(segment.start - handles.specpad + [0:length(data)-1] / handles.Fs, data);
set(handles.ToClassifyAxes,'YLim',[-5 5]);
set(hp,'ButtonDownFcn',{@PlotModeCallBack});
axis tight;

dataspan = [min(data) max(data)];

hline1 = line([segment.start segment.start],dataspan,'Color',[0 0 0],'LineWidth',3);
hline2 = line([segment.end segment.end],dataspan,'Color',[0 0 0],'LineWidth',3);

end

axes(handles.ToClassifySmallAxes);
ispecFull = uint8(zeros(handles.ispecheight,handles.ismaxwidth));

ispecFull = copy_into(ispecFull,handles.IconList{handles.NextIndex},1,1);

if length(ispecFull(1,:)) > handles.ismaxwidth
    ispecFull = ispecFull(:,1:handles.ismaxwidth);
end

if strcmp(get(handles.ZoomButton,'String'),'Zoom out')
    f = [handles.lowerfreq handles.upperfreq];
end

tsmall = [handles.movingwin(1),handles.ismaxwidth * handles.movingwin(2) - handles.movingwin(1)];
% 
% cmap = jet(256);
% cmap(1,:) = [1, 1, 1];
% 
% colormap(cmap);

% [0 (handles.ismaxwidth * (segment.end - segment.start))/length(iconS(1,:))]
ih = image(tsmall,f,flipud(ispecFull));
axis xy;
;

function PlotModeCallBack(src,eventdata)
% A Function for handling clicks to the axes
    handles = guidata(gcbo);
    
    handles.plotmodevalue = handles.plotmodevalue + 1;
    
    if handles.plotmodevalue > length(handles.plotmodes)
        handles.plotmodevalue = 1;
    end
    
    handles.plotmode = handles.plotmodes{handles.plotmodevalue};
    handles=get_and_plot(handles, handles.segments(handles.NextIndex));
    
    guidata(gcf,handles);
    
function handles=ConfigureSpecPlot(handles)
% Handles the gui configuration of the plotting
segment = handles.segments(handles.NextIndex);
if not(exist(segment.specfilename)) || handles.recompute
    precompute_spectra(handles,segment);
end

set(handles.ToClassifyPanel,'Title',['Segment ' num2str(handles.NextIndex) '/' num2str(handles.nsegments)])

segmentstatus = ['File: "' segment.wavfile '"; Segment length ' num2str(segment.end - segment.start,3)];

set(handles.SegmentText,'String',segmentstatus);
handles = get_and_plot(handles, segment);
guidata(gcbo,handles); 
;

function NextSpectraButton_Callback(hObject, eventdata, handles)
% Moves the segment viewer forward one segment
handles.NextIndex = handles.NextIndex + 1;
if handles.NextIndex == handles.nsegments
    set(handles.NextSpectraButton,'Enable','off');
end

if handles.NextIndex > 1
    set(handles.PreviousSpectraButton,'Enable','on');
end

handles=ConfigureClassSegment(handles);
handles=ConfigureSpecPlot(handles);
guidata(gcbo,handles); 
;

function PreviousSpectraButton_Callback(hObject, eventdata, handles)
% Moves the segment viewer backwards one segment
handles.NextIndex = handles.NextIndex - 1;

if handles.NextIndex == 1
    set(handles.PreviousSpectraButton,'Enable','off');
    set(handles.NextSpectraButton,'Enable','on');
end

if handles.NextIndex < handles.nsegments
    set(handles.NextSpectraButton,'Enable','on');
end

handles=ConfigureClassSegment(handles);
handles=ConfigureSpecPlot(handles);

set(handles.NextSpectraButton,'Enable','off');
set(handles.NextSpectraButton,'Enable','on');
guidata(gcbo,handles); 
;

function PrecomputeButton_Callback(hObject, eventdata, handles)
% Call back for the precompute button. This acts to load the file from the directory  

handles.directory = get(handles.DirectoryEditBox,'String');    

% fprintf('creating syllable list\n');
%handles.segments = segments;
handles.classes = [];
set(handles.PrecomputeButton, 'Enable', 'off' );
handles.NextIndex = 1;

handles = Load_InitialSegments(hObject,handles);
handles.recompute = logical(0);

if handles.nsegments >= 1
    handles=ConfigureClassSegment(handles);
    
    handles = load_configuration(handles,handles.configfile);
    
    if exist([handles.baseclassname '.dat']) && not(handles.configschanged)
        load('-mat','specicons');
        handles.IconList  = IconList;
        data = read_syllable_database(handles);
        handles = merge_syllable_database(handles,data);
    else
        handles = precompute_AllSpectra(handles);
    end
    
    %set(handles.ConfigureButton, 'Enable','off');
    handles=ConfigureSpecPlot(handles);
    handles=BrowseDirectory(handles);
    handles.precomputed = logical(1);
    set(handles.PlaySegmentButton, 'Enable', 'on' );
    set(handles.NextSpectraButton, 'Enable', 'on' );
    set(handles.ModePopupMenu,'Enable','on');
    set(handles.SaveButton,'Enable','on');
    set(handles.SaveItem,'Enable','on');
    set(handles.PrecomputeButton, 'Enable', 'off' );
    set(handles.CleanDirectoryItem, 'Enable', 'off' );
    set(handles.LoadDirectoryButton,'Enable','off');
    set(handles.LoadItem,'Enable','off');
    set(handles.ZoomButton,'Enable','on');
    set(handles.QuickModeButton,'Enable','on');
    set(handles.CompareToggleButton,'Enable','on');
    set(handles.ConfigureButton,'Enable','on');
    set(handles.ConfigureItem,'Enable','on');
    set(handles.AutoClassifyButton,'Enable','on');
    
    handles.fixed = logical(0); % turn off fixed scaling
    
    if not(strcmp(handles.segments(handles.NextIndex).class,'')) % Set the classification status
        set(handles.ClassifyButton,'String','Declassify');
    end
    
else
    set(handles.PrecomputeButton, 'Enable', 'on' );
end
guidata(gcbo,handles); 
;

% --- Executes on button press in PlaySegmentButton.

function PlaySegmentButton_Callback(hObject, eventdata, handles)
% Plays the current segment in the segment viewer
% hObject    handle to PlaySegmentButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
segment = handles.segments(handles.NextIndex);
data = wavread(segment.wavfile,[round(handles.Fs * segment.start),round(handles.Fs * segment.end)]);
wavplay(data,handles.Fs,'async');


function CurrentFilenameEdit_Callback(hObject, eventdata, handles)

function CurrentFilenameEdit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Functions for viewing spectra icon  %    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newwidth = SmallAxes_width(handles)
position1 = get(handles.ToClassifySmallAxes,'Position');
position2 = get(handles.ClassifiedAxes,'Position');

newwidth =  round(position1(3) *(handles.classified_width / position2(3))); 

function iconS = iconify_spec(S,height)
% Take a large spectra with high frequency bandwidth and reduce the height
% by pixel averaging
    Ssize = size(S);
    iconS = zeros(height,Ssize(2));
    
    % averaging of values to reduce size
    rf = floor(Ssize(1)/height);
    
    for i = 1:(height-1)
        for j = 1 : Ssize(2)
            iconS(i,j) = sum(S(((i-1)*rf)+1:i*rf,j))/rf;
        end
    end
    
    for j = 1 : Ssize(2) % take care of the last row by also pixel averaging
        iconS(height,j) = mean(S(rf*(height-1) : Ssize(1),j)); 
    end
    iconS = flipud(iconS);
    
    % Rescaling of values
    maxintense = max(max(iconS));
    minintense = min(min(iconS));
    
    iconS=uint8(1 + round(255 * (iconS-minintense)/(maxintense-minintense)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Functions for loading in images %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function destination = copy_into(destination,source,r,c)
% This copies an image into another image
% Two problems: one it can be optimized by removing the nested for loops
% two there is an indexing bug which returns a one pixel larger image
sized = size(destination);
sizes = size(source);

destination(r:r+sizes(1)-1,c:c+sizes(2)-1) = source(:,:);

% for i = 1:sizes(1)
%     for j = 1:sizes(2)
%         destination(r+i,c+j) = source(i,j);
%     end
% end
% ;

function positions = position_images(height,width,images_dim,xspacer,yspacer)
% This function returns a matrix consisting of two rows with the xy
% position for images.
% The function also assumes that the image's height is not restricted this
% allows for easier scrolling
% The function assumes that all images are of the same height
% height in pixels of the original
% width in pixels of the original
% images 
% xspacer in pixels for the horizontal space between images
% yspacer is the next height of the image
% image height is the fixed height of the images

number_images = length(images_dim(:,1));
imageheight = images_dim(1);

currentx = xspacer;
currenty = yspacer;

positions = zeros(number_images,2);

for i = 1:number_images 
    if (currentx + images_dim(i,2) + xspacer) > width % start a new row
        currentx = xspacer;
        currenty = currenty + imageheight + yspacer;
        positions(i,:) = [currenty currentx];
        currentx = currentx + images_dim(i,2) + xspacer;
    else
        positions(i,:) = [currenty currentx];
        currentx = currentx + images_dim(i,2) + xspacer; 
    end
    %positions(i,:) = [currenty currentx];
end

% function image_matrix = place_images_into(image_matrix, image_list, position_list)
% % Place images into a matrix
% number_images = length(image_list);
% 
% for i = 1:number_images
%     theimage = image_list{i};
%     image_matrix  = copy_into(image_matrix, theimage, position_list(i,1), position_list(i,2));
% end 
% ;

function image_dim = get_image_sizes(images)
% Returns an array of image sizes
number_images = length(images);
image_dim = zeros(number_images,2);
for i = 1:number_images
    image_dim(i,:) = size(images{i});
end
;

function image_matrix = place_images_into(image_matrix, image_list, position_list)
% Place images into a matrix
number_images = length(image_list);

for i = 1:number_images
    theimage = image_list{i};
    image_matrix  = copy_into(image_matrix, theimage, position_list(i,1), position_list(i,2));
end 
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for manipulating and plotting spectra icons    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IconList=get_SpecIcons(handles)
IconList = {};
status = 0;
for i = 1:handles.nsegments
    load('-mat',handles.segments(i).specfilename);
    IconList{i} = iconS;
end
save('specicons','IconList','-mat'); % temp code for saving


function handles = plot_classified_axes(handles, image_list, position_list)
% Low level drawing of the classified axes
handles.blank = logical(0);
handles.classmatrix = uint8(zeros(handles.classified_height,handles.classified_width));
axes(handles.ClassifiedAxes);
handles.classmatrix = place_images_into(handles.classmatrix,image_list,position_list);
set(handles.ClassifiedAxes,'XTick',[]);
set(handles.ClassifiedAxes,'YTick',[]);
handles.max_width = length(handles.classmatrix(1,:));

if  handles.max_width > handles.classified_width
    set(handles.SliderClassified,'Min',0);
    set(handles.SliderClassified,'Max',handles.max_width - handles.classified_width);
    set(handles.SliderClassified,'enable','on');
    handles.endx = handles.classified_width;
else
    set(handles.SliderClassified,'enable','off')
    handles.startx = 1;
    handles.endx = handles.classified_width;
end

classview = handles.classmatrix(:,handles.startx:handles.endx); % will cut overhang
handles.hiclass = image(classview);
set(handles.hiclass,'ButtonDownFcn',{@ClassifyAxesClickCallBack});
set(handles.ClassifiedAxes,'XTick',[]);
set(handles.ClassifiedAxes,'YTick',[]);

setrowbuttons(handles);

function handles = reposition_images(handles, image_list)
% this is a lower level function which is called to reposition the images.
% this would be called from higher level functions when images are added,
% deleted, or a new list of images needs to be loaded.

% initialize the handles for the images
handles.nimages = length(image_list);
handles.images_dim = get_image_sizes(image_list);

handles.positions = position_images(handles.classified_height,handles.classified_width,handles.images_dim,handles.xspacer,handles.yspacer);
handles.cnrows = length(unique(handles.positions(:,1)));

% Setup the first page view
handles.number_rows = floor(handles.classified_height / (handles.ispecheight + handles.yspacer));
handles.startpage = 1;
handles.endpage = 0;

for i = 1 : handles.number_rows
    handles.endpage = next_row_end(handles.positions,handles.endpage);
end
    
%guidata(gcbo,handles);
;

function nrow = which_row(positions,index)

nrow = 1;
i = 2;
while i <= index
    if not(positions(i,1) == positions(i-1,1))
        nrow = nrow + 1;
    end
    i = i + 1;
end

function cpositions = get_curr_position(handles)
% Setups the current view of the positions
    cpositions = handles.positions(handles.startpage:handles.endpage,:); % get the current view
    cpositions(:,1) = cpositions(:,1) - cpositions(1,1) + handles.yspacer;
    cpositions(:,2) = cpositions(:,2) - (handles.startx - 1);
;
    
function cposition = next_row_start(positions,cposition)
% Computes the position of the next row if the row based on the positions matrix
% it computes the position where the row starts
    npos = length(positions);
    i = cposition;
    while  (i <= npos) && positions(i,1) == positions(cposition,1)
        i = i + 1;
    end
    if i < npos % make sure the row numbers match
        cposition = i;
    end
    ;
    
function cposition = next_row_end(positions,cposition)
% Computes the position of the next row if the row based on the positions matrix
% it computes the last position before a new row starts
    npos = length(positions(:,1));
    if cposition < npos % not at the last row
        i = cposition + 1;
        while (i < npos) && (positions(i,1) == positions(cposition+1,1))
            i = i + 1;
        end
        if (positions(cposition + 1) == positions(npos))
            cposition = npos; % handle the condition that you are now at the last row 
        else
            cposition = i - 1; % make sure the row numbers match
        end
    else % handles the condition you are already at the last row
        cposition = npos;
    end
    ;
    
    
function handles = row_forward(handles)
% Moves the row forward in the classifiedaxes/browser view
    startpage = next_row_start(handles.positions,handles.startpage);
    endpage = next_row_end(handles.positions,handles.endpage);
    
    if not(endpage == handles.endpage) % indicates you are not at the last page
        handles.startpage = startpage;
        handles.endpage = endpage;
    end
%     guidata(gcbo,handles);
  ;

function cposition = previous_row_start(positions,cposition)
% Computes the position of the previous row if the row based on the positions matrix
% it computes the position where the row starts
    npos = length(positions(:,1));
    if cposition > 1
        i = cposition - 1;
        while (i > 1) && (positions(i,1) == positions(cposition-1,1))
            i = i - 1;
        end
        if (i > 1) && (cposition ~= 2)
            cposition = i + 1;
        else
            cposition = 1;
        end
    else
        cposition = 1; % in case things get missed up and neg index
    end
    ;   
  
function cposition = previous_row_end(positions,cposition)
    npos = length(positions(:,1));
    
    if cposition > 1
        i = cposition;
        while (i > 1) && (positions(i,1) == positions(cposition,1))
            i = i - 1;
        end
        if i > 1;
            cposition = i;
        else
            cposition = 1;
        end
    else
        cposition = 1;
    end
    ;
        
function nrows = number_of_rows(handles)
% Computes the number of rows in the current view
nrows = length(unique(handles.positions(handles.startpage:handles.endpage,1)));

function handles = row_backward(handles)
% Moves the row backwards in the classifiedaxes/browser view
    startpage = previous_row_start(handles.positions,handles.startpage);
    endpage = previous_row_end(handles.positions,handles.endpage);
    
    if  number_of_rows(handles) < handles.number_rows % indicates you are not at the last page
         handles.startpage = startpage;
         handles.endpage = length(handles.positions(:,1));
    elseif handles.startpage == 1;
        handles.startpage = 1;
        handles.endpage = handles.endpage;
    else
        handles.startpage = startpage;
        handles.endpage = endpage;
    end
%     guidata(gcbo,handles);
  ;
    
function handles = BrowseDirectory(handles)
% This function should be called only after precompute has been called
% it relies on their being a specicon file in the directory
%
% The function takes the current segments in the directory and loads their specicons
% into memory because the segments and icons are created in their order the
% order matches. This will need to be worked out better for the
% classification algorithms.

%load('-mat', 'specicons');
set_status(handles, ['Viewing all ', num2str(handles.nsegments),' segments']);
handles.mapindex = [1:handles.nsegments];
handles.image_list = handles.IconList;
handles = reposition_images(handles, handles.image_list);
handles.cpositions =  get_curr_position(handles);
handles = plot_classified_axes(handles, handles.image_list(handles.startpage:handles.endpage), handles.cpositions);
handles.mode='browse';

%guidata(gcbo,handles);


% --- Executes on button press in NextRowButton.
function NextRowButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextRowButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.cnrows > 1
    handles.startx = 1; % if you have been scrolling reset to defaults
    handles.endx = handles.classified_width;
    set(handles.SliderClassified, 'Value',0);
    handles = row_forward(handles);
    handles.cpositions =  get_curr_position(handles);
    handles = plot_classified_axes(handles, handles.image_list(handles.startpage:handles.endpage), handles.cpositions);
end
guidata(gcbo,handles);

% --- Executes on button press in PreviousRowButton.
function PreviousRowButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousRowButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.cnrows > 1
    handles.startx = 1;
    handles.endx = handles.classified_width;
    set(handles.SliderClassified, 'Value',0);
    handles = row_backward(handles);
    handles.cpositions =  get_curr_position(handles);
    handles = plot_classified_axes(handles, handles.image_list(handles.startpage:handles.endpage), handles.cpositions);
end
guidata(gcbo,handles);

% --- Executes on button press in LoadDirectoryButton.
function LoadDirectoryButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDirectoryButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%directoryname = uigetdir('./','Change directory');

directoryname = uigetdir;

if not(directoryname == 0)
    handles.directory = directoryname;
    set(handles.DirectoryEditBox, 'String',directoryname);
    cd(handles.directory); % this we will need to change
end
guidata(gcbo,handles);
;

function im_index = coordinate2index(handles,xpos,ypos)
% For the current images displayed tests if pointer position is in an image
% returns the index for that image based on a left to right ordering on
% that page
% First check if the coordinate is totally out of the range
if (xpos < 0) || (ypos < 0) || (ypos > handles.classified_height)  || (xpos > handles.classified_width)
    im_index = 0;
else
    im_index = 0;
    npos = length(handles.cpositions(:,1));
    i = 1;
    while (i <= npos) && (im_index < 1) 
        if (xpos >= handles.cpositions(i,2)) && (xpos <= (handles.cpositions(i,2) + handles.images_dim(handles.startpage + (i-1), 2)))
            if (ypos >= handles.cpositions(i,1)) && (ypos <= (handles.cpositions(i,1) + handles.images_dim(handles.startpage + (i-1),1)))
                im_index = i;
            end
        end
        i = i+1;
    end
end
;

function setnavigationbuttons(handles)
% Sets the navigation buttons based on where the pointer is

if not(handles.quickmode) 

if handles.NextIndex == handles.nsegments
    set(handles.PreviousSpectraButton,'Enable','off');
    set(handles.NextSpectraButton,'Enable','off');
elseif handles.NextIndex == handles.nsegments
    set(handles.NextSpectraButton,'Enable','off');
    set(handles.PreviousSpectraButton,'Enable','on');
elseif handles.NextIndex == 1
    set(handles.PreviousSpectraButton,'Enable','off');
    set(handles.NextSpectraButton,'Enable','on');
elseif (handles.NextIndex > 1) && (handles.NextIndex < handles.nsegments)
    set(handles.PreviousSpectraButton,'Enable','on');
    set(handles.NextSpectraButton,'Enable','on');

end

end

function setrowbuttons(handles)

if handles.startpage == 1
    set(handles.PreviousRowButton,'Enable','off');
elseif handles.startpage > 1
    set(handles.PreviousRowButton,'Enable','on');
end

if handles.endpage == length(handles.positions) % You are the last row
    set(handles.NextRowButton,'Enable','off');
elseif handles.cnrows <= handles.number_rows
    set(handles.NextRowButton,'Enable','off');
elseif handles.endpage < length(handles.positions)
    set(handles.NextRowButton,'Enable','on');
end
    
function DummyClassifyAxesClickCallBack(src,eventdata,handles)
% When the image is blank this allow you to select out of the class view
set(handles.ClassifyButton,'Enable','on');
set(handles.NewClassButton,'Visible','off');
set(handles.NextSpectraButton,'Enable','on');
set(handles.PreviousSpectraButton,'Enable','on');
set(handles.ModePopupMenu,'Enable','on');
set(handles.AutoClassifyButton,'Enable','on');
set(handles.QuickModeButton,'Enable','on');
set(handles.CompareToggleButton,'Enable','on');
handles.submode = 'select';
setnavigationbuttons(handles);
guidata(gcbo,handles);


function ClassifyAxesClickCallBack(src,eventdata)
% A Function for handling clicks to the axes
    handles = guidata(gcbo);
%     handles.mode
    %handles.submode
    %fprintf('\n');
    
    pos = get(handles.ClassifiedAxes,'CurrentPoint');
    cposition = coordinate2index(handles,pos(1,1),pos(1,2));
    
    if handles.quickmode
    if cposition == 0 % Selecting in the outside takes you out of classification mode
            set_status(handles,''); 
        else % You have selected an icon
            handles.startx = 1; % if you have been scrolling reset to defaults
            handles.endx = handles.classified_width;
            set(handles.SliderClassified,'Value',0);
            class = handles.classes(handles.mapindex(cposition + (handles.startpage - 1)));
            handles.segments(handles.NextIndex).class = class.name;
            handles.classes(handles.mapindex(cposition + (handles.startpage - 1))).nmembers = handles.classes(handles.mapindex(cposition + (handles.startpage - 1))).nmembers + 1;
            handles = jump_to_unclassified(handles);
            handles = configureclassview(handles,'select-class');
            if handles.lastsegment == handles.NextIndex % no more unclassified segments
                handles = quick_mode_exit(handles);
                set(handles.QuickModeButton,'Value',0);
                handles.quickmode = not(handles.quickmode);
            end
        end
    else % quick classify mode is off

    
    if strcmp(handles.mode,'browse') && (cposition > 0)
        handles.NextIndex = handles.mapindex((cposition - 1) + handles.startpage);
        handles=ConfigureClassSegment(handles);
        handles=ConfigureSpecPlot(handles);
    else
         %fprintf('%i\n', coordinate2index(handles,pos(1,1),pos(1,2)));
         %fprintf('%i, %i\n\n', pos(1,1),pos(1,2));
         ;
    end
    
    if strcmp('class-members',handles.mode) && strcmp('select',handles.submode) && (cposition > 0)
        handles.NextIndex = handles.mapindex((cposition - 1) + handles.startpage);
        handles=ConfigureClassSegment(handles);
        handles=ConfigureSpecPlot(handles);
    end
       
    %Show all members of a specific class
    if strcmp('class-view',handles.mode) && strcmp('select',handles.submode) && (cposition > 0)
        handles.startx = 1; % if you have been scrolling reset to defaults
        handles.endx = handles.classified_width;
        set(handles.SliderClassified,'Value',0);
%         set(handles.CompareToggleButton,'Visible','off');
        cindex = handles.mapindex(cposition + (handles.startpage - 1));
        class = handles.classes(cindex);
        handles.lastclass = cindex;
        handles = configureclassmembers(handles,class.name);
        handles.mode = 'class-members';
        handles.submode = 'select';
        handles = SetModePopupMenu(handles,'class members');
        set(handles.TypifyClassButton,'Visible','on');
        set(handles.RemoveClassButton,'Visible','off');
        set(handles.NextClassButton,'Visible','on');
        set(handles.RenameClassButton,'Visible','on');
        set_status(handles, ['Viewing ' num2str(length(handles.mapindex)),' members of ' class.name]);
    end
    
    if strcmp('class-view',handles.mode) && strcmp('select-class',handles.submode)
        if cposition == 0 % Selecting in the outside takes you out of classification mode
            set(handles.ClassifyButton,'Enable','on');
            set(handles.NewClassButton,'Visible','off');
            set_status(handles, ['Viewing all ', num2str(handles.nclasses),' classes']);
            
        else % You have selected an icon
            handles.startx = 1; % if you have been scrolling reset to defaults
            handles.endx = handles.classified_width;
            set(handles.SliderClassified,'Value',0);
            class = handles.classes(handles.mapindex(cposition + (handles.startpage - 1)));
            handles.segments(handles.NextIndex).class = class.name;
            handles.classes(handles.mapindex(cposition + (handles.startpage - 1))).nmembers = handles.classes(handles.mapindex(cposition + (handles.startpage - 1))).nmembers + 1;
            set(handles.ClassifyButton,'String','Declassify');
            handles = SetModePopupMenu(handles,'class members');
            handles = configureclassmembers(handles,class.name);
            handles.mode = 'class-members';
            handles.submode = 'select';
            set_status(handles,['Classified segment as ' class.name]);
        end
        set(handles.CompareToggleButton,'Enable','on');
        set(handles.QuickModeButton,'Enable','on');
        set(handles.AutoClassifyButton,'Enable','on');
        set(handles.ClassifyButton,'Enable','on');
        set(handles.NewClassButton,'Visible','off');
        set(handles.NextSpectraButton,'Enable','on');
        set(handles.PreviousSpectraButton,'Enable','on');
        set(handles.ModePopupMenu,'Enable','on');
        set(handles.SortText,'Visible','off');
        set(handles.SortPopupMenu,'Visible','off');
        handles.submode = 'select';
    end
    
    if strcmp(handles.mode,'class-members') && strcmp(handles.submode,'typify')
        if cposition == 0
            set(handles.TypifyClassButton, 'Enable','on');
        else
            handles.startx = 1; % if you have been scrolling reset to defaults
            handles.endx = handles.classified_width;
            set(handles.SliderClassified,'Value',0);
            sindex = handles.mapindex(cposition + (handles.startpage - 1));
            cindex = returnclassindex(handles,handles.segments(sindex).class);
            handles.classes(cindex).specfilename = handles.segments(sindex).specfilename;
            handles.classes(cindex).length = handles.segments(sindex).end - handles.segments(sindex).start;
            handles.classes(cindex).index = sindex;
            %load('-mat','specicons');
            handles.classes(cindex).iconS = handles.IconList{sindex};
            handles.subclass = 'xxx';
            set(handles.TypifyClassButton, 'Enable','on');
            set_status(handles,'');
        end
    end
    
    if strcmp('class-view', handles.mode) && strcmp('compare',handles.submode)
        if cposition == 0
            handles.mode = 'class-view';
            handles.submode = 'select';
            set(handles.CompareToggleButton,'Enable','on');
            set(handles.RemoveClassButton,'Enable','on');
            set(handles.NextSpectraButton,'Enable','on');
            set(handles.PreviousSpectraButton,'Enable','on');
            %set(handles.ModePopupMenu,'Enable','on');
        else
            handles.mode = 'comparison';
            if strcmp(handles.segments(handles.NextIndex).class,'')
                set(handles.ClassifyButton,'Enable','on');
            end
            handles.startx = 1; % if you have been scrolling reset to defaults
            handles.endx = handles.classified_width;
            set(handles.SliderClassified,'Value',0);
            cindex = handles.mapindex(cposition + (handles.startpage - 1));
            handles.image_list = {};
            handles.lastclass = cindex;
            handles.image_list{2} = handles.classes(cindex).iconS;
            %load('-mat','specicons');
            handles.image_list{1} = handles.IconList{handles.NextIndex};
            
            set_status(handles,['Comparing to ' handles.classes(cindex).name]);
            
            handles = reposition_images(handles, handles.image_list);
            handles.cpositions =  get_curr_position(handles);
            handles = plot_classified_axes(handles, handles.image_list(handles.startpage:handles.endpage), handles.cpositions);
            handles.submode = 'xxx';
        end
    end
 
    
    if strcmp(handles.mode,'class-view') && strcmp(handles.submode,'remove-class')
        if cposition == 0
            handles.submode = 'select';
            set(handles.RemoveClassButton, 'Enable','on');
            set_status(handles, ['Viewing all ', num2str(handles.nclasses),' classes']);
        else        
            cindex = handles.mapindex(cposition + (handles.startpage - 1));
            classname = handles.classes(cindex).name;
            
            answer = questdlg(['Remove class ' classname]);
            
            if strcmp(answer,'Yes')
                handles.startx = 1; % if you have been scrolling reset to defaults
                handles.endx = handles.classified_width;
                set(handles.SliderClassified,'Value',0);
                for i = 1:handles.nsegments
                    if strcmp(classname,handles.segments(i).class)
                        handles.segments(i).class = ''; 
                    end
                end
            
                if handles.nclasses > 1
                    handles.classes = [handles.classes(1:(cindex-1)) handles.classes((cindex+1):handles.nclasses)];
                    handles.nclasses = handles.nclasses - 1;
                    handles = configureclassview(handles,'select');
                    
                    set_status(handles, ['Viewing all ', num2str(handles.nclasses),' classes']);
                else
                    handles.classes = [];
                    handles.nclasses = 0;
                    handles.image_list = {};
                    set_status(handles,'');
                    set(handles.RemoveClassButton,'Visible','off');
                    handles.submode = 'xxx';
                    handles = blankaxes(handles);
                end
            
                if strcmp(handles.segments(handles.NextIndex).class,'')
                    set(handles.ClassifyButton,'String','Classify');
                end
                set(handles.RemoveClassButton, 'Enable','on');
            else
                handles.submode = 'select';
                set(handles.RemoveClassButton, 'Enable','on');
                set_status(handles, ['Viewing all ', num2str(handles.nclasses),' classes']);
            end
        end 
    end     
    
    end 
    
    if not(strcmp(handles.mode,'comparison'))
        setnavigationbuttons(handles);
    end
%     handles.mode
%     handles.submode
    guidata(gcf,handles);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%  Functions for setting up the classview  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = configureclassview(handles,submode)
% Configures the class-view for selecting classes
    handles.mode = 'class-view';
    handles.submode = submode;
%     if not(handles.quickmode)
%         set_status(handles, ['Select a class']);
%     end
    handles.mapindex = [1 : handles.nclasses];
    
    if strcmp(handles.sortclass, 'popularity')
        handles.mapindex = sortindexbypop(handles);
    elseif strcmp(handles.sortclass,'length');
        handles.mapindex = sortindexbylength(handles);
    else
        handles.mapindex = [1 : handles.nclasses];
    end
    handles = cindex2imagelist(handles);
    handles = reposition_images(handles, handles.image_list);
    
    if (strcmp(handles.sortclass,'length') && strcmp(handles.submode,'select-class')) || (strcmp(handles.sortclass,'length') && strcmp(handles.submode,'compare'))
        % jump to the segment with the closest size match
        i = 1;
        while (i <= handles.nclasses) && (handles.classes(handles.mapindex(i)).length >= (handles.segments(handles.NextIndex).end - handles.segments(handles.NextIndex).start))
            i = i + 1;
        end
        cnrow = which_row(handles.positions,i-1); % get the current row
        
        if cnrow > handles.number_rows % The closes size segment is not in view
            for i = 1:(cnrow - handles.number_rows) % matching size row is last
                handles = row_forward(handles);
            end
            
            if not(handles.endpage == length(handles.positions)) % if not at the last row position so that larger and smaller rows match 
                nrows = floor(handles.number_rows / 2);
                for i = 1:nrows
                    handles = row_forward(handles);
                end
            end
            
            
        end
    end
    
    handles.cpositions = get_curr_position(handles);
    handles = plot_classified_axes(handles, handles.image_list(handles.startpage:handles.endpage), handles.cpositions);
    %guidata(gcbo,handles);
;

function handles = SetModePopupMenu(handles,viewstring)
popmodes = get(handles.ModePopupMenu,'String');

find_index = 0;
i = 1;
while (i <= length(popmodes)) && not(strcmp(popmodes(i),viewstring))
    i = i + 1;
end

if i <= length(popmodes) % Don't do anything if the string cannot be found
    set(handles.ModePopupMenu,'Value',[i]);
end

% --- Executes on selection change in ModePopupMenu.
function ModePopupMenu_Callback(hObject, eventdata, handles)
% Configure call back
% hObject    handle to ModePopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ModePopupMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ModePopupMenu

handles.startx = 1; % if you have been scrolling reset to defaults
handles.endx = handles.classified_width;
set(handles.SliderClassified,'Value',0);
%set(handles.CompareToggleButton,'Visible','off');
set(handles.RenameClassButton, 'Visible', 'off');

nmode = get(hObject,'Value');
modeview = get(hObject,'String');


if strcmp(modeview(nmode),'all')
    handles.mode = 'browse';
    handles.submode = 'select';
    handles=BrowseDirectory(handles);
    set(handles.RemoveClassButton,'Visible','off');
    set(handles.TypifyClassButton, 'Visible', 'off');
    set(handles.NextClassButton,'Visible','off');

elseif strcmp(modeview(nmode),'class view')
    set(handles.NextClassButton,'Visible','off');
    if handles.nclasses >= 1
        handles = configureclassview(handles,'select');
        handles.mode = 'class-view';
        handles.submode = 'select';
        set(handles.RemoveClassButton,'Visible','on');
        set(handles.RemoveClassButton,'Enable','on');
        set(handles.CompareToggleButton,'Visible','on');
        set(handles.CompareToggleButton,'Enable','on');
        set(handles.TypifyClassButton, 'Visible', 'off');
        set_status(handles,['Viewing all ' num2str(handles.nclasses) ' classes'])
    else % empty axes
        handles = blankaxes(handles);
        handles.mode = 'class-view';
        handles.submode = 'xxx';
        set_status(handles,'');
    end
elseif strcmp(modeview(nmode),'unclassified') || (strcmp(handles.segments(handles.NextIndex).class,'') && strcmp(modeview(nmode),'class members'))
    handles = configureclassmembers(handles,'');
    set_status(handles,['A total of ' num2str(length(handles.mapindex)) ' unclassified segments ']);
    handles = SetModePopupMenu(handles,'unclassified');
    handles.mode = 'browse';
    handles.submode = 'select';
    set(handles.RemoveClassButton, 'Visible', 'off');
    set(handles.TypifyClassButton,'Visible','off');
    set(handles.NextClassButton,'Visible','off');
elseif strcmp(modeview(nmode),'class members')
    handles = configureclassmembers(handles,handles.segments(handles.NextIndex).class);
    set_status(handles,['Viewing ' num2str(length(handles.mapindex)) ' members of ' handles.segments(handles.NextIndex).class]);
    handles.mode = 'class-members';
    handles.submode = 'select';
    handles.lastclass = get_class_index(handles,handles.segments(handles.NextIndex).class);
    set(handles.TypifyClassButton,'Visible', 'on');
    set(handles.RemoveClassButton, 'Visible', 'off');
    set(handles.RenameClassButton,'Visible','on');
    set(handles.NextClassButton,'Visible','on');
end
guidata(gcbo,handles);

% --- Executes during object creation, after setting all properties.
function ModePopupMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModePopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function set_status(handles, statusstring)
    set(handles.SegInfoText,'String',statusstring);
    set(handles.SegInfoText,'Visible','on');
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for filtering by class type    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = configureclassmembers(handles,classname)

handles.mapindex = select_class(handles,handles.segments,classname);

if length(handles.mapindex) > 0
    handles.image_list = {};
    %load('-mat','specicons');
        
    for i = 1:length(handles.mapindex)
        handles.image_list(i) = handles.IconList(handles.mapindex(i));
    end
    handles = reposition_images(handles, handles.image_list);
    handles.cpositions =  get_curr_position(handles);
    %
    %handles.mode = 'class-members';
    handles = plot_classified_axes(handles, handles.image_list(handles.startpage:handles.endpage), handles.cpositions);
%handles = SetModePopupMenu(handles,'class members');
else
    handles = blankaxes(handles);
    %handles.mode = 'class-view';
    %handles.submode = 'select';
    set_status(handles,''); 
end
;

function indexfilter = select_class(handles,segments,classname)
% Returns an index array of original addresses of segments which are members of a
% specified class
indexfilter = [];
nclassmembers = 0;
for i = 1:handles.nsegments
    if strcmp(segments(i).class,classname)
        nclassmembers = nclassmembers + 1;
        indexfilter(nclassmembers) = i;
    end
end

function cindex = get_class_index(handles,classname);
cindex = 0;
i = 1;
while (i <= handles.nclasses) && not(strcmp(handles.classes(i).name,classname))
    i = i + 1;
end
cindex = i;

function value = mapind(index)
% Map the index value back to its original value
value = handles.mapindex(index);

%%&& not(strcmp(segments(i).class,''));
function write_syllable_database(handles)
filename = [handles.baseclassname '.dat'];
fid = fopen(filename,'wt');

classes = handles.classes;

if fid > -1
    for i = 1:handles.nsegments
        segment = handles.segments(i);
        wavfile = segment.wavfile;
        specfilename = segment.specfilename;
        seg = [ num2str(segment.start) '\t' num2str(segment.end)];
        classname = [ segment.class ];
        typify = '';
        nclasses = length(classes);
        j = 1;
        while (nclasses > 0) && (j <= nclasses) && not(strcmp(segment.specfilename,classes(j).specfilename))
            j = j + 1;
        end
        if j <= length(classes) % found a match
            classes = [classes(1:j-1) classes(j+1:nclasses)]; % shorten the classes
            typify = '*'; % indicates that this is typological class
        end
        fprintf(fid,[specfilename '\t' wavfile '\t' seg '\t' classname '\t' typify '\n']);
    end
    fclose(fid);
end
;

function data = read_syllable_database(handles)
filename = [handles.baseclassname '.dat'];
fid = fopen(filename,'rt');
data = textscan(fid,'%s %s %n %n %s %s', 'delimiter','\t');
data = [data(1), data(5), data(6)]; % throw out extra stuff which will be useful for external analysis
fclose(fid);
;

function handles = merge_syllable_database(handles,data)
% Merge the syllable list with the loaded database

%load('-mat','specicons');
specfilenames = data{1};
classnames = data{2};
typifies = data{3};
ndata = length(specfilenames);
%classnum = 1;
handles.classes = [];
handles.nclasses = 0;
for i = 1:ndata
    j=1; % allow for no matches and allow for 
    while (j < handles.nsegments) && not(strcmp(specfilenames(i),handles.segments(j).specfilename))
        j = j + 1;
    end
    if j <= handles.nsegments
        handles.segments(j).class = classnames{i};
        if  strcmp(typifies(i),'*') % this is the type class
            %classnum = classnum + 1;
            class.specfilename = specfilenames{i};
            class.name = classnames{i};
            class.index = j;
            class.length = handles.segments(j).end - handles.segments(j).start;
            class.iconS = handles.IconList{j};
            class.nmembers = 0; % will update shortly
            handles.nclasses = handles.nclasses + 1;
            handles.classes = [handles.classes class];
        end    
    end
end

% Now that we have the classes defined update the number of members
for i = 1:handles.nclasses
    for j = 1 : handles.nsegments
        if strcmp(handles.classes(i).name,handles.segments(j).class)
            handles.classes(i).nmembers = handles.classes(i).nmembers + 1;
        end
    end
end


% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

write_syllable_database(handles); % save the database
save_configuration(handles,handles.configfile); % save the current configuration


% --- Executes on slider movement.
function SliderClassified_Callback(hObject, eventdata, handles)
% hObject    handle to SliderClassified (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

xposition = round(get(hObject,'Value'));
%xposition
handles.startx = 1 + xposition;
handles.endx  = xposition + handles.classified_width;
handles.cpositions =  get_curr_position(handles);
classview = handles.classmatrix(:,handles.startx:handles.endx); % will cut overhang
axes(handles.ClassifiedAxes);
handles.hiclass = image(classview);

set(handles.hiclass,'ButtonDownFcn',{@ClassifyAxesClickCallBack});
set(handles.ClassifiedAxes,'XTick',[]);
set(handles.ClassifiedAxes,'YTick',[]);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SliderClassified_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderClassified (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in SortPopupMenu.
function SortPopupMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SortPopupMenu (see GCBO) eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns SortPopupMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SortPopupMenu

sortlist = get(hObject,'String');
sortn = get(hObject,'Value');
sortby = sortlist(sortn);

if strcmp(sortby,'original')
    handles.sortclass = 'original';
elseif strcmp(sortby,'by length')
    handles.sortclass = 'length';
elseif strcmp(sortby,'by popularity')
    handles.sortclass = 'popularity';
end

handles = configureclassview(handles,handles.submode);

guidata(gcbo,handles);

% --- Executes during object creation, after setting all properties.
function SortPopupMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SortPopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function segments = rename_segments(segments,oldname,newname)

nsegments = length(segments);
for i = 1 : nsegments
    if strcmp(segments(i).class,oldname)
        segments(i).class = newname;
    end
end

function cln = doesclassexist(classname,classes)
% tests whether the current class name already exists if does not exist
% returns 1 if it does exist returns 0

i = 1;
nclasses = length(classes);

while (i <= nclasses) && (strcmp(classname,classes(i).name))
    i = i + 1;
end

if i <= nclasses
    cln = i;
else
    cln = 0;
end

% --- Executes on button press in RenameClassButton.
function RenameClassButton_Callback(hObject, eventdata, handles)
% hObject    handle to RenameClassButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

class = handles.classes(handles.lastclass);
answer = inputdlg({'Class name'},'Edit class name',1,{class.name});

% check the answer
sizeanswer = size(answer);

secondanswer = '';

if sizeanswer(1)
    if not(strcmp(class.name,answer{1}))
        classnexists = doesclassexist(answer{1},handles.classes);
        if classnexists
            secondanswer = questdlg('Class name already exists. Do you want to merge the two classes.')
        end
        if strcmp(secondanswer,'Yes') || strcmp(secondanswer,'')
            segments = rename_segments(handles.segments,class.name,answer{1});
            class.name = answer{1};
            cindex = handles.lastclass;
            handles.classes(cindex) = class;
            handles.segments = segments;
            
            if strcmp(secondanswer,'Yes') % functionality for merging two classes
                handles.classes(classnexists).nmembers = handles.classes(classnexists).nmembers + handles.classes(cindex).nmembers;
                handles.classes = [handles.classes(1:(cindex-1)) handles.classes((cindex+1):handles.nclasses)];
                handles.nclasses = handles.nclasses - 1;
                class.nmembers = handles.classes(classnexists).nmembers; % this is parasitic code
            end
        end
        handles = configureclassmembers(handles,class.name);
        set_status(handles, ['Viewing ' num2str(class.nmembers),' members of ' class.name]);
    end 
end
guidata(gcbo,handles);


function spectras = subsamplespectra(spectra,lowerfreq,upperfreq,freqrange)
% Returns a subsampled frequency of the spectra

freqsamples = length(spectra(:,1));
freqratio = freqsamples / (freqrange(2) - freqrange(1));

lowerfreqsamp = round(lowerfreq * freqratio) + 1;
upperfreqsamp = round(upperfreq * freqratio) + 1;

if lowerfreqsamp < 1 % make sure we are not out of range
    lowerfreqsamp = 1;
end

if upperfreq > freqsamples
    upperfeqsamp = freqsamples;
end

spectras = spectra(lowerfreqsamp:upperfreqsamp,:);

function handles = generate_subsamples_icons(handles)
segments = handles.segments;
IconListf = {};
hw = waitbar(0,'Zooming spectra. . .');
for i = 1:handles.nsegments
    load('-mat',segments(i).specfilename);
    Ssub = subsamplespectra(S,handles.lowerfreq,handles.upperfreq,handles.fpass);
    IconListf{i} = iconify_spec(Ssub,handles.ispecheight);
    waitbar(i/handles.nsegments);
end
close(hw);
handles.IconListf = IconListf;


function handles = ZoomSpectra(handles,status)
    if strcmp(status,'Zoom in')
        if (length(handles.IconListf) == 0) || (handles.rezoom) 
            handles = generate_subsamples_icons(handles);
            handles.rezoom = logical(0);
        end
        handles.FullIconList = handles.IconList;
        set(handles.ZoomButton,'String','Zoom out');
        handles.IconList = handles.IconListf;
    elseif strcmp(status,'Zoom out');
        
        if handles.rezoom
            handles = generate_subsamples_icons(handles);
            handles.rezoom = logical(0);
            set(handles.ZoomButton,'String','Zoom out');
            handles.IconList = handles.IconListf;
        else
            handles.IconList = handles.FullIconList;
            set(handles.ZoomButton,'String','Zoom in');
        end
    end

    for j = 1:handles.nclasses % switch over class icons
        handles.classes(j).iconS = handles.IconList{handles.classes(j).index};  
    end

    nimages = length(handles.mapindex);

    % this is kind of ugly
    if not(strcmp(handles.mode,'class-view')) && not(strcmp(handles.mode,'comparison')) % update images
        for i = 1:nimages
            handles.image_list{i} = handles.IconList{handles.mapindex(i)};
        end
    elseif strcmp(handles.mode,'comparison')
        handles.image_list{1} = handles.IconList{handles.NextIndex};
        handles.image_list{2} = handles.IconList{handles.classes(handles.lastclass).index};
    else
        for i = 1:nimages
            handles.image_list{i} = handles.IconList{handles.classes(handles.mapindex(i)).index};
        end
    end

    handles=get_and_plot(handles,handles.segments(handles.NextIndex));
    handles = plot_classified_axes(handles, handles.image_list(handles.startpage:handles.endpage), handles.cpositions);



% --- Executes on button press in ZoomButton.
function ZoomButton_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

status = get(handles.ZoomButton,'String');

handles = ZoomSpectra(handles,status);

guidata(gcbo,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for computing class statistics      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [meanl sdl] = stat_lengths(handles)
lengthsarray = [];
for i = 1:length(handles.mapindex)
    lengthsarray(i) = handles.segments(i).end - handles.segments(i).start;
end
meanl = mean(lengthsarray);
sdl = sd(lengthsarray);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for quick classify mode                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = quick_mode_exit(handles)
    set(handles.SkipButton,'Visible','off');
    set(handles.UndoButton,'Visible','off');
    set(handles.NewQuickButton,'Visible','off');
    
    set(handles.ViewText,'Visible','on');
    set(handles.ModePopupMenu,'Visible','on');
    set(handles.SegInfoText,'Visible','on');
    
    set(handles.ClassifyButton,'Enable','on');
    set(handles.RemoveClassButton,'Enable','on');
    set(handles.RemoveClassButton,'Visible','on');
    set(handles.CompareToggleButton,'Enable','on');
    set(handles.AutoClassifyButton,'Enable','on');
    set(handles.SortPopupMenu,'Visible','off');
    set(handles.ModePopupMenu,'Value',3);
    
    handles.submode='select';
    
    set_status(handles,['Viewing all ' num2str(handles.nclasses) ' classes']);
    
    setnavigationbuttons(handles);


% --- Executes on button press in QuickModeButton.
function QuickModeButton_Callback(hObject, eventdata, handles)
% hObject    handle to QuickModeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if not(handles.quickmode) % turn on quick mode
    
    handles.quickmode = 1;
    
    % Turn off top header
    set(handles.ViewText,'Visible','off');
    set(handles.ModePopupMenu,'Visible','off');
    set(handles.SegInfoText,'Visible','off');
    
    % Make visible quick mode buttons
    set(handles.SkipButton,'Visible','on');
    set(handles.SkipButton,'Enable','on');
    set(handles.UndoButton,'Visible','on');
    set(handles.UndoButton,'Enable','on');
    set(handles.NewQuickButton,'Visible','on');
    set(handles.NewQuickButton,'Visible','on');
    
    set(handles.SortPopupMenu,'Visible','on');
    
    % Disable regular mode functions
    set(handles.RemoveClassButton,'Visible','off');
    set(handles.NextClassButton,'Visible','off');
    set(handles.NextSpectraButton,'Enable','off');
    set(handles.PreviousSpectraButton,'Enable','off');
    set(handles.ClassifyButton,'Enable','off');
    set(handles.TypifyClassButton,'Visible','off');
    set(handles.RenameClassButton,'Visible','off');
    set(handles.CompareToggleButton,'Enable','off');
    set(handles.AutoClassifyButton,'Enable','off');
    
    % Will need to disable buttons underneath classified axes
    if not(strcmp(handles.segments(handles.NextIndex).class,'')) 
        handles = jump_to_unclassified(handles);
        handles.lastsegment = handles.NextIndex;
    end
    
    if handles.nclasses >= 1
        handles = configureclassview(handles,'select-class');
        handles.mode = 'class-view';
        handles.submode = 'select-class';
    else % empty axes
        handles = blankaxes(handles);
        handles.mode = 'class-view';
        handles.submode = 'xxx';
    end
    
    
else
    handles.quickmode = 0;
    handles = quick_mode_exit(handles);
     
    % Will need to intellegently enable buttons underneath classified axes
end

guidata(gcbo,handles);

function handles=jump_to_unclassified(handles)
% Jumps to the next unclassified segment
currindex = handles.NextIndex;
i = currindex;
while (mod(i,handles.nsegments)+1 ~= currindex) && not(strcmp(handles.segments(mod(i,handles.nsegments) + 1).class,''))
    i = i + 1;
end

if not(mod(i,handles.nsegments)+1 == currindex) % there are some unclassified segments
    handles.NextIndex = mod(i,handles.nsegments) + 1;
    handles=ConfigureClassSegment(handles);
    handles=ConfigureSpecPlot(handles);
    set(handles.ClassifyButton,'Enable','off');
end
%get(handles.QuickModeButton,'Value')
handles.lastsegment = currindex;

% --- Executes on button press in SkipButton.
function SkipButton_Callback(hObject, eventdata, handles)
% hObject    handle to SkipButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = jump_to_unclassified(handles);

if handles.lastsegment == handles.NextIndex % no more unclassified segments
    handles = quick_mode_exit(handles);
    set(handles.QuickModeButton,'Value',0);
    handles.quickmode = not(handles.quickmode);
end

handles = configureclassview(handles,'select-class');

guidata(gcbo,handles);
% --- Executes on button press in UndoButton.
function UndoButton_Callback(hObject, eventdata, handles)
% hObject    handle to UndoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.NextIndex = handles.lastsegment;
handles=ConfigureClassSegment(handles);
handles=ConfigureSpecPlot(handles);

handles = quick_mode_exit(handles);

set(handles.QuickModeButton,'Value',0);
handles.quickmode = not(handles.quickmode);

guidata(gcbo,handles);

% --- Executes on button press in NewQuickButton.
function NewQuickButton_Callback(hObject, eventdata, handles)
% hObject    handle to NewQuickButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = add_new_class(handles,handles.segments(handles.NextIndex));
handles.segments(handles.NextIndex).class = handles.classes(handles.nclasses).name;

set(handles.ClassifyButton,'Enable','off');

handles = jump_to_unclassified(handles);
handles = configureclassview(handles,'select-class');

if handles.lastsegment == handles.NextIndex % no more unclassified segments
    handles = quick_mode_exit(handles)
    set(handles.QuickModeButton,'Value',0);
    handles.quickmode = not(handles.quickmode);
end

guidata(gcbo,handles);



% --- Executes on button press in CompareToggleButton.
function CompareToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to CompareToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CompareToggleButton

state = get(hObject,'Value');

if handles.nclasses > 0

 if state % enter compare mode
    handles.mode = 'class-view';
    set(handles.RemoveClassButton,'Enable','off');
    set(handles.RenameClassButton,'Visible','off');
    set(handles.NextClassButton,'Visible','off');
    set(handles.ModePopupMenu,'Enable','off');
    set(handles.ClassifyButton,'Enable','off');
    set(handles.QuickModeButton,'Enable','off');
    set_status(handles,'Select a class to compare segment against');
    set(handles.NextSpectraButton,'Enable','off');
    set(handles.PreviousSpectraButton,'Enable','off');
    set(handles.TypifyClassButton, 'Visible', 'off');
    set(handles.AutoClassifyButton,'Enable','off');
    
    handles = configureclassview(handles,'compare');
 else
    if strcmp(handles.mode,'comparison') % go back to the class view
        set(hObject,'Value',1);
        handles.startx = 1; % if you have been scrolling reset to defaults
        handles.endx = handles.classified_width;
        
        set(handles.SliderClassified, 'Value',0);
        handles = configureclassview(handles,'compare');
        handles.mode = 'class-view';
        set_status(handles,'Select a class to compare segment against');
    elseif strcmp(handles.mode,'class-view') % exit the compare mode entirely
        set_status(handles, ['Viewing all ', num2str(handles.nclasses),' classes']);
        set(handles.ModePopupMenu,'Enable','on');
        set(handles.QuickModeButton,'Enable','on');
        setnavigationbuttons(handles);
        SetModePopupMenu(handles,'class view');
        set(handles.RemoveClassButton,'Enable','on');
        set(handles.RemoveClassButton,'Visible','on');
        set(handles.ClassifyButton,'Enable','on');
        set(handles.AutoClassifyButton,'Enable','on');
        handles.submode = 'select';
    end
 end

else
    set(hObject,'Value',0);
end
guidata(hObject, handles); 

function handles = recompute_classifiedaxes(handles)
% Recomputes the width and height based on a new size for the classified
% axes

    classaxpos=get(handles.ClassifiedAxes,'Position');
    
    set(handles.SliderClassified,'Value',0);
    
    handles.startx = 1;
    
    handles.classified_width = round(handles.classified_width_density * classaxpos(3));
    handles.classified_height = round(handles.classified_height_density * classaxpos(4));
    
    oldpos = handles.startpage;
    handles = reposition_images(handles, handles.image_list); % try to keep rows matched
    
    cnrow = which_row(handles.positions,oldpos);
    for i = 1:(cnrow-1)
           handles = row_forward(handles);
    end
    
    handles.cpositions =  get_curr_position(handles);
    handles = plot_classified_axes(handles, handles.image_list(handles.startpage:handles.endpage), handles.cpositions);
    
function reposy(guielement,deltay)
    oldpos = get(guielement,'Position');
    set(guielement,'Position',[oldpos(1), oldpos(2) + deltay, oldpos(3), oldpos(4)]);

function reposelementsy(handles,deltay)
    reposy(handles.PrecomputeButton,deltay);
    reposy(handles.DirectoryEditBox,deltay);
    reposy(handles.LoadDirectoryButton,deltay);
    reposy(handles.SaveButton,deltay);
    reposy(handles.ConfigureButton,deltay);
    reposy(handles.ViewText,deltay);
    reposy(handles.ModePopupMenu,deltay);
    reposy(handles.SegInfoText,deltay);
    reposy(handles.ToClassifyPanel,deltay);
    reposy(handles.NewQuickButton,deltay);
    reposy(handles.UndoButton,deltay);
    reposy(handles.SkipButton,deltay);
    reposy(handles.SortText,deltay);
    reposy(handles.SortPopupMenu,deltay);
    reposy(handles.NextClassButton,deltay);
    
function ResizeFcn(h, eventdata, handles, varargin)

    handles = guidata(gcbo);
    
    originalsize = handles.originalsize;
    newsize = get(h,'Position');
    classaxpos = get(handles.ClassifiedAxes,'Position');
     
    if handles.originalsize(3) > newsize(3) % if form has smaller width bounce back
        newsize(3) = originalsize(3);
        classaxpos(3) = handles.originalaxessize(3);
        
        sliderpos = get(handles.SliderClassified,'Position'); % Update the slider
        sliderpos(3) = classaxpos(3);
        set(handles.SliderClassified,'Position',sliderpos);
        
    else % form is larger
        deltax = newsize(3) - handles.prevsize(3);
        classaxpos(3) = classaxpos(3) + deltax; % update width of the axes
        
        sliderpos = get(handles.SliderClassified,'Position'); % Update the slider
        sliderpos(3) = sliderpos(3) + deltax;
        set(handles.SliderClassified,'Position',sliderpos);
        
    end

    if originalsize(4) > newsize(4) % if form has a smaller height bounce back
        newsize(2) = newsize(2) + newsize(4) - originalsize(4);
        
        deltay =  originalsize(4) - handles.prevsize(4);
        newsize(4) = originalsize(4);
        classaxpos(4) = classaxpos(4) + newsize(4) - handles.prevsize(4);
        reposelementsy(handles,deltay);
        
    else % form is larger we need to reposition elements
        deltay = newsize(4) - handles.prevsize(4);
        reposelementsy(handles,deltay);
        classaxpos(4) = classaxpos(4) + newsize(4) - handles.prevsize(4);
    end
    
    set(h,'Position',newsize);
    set(handles.ClassifiedAxes,'Position',classaxpos);
    
    if handles.fixed || handles.blank % allow fixed resizing
        % do nothing except update the density measurements
        handles.classified_width_density = handles.classified_width / classaxpos(3);
        handles.classified_height_density = handles.classified_height / classaxpos(4);
    else % reposition images
        handles = recompute_classifiedaxes(handles);
    end
        
    handles.prevsize = newsize;

    guidata(gcbo,handles);

function truth = truthrange(range)
    truth = range(1) && range (2);

% --- Executes on button press in ConfigureButton.
function ConfigureButton_Callback(hObject, eventdata, handles)
% hObject    handle to ConfigureButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlenb  vcccs    structure with handles and user data (see GUIDATA)


    configs = configure_classify(handles.lowerfreq,handles.upperfreq,handles.classified_height,handles.classified_width,handles.Fs,handles.movingwin,handles.tapers,handles.fpass,handles.fixed);
    if not(isempty(configs) || length(configs) == 1)  % if okay is pressed
        lowerfreq= configs{1};
        upperfreq= configs{2};
        classified_height= configs{3};
        classified_width= configs{4};
        Fs=configs{5};
        movingwin= configs{6};
        tapers= configs{7};
        fpass = configs{8};
        
        rezoom = 0;
    
        if not(handles.lowerfreq == lowerfreq)
            handles.lowerfreq = lowerfreq;
            rezoom = 1;
        end
    
        if not(handles.upperfreq==upperfreq)
            handles.upperfreq=upperfreq;
            rezoom=1;
        end
    
        redraw = 0;
        if not(handles.classified_height==classified_height)
            handles.classified_height=classified_height;
            redraw = 1;
        end
    
        if not(handles.classified_width==classified_width)
            handles.classified_width=classified_height;
            redraw = 1;
        end
    
        recompute = 0; % check to see if the spectra needs to be recomputed
    
        if not(handles.Fs==Fs)  
            handles.Fs = Fs;
            recompute = 1;
        end
    
        if not(truthrange(handles.movingwin == movingwin))
            handles.movingwin = movingwin;
            recompute = 1;
        end
    
        if not(truthrange(handles.tapers==tapers))
            handles.tapers=tapers;
            recompute = 1;
        end
    
        if not(truthrange(handles.fpass==fpass))
            handles.fpass=fpass;
            recompute = 1;
        end
    
        
        if recompute % The spectra need to be recomputed
            status = questdlg('Recompute spectra with changed parameters');
            if not(isempty(status))
                if strcmp(status,'Yes')
                     handles = precompute_AllSpectra(handles); 
                     save_configuration(handles,handles.configfile); % Save configuration no matter what
                    
                     handles=ConfigureSpecPlot(handles); % automatically force into browse mode
                     set(handles.RemoveClassButton,'Visible','off');
                     set(handles.TypifyClassButton, 'Visible', 'off');
                     set(handles.NextClassButton,'Visible','off');
                     
                     handles=BrowseDirectory(handles);
                     set(handles.ModePopupMenu,'Value',1);
                     
                     if handles.nclasses >= 1 % sync changes to the classview
                        for j = 1:handles.nclasses
                            handles.classes(j).iconS = handles.IconList{handles.classes(j).index};
                        end
                     end
                     
                     set(handles.RemoveClassButton,'Visible','off'); % update buttons
                     set(handles.TypifyClassButton, 'Visible', 'off');
                     set(handles.NextClassButton,'Visible','off');                     
                     
                     redraw = logical(1);
                end
            end
        end
        
        if rezoom
            if strcmp(get(handles.ZoomButton,'String'),'Zoom out')
                handles.rezoom = logical(1);
                handles = ZoomSpectra(handles,'Zoom out');
            else
                handles.rezoom = logical(1); % Otherwise wait to user hits the rezoom button
            end
        end
        
        if not(handles.configschanged) % handles the case were the configurations were changed already
            if not(recompute)
                handles.configschanged = logical(1);
            end
        end
        
        if redraw % Redraw the axes
            classaxpos = get(handles.ClassifiedAxes,'Position');
            handles.classified_width_density = handles.classified_width / classaxpos(3);
            handles.classified_height_density = handles.classified_height / classaxpos(4);
            handles = recompute_classifiedaxes(handles);
         end
        
        handles.fixed = configs{9};
        
end
guidata(gcbo,handles);

function status = save_configuration(handles,filename)
    
    Fs = handles.Fs;
    movingwin = handles.movingwin;
    tapers = handles.tapers;
    fpass = handles.fpass;
    
    try
        save(filename,'Fs','movingwin','tapers','fpass','-mat');
    catch
        status = 1;
    end
    
function handles = load_configuration(handles,filename)
    try
        load('-mat',filename)
        if handles.Fs == Fs || handles.movingwin == movingwin || handles.tapers == tapers || handles.fpass == fpass
            handles.configshavechanged = logical(1);
        end 
        
        handles.Fs = Fs;
        handles.movingwin = movingwin;
        handles.tapers = tapers;
        handles.fpass = fpass;
    catch
        handles = handles;
    end


% --- Executes on button press in NextClassButton.
function NextClassButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextClassButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.lastclass = handles.lastclass + 1;

if handles.lastclass > handles.nclasses
    handles.lastclass = 1;
end

class = handles.classes(handles.lastclass);
handles = configureclassmembers(handles,class.name);
set_status(handles,['Viewing ' num2str(length(handles.mapindex)) ' members of ' class.name]);
handles.mode = 'class-members';
handles.submode = 'select';

guidata(gcbo,handles);

function handles = features_segments(handles)

h=waitbar(0,'Computing Data Features. . . ');
for i = 1:handles.nsegments
    startf = handles.segments(i).start;
    endf = handles.segments(i).end;
    cepestral = cepsfromspectra(handles.segments(i).specfilename,handles.ncepestral);
    
    % Below shows how additional data features can be included in the
    % exporting and classifying using the software
    
    % The commented code shows how the wave file can be read in to
    % calculate additional data features
    
    %     segment = handles.segments(i);
    %     data = wavread(segment.wavfile,round(handles.Fs * [segment.start segment.end]));
        
    waitbar(i/handles.nsegments);
    
    handles.segments(i).features = cepestral';
    
    % Three additional data features can be added for example as following
    %
    % additional_features = compute_features(data,3);
    % handles.nfeatures = handles.nfeatures + 3;
    % handles.segments(i).features = [handles.segments(i).features additional_features'];
end
close(h);

function coefs = cepsfromspectra(specfilename,ncepestral)
  load('-mat',specfilename);
  sbase = (mean(min(exp(Spre))) + mean(min(exp(Spost)))) / 2; % might not be a good
                                           % baseline because of
                                           % the intensity of the
                                           % sound around the syllable that
                                           % is why I took the average
                                           % minimum value
  
  %Sdiff = S - sbase;           
  
  Sdiff = exp(S) - sbase; % arithmetic average
  spectra = log(mean(Sdiff')); % average across time to generate spectra
    
  %spectra=mean(Sdiff');
  
  cepstrum = real(ifft(spectra));
  
  % alternative methods for computing the cepstrum
  %real(fft(hamming(length(spectra) .* spectra )));
  %cepstrum = real(fft(spectra))
  
  coefs = cepstrum(1:ncepestral); % this I understand is the
                                   % correct way 
  
function classifymatrix = generate_classify_mat(handles)
    classifymatrix = zeros(handles.nsegments,handles.nfeatures+1);
    for i = 1:handles.nsegments
        seglength = handles.segments(i).end - handles.segments(i).start;
        classifymatrix(i,1:1+handles.nfeatures) = [seglength handles.segments(i).features(1:handles.nfeatures)'];
    end
;
    
function handles = generate_classes_auto(handles,classification)
    nclasses = length(unique(classification));
    classesfound = zeros(nclasses,1); % used to check which classes exist
    handles.classes = [];
    handles.nclasses = 0;
    j = 1;
    for i = 1:handles.nsegments
        if classification(i) == 0 % segment has not been classified
            handles.segments(i).class = '';
        else 
            if isempty(classesfound(find(classesfound == classification(i)))) % a new class is found
                class.specfilename = handles.segments(i).specfilename;
                class.name = newclassname(handles);
                class.index = i;
                class.length = handles.segments(i).end - handles.segments(i).start;
                class.iconS = handles.IconList{i};
                class.nmembers = 1;
                handles.segments(i).class = class.name;
                handles.classes = [handles.classes class];
                classesfound(j) = classification(i);
                j = j + 1;
                handles.nclasses = handles.nclasses + 1;
            else % segment belongs to a class that already exists
                classindex = find(classification(i) == classesfound);
                handles.segments(i).class = handles.classes(classindex).name;
                handles.classes(classindex).nmembers = handles.classes(classindex).nmembers + 1; 
            end
        end
    end
   ;
        
% --- Executes on button press in AutoClassifyButton.
function AutoClassifyButton_Callback(hObject, eventdata, handles)
% hObject    handle to AutoClassifyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if length(handles.segments(1).features) == 0;  % compute cepestral coefficients if not computed before
    handles = features_segments(handles);
end

matrix2classify = generate_classify_mat(handles);

classification = auto_classify(matrix2classify,handles.nfeatures);

if length(classification) > 1  %&& not(classification == 0)
    if not(isempty(handles.classes))
        answer = questdlg('A classification already exists. Do you want to replace the current classification?');
        if strcmp(answer,'Yes')
            handles = generate_classes_auto(handles,classification);
        end
    else
        handles = generate_classes_auto(handles,classification);
        
    end
    
    handles=configureclassview(handles,'select');
    set(handles.ClassifyButton,'String','Declassify'); % all classes are now classified
    set(handles.RemoveClassButton,'Visible','on');
    set(handles.RemoveClassButton,'Enable','on');
    set(handles.CompareToggleButton,'Visible','on');
    set(handles.CompareToggleButton,'Enable','on');
    set(handles.TypifyClassButton, 'Visible', 'off');
    set_status(handles,['Viewing all ' num2str(handles.nclasses) ' classes'])
%     
%     set(handles.ModePopupMenu,'Value',3);
%     
%     ModePopupMenu_Callback(handles.ModePopupMenu,eventdata, handles);
%     handles.mode = 'class-view';
%     handles.submode = 'select';
    
end
guidata(gcbo,handles);


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SaveItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SaveButton_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function LoadItem_Callback(hObject, eventdata, handles)
% hObject    handle to LoadItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PrecomputeButton_Callback(hObject, eventdata, handles)

%guidata(gcbo,handles)

% --------------------------------------------------------------------
function HelpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to HelpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function HelpItem_Callback(hObject, eventdata, handles)
% hObject    handle to HelpItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

web('classify_spectra_help.html');

% --------------------------------------------------------------------
function AboutItem_Callback(hObject, eventdata, handles)
% hObject    handle to AboutItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('Classify_spectra (version 0.2) is being developed by the Mitra Lab at the Cold Spring Harbor Laboratory.','About classify_spectra');


% --------------------------------------------------------------------
function ConfigureItem_Callback(hObject, eventdata, handles)
% hObject    handle to ConfigureItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ConfigureButton_Callback(hObject, eventdata, handles)

function coefs = cepcoefs(data,p)
% Computes the cepstral coefficients from the data returns p coefficients
coefs = [];
if p > 0
    y = fft(hamming(length(data)) .* data);
    coefs = aryule(ifft(log(abs(y))),p);
end

function generate_output(handles,absolutetime,filename);
    
    filebasetosave = filename(1:length(filename)-4);

    matrixoutput = cell(handles.nsegments,1);

    for i = 1:handles.nsegments
        segment = handles.segments(i);
        
        if absolutetime
            dstr = regexp(segment.wavfile,'[0-9]+\-[0-9]+\-[0-9]+','match');
            dstr = dstr{1}; % take first match only
            tstr = regexp(segment.wavfile,'[0-9][0-9][0-9][0-9][0-9][0-9]','match');
            tstr = tstr{1};
            tstr = [tstr(1:2) ':' tstr(3:4) ':' tstr(5:6)];
            segmentstart = datenum([dstr ' ' tstr]) + segment.start;
            segmentstart = datevec(segmentstart);
        else
            segmentstart = segment.start;
        end
            matrixoutput{i} = {segment.wavfile,segment.class,segmentstart,segment.end - segment.start,segment.features};
    end
    
    save([filebasetosave '.mat'],'-mat','matrixoutput','-mat');
    
    delimiter = '\t';
    fp = fopen([filebasetosave '.txt'],'wt');
    
    % Generate header for file
    
    fprintf(fp,'filename');fprintf(fp,delimiter);
    fprintf(fp,'class');fprintf(fp,delimiter);
    
    if absolutetime
        fprintf(fp,'year');fprintf(fp,delimiter);
        fprintf(fp,'month');fprintf(fp,delimiter);
        fprintf(fp,'day');fprintf(fp,delimiter);
        fprintf(fp,'hour');fprintf(fp,delimiter);
        fprintf(fp,'minute');fprintf(fp,delimiter);
        fprintf(fp,'second');fprintf(fp,delimiter);
    else
        fprintf(fp,'start');fprintf(fp,delimiter);
    end
    
    fprintf(fp,'length');fprintf(fp,delimiter);
    
    for i = 1:handles.nfeatures
        fprintf(fp,['d' num2str(i)]);
        if i < handles.nfeatures;
            fprintf(fp,delimiter);
        else
            fprintf(fp,'\n');
        end
    end
    
    
    % Generate main data
   
    for i = 1:handles.nsegments
        segment = handles.segments(i);
        fprintf(fp,['"' segment.wavfile '"']);fprintf(fp,delimiter);
        fprintf(fp,['"' segment.class '"']);fprintf(fp,delimiter);
        
        if absolutetime
            time = matrixoutput{i}{3};
            
            fprintf(fp,num2str(time(1)));fprintf(fp,delimiter);
            fprintf(fp,num2str(time(2)));fprintf(fp,delimiter);
            fprintf(fp,num2str(time(3)));fprintf(fp,delimiter);
            fprintf(fp,num2str(time(4)));fprintf(fp,delimiter);
            fprintf(fp,num2str(time(5)));fprintf(fp,delimiter);
            fprintf(fp,num2str(time(6)));

        else
            fprintf(fp,num2str(segment.start));
        end
        
        fprintf(fp,delimiter); 
        fprintf(fp,num2str(segment.end-segment.start));fprintf(fp,delimiter);
        
        for j = 1:handles.nfeatures
            fprintf(fp,num2str(segment.features(j)));
            if j < handles.nfeatures % handle eol formatting
                fprintf(fp,delimiter);
            else
                if i < handles.nsegments % handle eof formatting
                    fprintf(fp,'\n');    
                end
            end
        end
    end
    fclose(fp);
% --------------------------------------------------------------------
function ExportDataItem_Callback(hObject, eventdata, handles)
% hObject    handle to ExportDataItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = uiputfile('*.txt','File to write exported data to');

if not(filename == 0)
    if length(handles.segments(1).features) == 0
        handles = features_segments(handles);
    end
    generate_output(handles,1,filename);
end


% --------------------------------------------------------------------
function CleanDirectoryItem_Callback(hObject, eventdata, handles)
% hObject    handle to CleanDirectoryItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answer = questdlg('Clean the current directory. This will remove all files generated by classify_spectra including the file which includes the classification');
if strcmp(answer,'Yes');
    specfiles = dir('*.spec');
    
    for i = 1:length(specfiles)
        delete(specfiles(i).name);
    end
    
    
    if exist('class_spec.conf')
        delete('class_spec.conf');
    end
    
    if exist([handles.baseclassname '.dat'])
        delete([handles.baseclassname '.dat']);
    end
   
    
end