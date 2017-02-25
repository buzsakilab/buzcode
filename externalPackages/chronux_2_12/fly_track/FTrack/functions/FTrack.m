function varargout = FTrack(varargin)
% FTRACK 
% For all your fly-tracking needs! . See documentation for usage details.

% Last Modified by GUIDE v2.5 26-Nov-2007 18:07:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FTrack_OpeningFcn, ...
                   'gui_OutputFcn',  @FTrack_OutputFcn, ...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before FTrack is made visible.
function FTrack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FTrack (see VARARGIN)

% Choose default command line output for FTrack
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
disp('Welcome to FTrack!')
warning off all

% UIWAIT makes FTrack wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Outputs from this function are returned to the command line.
function varargout = FTrack_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in load_video.
function load_video_Callback(hObject, eventdata, handles)
% hObject    handle to load_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA


%The following allows the user to select multiple videos to be tracked.  
%The videos will be tracked sequentially, not in parallel. The filenames
%will be displayed in the message center.

[filename, pathname] = uigetfile({'*.avi;*.mpg;*.mp2','Video Files (*.avi,*.mpg,*.mp2)'}, 'Pick a video', 'MultiSelect','on');
if iscell(filename)
    NFiles = length(filename);
else
    NFiles = 1;
    filename = {filename};
end

if isequal(filename,0) || isequal(pathname,0)
    disp('File select canceled')
else
    for i = 1:NFiles
        disp(['Video selected: ', fullfile(pathname, filename{i})])
        handles.filename{i} = fullfile(pathname, filename{i}); 
    end
end

guidata(gcbo,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in viewframe.
function viewframe_Callback(hObject, eventdata, handles)
% hObject    handle to viewframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This is the single frame viewer.  It simply gives the user some info about
%the video (displayed in the message center), and also shows the user a
%single frame of the movie so that they can make sure it looks correct and
%they can get some parameters off of the frame if need be (such as fly size
%or a pixel/cm calilbration).

filename = handles.filename;
if iscell(filename)
    NFiles = length(filename);
else
    NFiles = 1;
    filename = {filename};
end

for i = 1:NFiles
    video = videoReader(filename{i});
    seek(video,1);
    info = getinfo(video);
    img = getframe(video);
    figure
    imagesc(img)
    title(filename{i})
    disp(['Video name:',filename{i}])
    disp(['Video dimensions [Width, Height]:  [', num2str(info.width),' , ' num2str(info.height),']'])
    disp(['Number of frames:  ',num2str(info.numFrames)])
    disp(['Frame rate:  ',num2str(info.fps),' frames/s'])
    
    handles.VideoInfo(i) = info;
end

disp('Frame Viewer Finished')
guidata(gcbo, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in arena_dims.
function arena_dims_Callback(hObject, eventdata, handles)
% hObject    handle to arena_dims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = handles.filename;
if iscell(filename)
    NFiles = length(filename);
else
    NFiles = 1;
    filename = {filename};
end

p = [];
q = [];
for i = 1:NFiles
    video = videoReader(filename{i});
    info = getinfo(video);
    seek(video,1);
    img = getframe(video);
    sz = size(img);
    figure
    colormap gray
    imagesc(img)
    axis image
    title({filename{i}; 'Left click to select points on boundary of arena. Press return when done.'})
 
    % Choose points on arena boundary
    [p q] = ginput;

    %Find center and radius of arena

    ellipse = fit_ellipse(p',q','y');
    semiminor = min(ellipse.a,ellipse.b);
    semimajor = max(ellipse.a,ellipse.b);
    ellipse.epsilon = sqrt(1-semiminor^2/semimajor^2);
    ellipse.psi = asin(ellipse.epsilon);
    ellipse.semiminor = semiminor;
    ellipse.semimajor = semimajor;

    rotated_ellipse = ellipse.rotated_ellipse;
    ellipse.points_selected = [p q];
    title('Arena w/ Boundary')
    ellipse.boundaries = [ rotated_ellipse(1,:)', info.height-rotated_ellipse(2,:)'];
       
    info = handles.VideoInfo(i);
    figure
    title('Masking arena.  Please wait...')
    drawnow
    disp('Masking arena...')
    mask = zeros(info.height, info.width);
    for k=1:info.height
        for j=1:info.width
            rot = inv(ellipse.R)*[j;k];
            testpoint = (rot(1)-ellipse.X0).^2/(ellipse.a).^2+(rot(2)-ellipse.Y0).^2/(ellipse.b).^2;
            if (testpoint <= 1.01) %giving a little lee-way near boundary
                mask(k,j) = 1;
            end
        end
    end
    ellipse.mask =  double(mask);
    handles.arena{i} = ellipse;
    
    imagesc(double(img(:,:,1)).*double(mask))
    axis image
    colormap gray
    title('Portion to be tracked')
    
end

%setting things up to only consider pixels within selected arena...

disp('Arena Dimensions Calculated')
disp(ellipse)
guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in input_params.
function input_params_Callback(hObject, eventdata, handles)
% hObject    handle to input_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Opens up a dialog box so that the user can input relevant parameters for
%each video that is being tracked.  A separate box will open for each
%video. 

filename = handles.filename;
if iscell(filename)
    NFiles = length(filename);
else
    Nfiles = 1;
    filename = {filename};
end
UserIn = [];

for i = 1:NFiles
    first = strvcat(filename{i}, '          ','Output directory');
    prompt = {first,'Start Frame (Remember that first frame is indexed 0)','End Frame:', 'Initial background size:', 'Arena radius (cm)', 'Bounding box half-size (in pixels)','Background weight (0.9 < a < 1)'};
    dlg_title = ['Input Paramters'];
    num_lines = 1;
    defaults = {'C:\Documents and Settings\liam\My Documents\LabVIEW Data\FTrack output','0','999', '100','7.5','10','0.9'};
    options.Resize='on';
    options.WindowStyle='normal';
    answer = inputdlg(prompt,dlg_title,num_lines,defaults, options);
    UserIn = [UserIn answer];
    disp(['Input parameters have been entered for ', filename{i}])
    
    InputData(i).OutputPath = UserIn{1,i};
    InputData(i).StartFrame = str2double(UserIn(2,i));
    InputData(i).EndFrame = str2double(UserIn(3,i));
    InputData(i).NBackFrames = str2double(UserIn(4,i));
    InputData(i).ArenaRadius = str2double(UserIn(5,i));
    InputData(i).sqrsize = str2double(UserIn(6,i));
    InputData(i).alpha = str2double(UserIn(7,i));

end

handles.InputData = InputData;
 
guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radio Buttons to select fly finding option

% --- Executes during object creation, after setting all properties.
function find_opt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to find_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function find_opt_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to find_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Radio buttons to select whether fly is white on black or black on white

% --- Executes during object creation, after setting all properties.
function neg_opt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neg_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --------------------------------------------------------------------
function neg_opt_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to neg_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in TrackStart.
function TrackStart_Callback(hObject, eventdata, handles)
% hObject    handle to TrackStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This part of the calls FlyTracker to do the actual tracking.

%grab the filename and user input parameters 
filename = handles.filename;
InputData = handles.InputData;

if iscell(filename)
    NVideos = length(filename);
else
    NVideos = 1;
    filename = {filename};
end

neg_opt = get(get(handles.neg_opt,'SelectedObject'), 'Tag');

% Call FlyTracker
for i=1:NVideos
     
    %video parameters, just because...
    info = handles.VideoInfo(i);
    FrameRate = info.fps;
    FrameRange = [InputData(i).StartFrame:InputData(i).EndFrame];
    
    %Track the fly!
    [x, y, orientation] = FlyTracker(filename{i}, FrameRange,...
                       InputData(i).NBackFrames, neg_opt, InputData(i).sqrsize,...
                        InputData(i).alpha, handles.arena{i});
                    
    %Time vector
    if (info.numFrames == InputData(i).EndFrame)
        InputData(i).EndFrame = InputData(i).EndFrame-1;
    elseif (info.numFrames < InputData(i).EndFrame)
        InputData(i).EndFrame = info.numFrames-1;
    end
    
    t = [InputData(i).StartFrame:InputData.EndFrame(i)]/FrameRate;
        
    %Save data to handles so the rest of the GUI can access it.  
    handles.x = x;
    handles.y = y;
    handles.orientation = orientation;
    handles.t = t;

    guidata(gcbo, handles);  
    %Also, save variables as a .mat file so that the user will be able to
    %load in the tracked data later on.
    
    out = SaveFiles(i,handles);
    
end
disp('Tracking Complete.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in view_traj.
function view_traj_Callback(hObject, eventdata, handles)
% hObject    handle to view_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = ...
     uigetfile({'*.mat'},'Select Raw Trajectory Data');

if isequal(filename,0) || isequal(pathname,0)
    disp('File select canceled')
    return;
else
    fullname = fullfile(pathname, filename);
end


load(fullname)
figure
set(gcf,'Name',fullname)
subplot(2,2,1)
plot(t,x)
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('cm')
title('Raw x data')
subplot(2,2,3)
plot(t,y)
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('cm')
title('Raw y data')
subplot(2,2,[2 4])
plot(x,y)
xlabel('x position (cm)')
ylabel('y position (cm)')
title('Raw Trajectory')
axis equal

handles.filename = filename;
handles.x = x;
handles.y = y;
handles.orientation = orientation;
handles.t = t;
handles.VideoInfo = VideoInfo;
handles.InputData = InputData;
handles.arena{1} = arena;

guidata(gcbo, handles);  
     



% --- Executes on button press in clean_x.
function clean_x_Callback(hObject, eventdata, handles)
% hObject    handle to clean_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = handles.x;
filename = handles.filename;

%open up plot to examine data
figure
plot(x)
title('Zoom in if necessary.  Press any key to continue.')
xlabel('Frame')
ylabel('cm')
pause

% This loop runs until the user hits return to exit. We wait for the user
% to click four points on the graph, and then run the CleanData function on
% the region defined by those four points.

while(1)
title('Define region of data to clean: left, right, baseline, threshold. Hit Return to exit.')
[p q] = ginput(4);
if isempty(p)
    close;
    disp('X data has been cleaned.')
    out = SaveFiles(1,handles);
    return;
else
    rnge = floor(p(1)):floor(p(2));
end

if (q(3) > q(4))
    choice = 'below';
elseif (q(3) < q(4));
    choice = 'above';
end
epsilon = q(4);

x = CleanData(x, rnge, choice, epsilon);
handles.x = x;
guidata(gcbo, handles)
ax = gca;
xlim_temp = get(ax, 'XLim');
ylim_temp = get(ax, 'YLim');

%Show plot of clean data
figure
plot(x)
xlim(xlim_temp);
ylim(ylim_temp);
xlabel('Frame')
ylabel('cm')
title('Zoom in if necessary.  Press any key to continue.')
pause
end



% --- Executes on button press in clean_y.
function clean_y_Callback(hObject, eventdata, handles)
% hObject    handle to clean_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

y = handles.y;
filename = handles.filename;
%open up plot to examine data
figure(2)
plot(y)
xlabel('Frame')
ylabel('cm')
title('Zoom in if necessary.  Press any key to continue.')
pause

% This loop runs until the user hits return to exit. We wait for the user
% to click four points on the graph, and then run the CleanData function on
% the region defined by those four points.

while(1)
title('Define region of data to clean: left, right, baseline, threshold. Hit Return to exit.')
[p q] = ginput(4);
if isempty(p)
    close;
    disp('Y data has been cleaned.')
    out = SaveFiles(1,handles);
    return;
else
    rnge = floor(p(1)):floor(p(2));
end

if (q(3) > q(4))
    choice = 'below';
elseif (q(3) < q(4));
    choice = 'above';
end
epsilon = q(4);

y = CleanData(y, rnge, choice, epsilon);
handles.y = y;
guidata(gcbo, handles)
ax = gca;
xlim_temp = get(ax, 'XLim');
ylim_temp = get(ax, 'YLim');

%Show plot of clean data
figure(2)
plot(y)
xlim(xlim_temp);
ylim(ylim_temp);
xlabel('Frame')
ylabel('cm')
title('Zoom in if necessary.  Press any key to continue.')
pause
end

return;

function out = SaveFiles(i,handles)
    
    if iscell(handles.filename)
        filename = handles.filename{i};
    else
        filename = handles.filename;
    end
    
    OutputPath = handles.InputData(i).OutputPath; 
    x = handles.x;
    y = handles.y;
    t = handles.t;
    orientation = handles.orientation;
    InputData = handles.InputData(i);
    VideoInfo = handles.VideoInfo(i);
    arena =  handles.arena{i};
    
    [temp, name, ext, versn] = fileparts(filename);
    
    name_mat = strcat(name,'.mat');
    name_xy=strcat(name,'.xy');
    name_ori=strcat(name,'.ori');
    save_filename = fullfile(OutputPath,name);
    save_filename_xy = fullfile(OutputPath,name_xy);
    save_filename_ori = fullfile(OutputPath,name_ori);
    ori=handles.orientation(1,:)';
    xy=[handles.x' handles.y'];
    save(save_filename,'x','y','t','orientation','InputData', 'VideoInfo','arena')
    disp(['Saved ',save_filename])
    save(save_filename_xy,'xy','-ascii','-tabs');
    disp(['Saved ',save_filename_xy])
    save(save_filename_ori,'ori','-ascii','-tabs');   
    disp(['Saved ',save_filename_ori])
    out = 1;
    
return;


% --- Executes on button press in tilt_correct.
function tilt_correct_Callback(hObject, eventdata, handles)
% hObject    handle to tilt_correct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Note.  This is a slightly different version of tilt correction than
%described in the PLoS paper.  I think this is more robust and general.

x = handles.x;
y = handles.y;
t = handles.t;
ellipse = handles.arena{1};
radius_in_cm = handles.InputData.ArenaRadius;
video = videoReader(handles.VideoInfo.url);
info = getinfo(video);
height = info.height;

psi = ellipse.psi;
phi = ellipse.phi;

%redefining coordinate system so center of arena is at (0,0)
rotated_ellipse(:,1) = ellipse.boundaries(:,1)-ellipse.X0_in;
rotated_ellipse(:,2) = ellipse.boundaries(:,2)-(height-ellipse.Y0_in);

M = [cos(phi) sin(phi);-sin(phi) cos(phi)];  %rotation matrix.

%stretch ellipse only in shortest dimension to length of longest.  
if (ellipse.a > ellipse.b)
    L = [1 0; 0 ellipse.a/ellipse.b]; 
else
    L = [ellipse.b/ellipse.a 0; 0 1];
end

T = L*M;  %transformation matrix.  A rotation and a stretch to turn an ellipse into a circle.

disp('Correcting for camera tilt.  This may take a few minutes, so be patient!')

z = zeros(1,length(rotated_ellipse(:,1)));

%Transform boundary line
for j = 1:length(rotated_ellipse(:,1))
    temp = T*[rotated_ellipse(j,1); rotated_ellipse(j,2)];
    xe(j) = temp(1);
    ye(j) = temp(2);      
end

%fit new ellipse 
z = zeros(1,length(x));
new_ellipse = fit_ellipse(xe',ye','n');

semiminor = min(new_ellipse.a,new_ellipse.b);
semimajor = max(new_ellipse.a,new_ellipse.b);
new_ellipse.epsilon = sqrt(1-semiminor^2/semimajor^2);
new_ellipse.psi = asin(new_ellipse.epsilon);
new_ellipse.semiminor = semiminor;
new_ellipse.semimajor = semimajor;
new_ellipse.boundaries = [xe' ye'];

%Now transform trajectory
for j = 1:length(x)   
        temp = T*[x(j)-ellipse.X0_in; y(j)-(height-ellipse.Y0_in)];
        xp(j) = temp(1);
        yp(j) = temp(2);              
end


if (abs(new_ellipse.b-new_ellipse.a) <= 2) %2 pixel accuracy seems adequate
    disp('Tilt corrected')
    pixels_in_cm = semimajor/radius_in_cm;

    x = (xp)/pixels_in_cm;
    y = (yp)/pixels_in_cm;
    
    figure
    plot(x,y,new_ellipse.boundaries(:,1)/pixels_in_cm,new_ellipse.boundaries(:,2)/pixels_in_cm,'r')
    title('Transformed trajectory with boundary')

    handles.arena{1}=new_ellipse;
    handles.x = x;
    handles.y = y;
    handles.InputData.pixels_in_cm =  pixels_in_cm;
    out = SaveFiles(1,handles);
else
    disp('Something went wrong. Arena is still an ellipse. ')
end


guidata(gcbo, handles)


 
 
 
 

