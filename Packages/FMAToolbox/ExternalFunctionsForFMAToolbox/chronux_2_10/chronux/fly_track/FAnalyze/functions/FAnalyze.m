function varargout = FAnalyze(varargin)
% FANALYZE 
% For all your trajectory analysis needs! . See documentation for usage details.

%Written by Dan Valente
%November 2007

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FAnalyze_OpeningFcn, ...
                   'gui_OutputFcn',  @FAnalyze_OutputFcn, ...
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


% --- Executes just before FAnalyze is made visible.
function FAnalyze_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FAnalyze (see VARARGIN)

% Choose default command line output for FAnalyze
handles.output = hObject;
handles.called = 0;
disp('Welcome to FAnalyze!')
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FAnalyze wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FAnalyze_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ws_vars.
function ws_vars_Callback(hObject, eventdata, handles)
% hObject    handle to ws_vars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ws_vars contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ws_vars


% --- Executes during object creation, after setting all properties.
function ws_vars_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ws_vars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bins_Callback(hObject, eventdata, handles)
% hObject    handle to bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bins as text
%        str2double(get(hObject,'String')) returns contents of bins as a double


% --- Executes during object creation, after setting all properties.
function bins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in view_dist.
function view_dist_Callback(hObject, eventdata, handles)
% hObject    handle to view_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.called = handles.called+1;
index_selected = get(handles.ws_vars,'Value');
bins = get(handles.bins,'String');
phase_opt = get(get(handles.phase_space,'SelectedObject'), 'Tag');

if (length(index_selected) == 1)
   
        var1 = get_var_names(handles);
        u = strfind(var1,'_');
        cur_zone1 = var1(u(1)+1:u(2)-1);
        cur_zone1 = strmatch(cur_zone1, handles.zone_names,'exact');
        cur_seg1 = var1(u(2)+1:end);
        cur_seg1 = strmatch(cur_seg1, handles.zone{cur_zone1}.seg_label,'exact');
        cur_var1 = var1(1:strfind(var1,'_')-1);

        if strcmp(phase_opt,'phase1D')
            P.label = var1;
            P.phase_opt = phase_opt;
            [P.data P.bins] = eval(['ProbDist1D(handles.zone{' num2str(cur_zone1) '}.' cur_var1 ...
                '{' num2str(cur_seg1) '}.data ,' bins ')']);
        elseif strcmp(phase_opt,'phase2D')
            P.label = var1;
            P.phase_opt = phase_opt;
            [P.data P.bins] = eval(['ProbDist2D(handles.zone{' num2str(cur_zone1) '}.' cur_var1...
                '{' num2str(cur_seg1) '}.data ,' bins ')']);
        end

        %plot the distribution
        figure
        plot(P.bins,P.data)
    
    
elseif (length(index_selected) == 2)
    
    % should only look at joint distribution for same speed segments in
    % same zones right now.   The JointDist command expects equal length
    % vectors.  The following variables are just calculated in case we
    % modify JointDist down the road to handle vectors of different
    % lengths.
    
    if numel(str2num(bins)) == 2 | numel(str2num(bins))== 0
        [var1 var2] = get_var_names(handles);
        u1 = strfind(var1,'_');
        u2 = strfind(var2,'_');
        cur_zone1 = var1(u1(1)+1:u1(2)-1);
        cur_zone1 = strmatch(cur_zone1, handles.zone_names,'exact');
        cur_zone2 = var2(u2(1)+1:u2(2)-1);
        cur_zone2 = strmatch(cur_zone2, handles.zone_names,'exact');
        cur_seg1 = var1(u1(2)+1:end);
        cur_seg1 = strmatch(cur_seg1, handles.zone{cur_zone1}.seg_label,'exact');
        cur_seg2 = var2(u2(2)+1:end);
        cur_seg2 = strmatch(cur_seg2, handles.zone{cur_zone2}.seg_label,'exact');
        cur_var1 = var1(1:strfind(var1,'_')-1);
        cur_var2 = var2(1:strfind(var2,'_')-1);

        P.label = [var1 '_' var2];
        P.phase_opt = 'N/A';
        [P.data P.bins] = eval(['JointDist(handles.zone{' num2str(cur_zone1) '}.' cur_var1...
                '{' num2str(cur_seg1) '}.data , handles.zone{' num2str(cur_zone2) '}.' cur_var2...
                '{' num2str(cur_seg2) '}.data ,' bins ')']);

        figure
        imagesc(P.bins{2},P.bins{1},log(P.data))
    else
        errordlg('For a joint distribution, you must define bins or number of bins for BOTH directions!','Missing Bin Parameter')
        return;
    end
elseif length(index_selected) > 2
        errordlg('You must select only one or two variables, no more than two!!!',...
                'Incorrect Selection','modal')
end

if exist('P')
handles.P{handles.called} = P;
else
    handles.P{handles.called} = [];
end

guidata(gcbo,handles);


% --- Executes on button press in segment_speed.
function segment_speed_Callback(hObject, eventdata, handles)
% hObject    handle to segment_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

noise_thresh = str2num(get(handles.noise_thresh,'String'));
NumZones = length(handles.zone_names);
list1 = {'x','y','r','theta','vx','vy','v','vtheta','tau','kappa','beta'};
fps = handles.fps;

disp('Segmenting speed.  Input parameters and wait...')
%based on how many spatial zones, ask how we want to segment speed.
%Remember that first zone is always the full arena.
n=0;
for i=1:NumZones
    prompt1 = {['How many speed segments in ' handles.zone{i}.zone_label{1}]};
    dlg_title1 = 'Segment Speed';
    num_lines1 = 1;
    def1 = {'2'};
    num_segs = str2double(inputdlg(prompt1,dlg_title1,num_lines1,def1));
    nameB4 = 'dummy';
    seg_names = [];
    speed_thresh = noise_thresh;
    for j=1:num_segs-1;
        prompt2 = {['Threshold between Segments ' num2str(j), ' and ' num2str(j+1) '.'],...
            'Names for these zones:', ' '};
        dlg_title2 = 'Please Input Speed Segmentation Data';
        num_lines2 = 1;
        def2 = {'0.75','NZS','FSS'};
        answer = inputdlg(prompt2,dlg_title2,num_lines2,def2);
        speed_thresh = [speed_thresh str2double(answer(1))];
        nameA = answer(2);
        nameB = answer(3);
        check = strcmp(nameA,nameB4);
        if (j~=1 & ~check)
            errordlg([handles.zone{i}.zone_label{1} ' Segment ' num2str(j)...
                ' name must be same as previous ' handles.zone{i}.zone_label{1} ' Segment '...
                 num2str(j) ' name! Rename!'])
             return;
        else
            seg_names = [seg_names nameA nameB];
        end
        nameB4 = nameB;
    end
    
    if (num_segs == 1 | num_segs == 2)
        seg_names = seg_names;
    elseif (num_segs== 3)
        seg_names = [seg_names(1:2) seg_names(end)];
    else
        seg_names = [seg_names(1:2) seg_names(end-1:end)];
    end
    %check speed thresholds in that zone
    [s indx] = sort(speed_thresh,'ascend');
    if (indx ~= [1:length(speed_thresh)])
        errordlg(['Speed thresholds for successive segments should be larger than previous segments.  Try Again!'])
             return;
    elseif (~isempty(find(s < noise_thresh)))
        errordlg(['None of the thresholds are permitted to be below the noise threshold. Try Again!'])
             return;
    end
    
    seg_names = [{'all'} {'stops'} seg_names];
    
  
        for q=1:length(seg_names)
            for  j=1:length(list1)
            if (strcmp(seg_names{q},'all') & strcmp(list1{j},'kappa')) | ...
                    (strcmp(seg_names{q},'all') & strcmp(list1{j},'tau'))
               temp_list{j,q+n*length(seg_names)} = [];
            elseif (strcmp(seg_names{q},'stops') & strcmp(list1{j},'kappa')) | ...
                    (strcmp(seg_names{q},'stops') & strcmp(list1{j},'beta'))
               temp_list{j,q+n*length(seg_names)} = [];
            elseif ~strcmp(seg_names{q},'all') & strcmp(list1{j},'beta')
                    temp_list{j,q+n*length(seg_names)} = [];
            else
                temp_list{j,q+n*length(seg_names)} = [list1{j} '_' handles.zone{i}.zone_label{1} '_' seg_names{q}];
            end
        end
    end
    
    %calculate stuff for this particular zone
    t = handles.zone{i}.t{1}.data;
    x = handles.zone{i}.x{1}.data;
    y = handles.zone{i}.y{1}.data;
    r = handles.zone{i}.r{1}.data;
    theta = handles.zone{i}.theta{1}.data;
    vx = handles.zone{i}.vx{1}.data;
    vy = handles.zone{i}.vy{1}.data;
    v = handles.zone{i}.v{1}.data;
    vtheta =handles.zone{i}.vtheta{1}.data;
      
    
    stops_indx = find(v < speed_thresh(1));
    moves_indx = find(v >= speed_thresh(1));
    
    v(stops_indx) = 0; vx(stops_indx) = 0; vy(stops_indx) = 0; vtheta = atan2(vy,vx);
    
    
    handles.zone{i}.seg_label = seg_names;
    handles.zone{i}.t{2}.data = t(stops_indx);
    handles.zone{i}.x{2}.data = x(stops_indx);
    handles.zone{i}.y{2}.data = y(stops_indx);    
    handles.zone{i}.r{2}.data = r(stops_indx);    
    handles.zone{i}.theta{2}.data = theta(stops_indx);    
    handles.zone{i}.vx{2}.data = vx(stops_indx);    
    handles.zone{i}.vy{2}.data = vy(stops_indx);   
    handles.zone{i}.v{2}.data = v(stops_indx);    
    handles.zone{i}.vtheta{2}.data = vtheta(stops_indx);
    handles.zone{i}.tau{2}.data = FindDuration(stops_indx);
    handles.zone{i}.kappa{2}.data = 'N/A';
    
    
    for j=2:length(speed_thresh)
        
        indxA =  find(v >= speed_thresh(j-1) & v < speed_thresh(j));
        indxB =  find(v >= speed_thresh(j));
        
        handles.zone{i}.t{j+1}.data = t(indxA);
        handles.zone{i}.x{j+1}.data = x(indxA);     
        handles.zone{i}.y{j+1}.data = y(indxA);        
        handles.zone{i}.r{j+1}.data = r(indxA);       
        handles.zone{i}.theta{j+1}.data = theta(indxA);       
        handles.zone{i}.vx{j+1}.data = vx(indxA);       
        handles.zone{i}.vy{j+1}.data = vy(indxA);     
        handles.zone{i}.v{j+1}.data = v(indxA);      
        handles.zone{i}.vtheta{j+1}.data = vtheta(indxA);
        handles.zone{i}.tau{j+1}.data = FindDuration(indxA);
        handles.zone{i}.kappa{j+1}.data = CalcCurvature(v, vtheta, indxA, 1/fps);
        
        handles.zone{i}.t{j+2}.data = t(indxB);
        handles.zone{i}.x{j+2}.data = x(indxB);    
        handles.zone{i}.y{j+2}.data = y(indxB);     
        handles.zone{i}.r{j+2}.data = r(indxB);       
        handles.zone{i}.theta{j+2}.data = theta(indxB);      
        handles.zone{i}.vx{j+2}.data = vx(indxB);     
        handles.zone{i}.vy{j+2}.data = vy(indxB);       
        handles.zone{i}.v{j+2}.data = v(indxB);    
        handles.zone{i}.vtheta{j+2}.data = vtheta(indxB);
        handles.zone{i}.tau{j+2}.data = FindDuration(indxB);
        handles.zone{i}.kappa{j+2}.data = CalcCurvature(v, vtheta, indxB, 1/fps);

    end
    handles.zone{i}.beta{1}.data = CalcReorientAngle(vtheta, moves_indx);
    n=n+1;
end


temp3 = {};
sz = size(temp_list);
for q=1:sz(2)
    temp3 = [temp3 temp_list{:,q}];
end

disp('Speed has been segmented.')
update_listbox(handles, temp3)
guidata(gcbo,handles);

% --- Executes on button press in load_traj.
function load_traj_Callback(hObject, eventdata, handles)
% hObject    handle to load_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile({'*.mat'}, 'Select the .mat file containing the trajectory', 'MultiSelect','off');

if isequal(filename,0) || isequal(pathname,0)
    disp('File select canceled')
    return;
else
    disp(['File selected: ', fullfile(pathname, filename)])
    load(fullfile(pathname,filename));
    handles.filename = fullfile(pathname, filename); 
    
end


handles.zone{1}.t{1}.data = t;
handles.fps = 1/(t(2)-t(1));
handles.zone{1}.x{1}.data = x;
handles.zone{1}.y{1}.data = y;

guidata(gcbo,handles);

% --- Executes on button press in smooth.
function smooth_Callback(hObject, eventdata, handles)
% hObject    handle to smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%grab data 
t = handles.zone{1}.t{1}.data;
fps = handles.fps;
x = handles.zone{1}.x{1}.data;
y = handles.zone{1}.y{1}.data;
handles.P = [];

%grab smoothing parameters as input by the user
n = str2double(get(handles.n,'String'));
dn = str2double(get(handles.dn,'String'));

disp('Smoothing data.  Please wait...')
%smooth using runline
x = runline(x, n, dn);
y = runline(y, n, dn);

%Calculate polar coords and velocity.
r = sqrt(x.^2+y.^2);
theta = atan2(y,x);

vx = fps*gradient(x);
vy = fps*gradient(y);
v = sqrt(vx.^2+vy.^2);
vtheta = atan2(vy,vx);

clear handles.x;
clear handles.y;
%save smooth data to handles to be used by other GUI functions
handles.zone_names{1} = 'Full Arena';
handles.zone{1}.zone_label = {'Full Arena'};
handles.zone{1}.seg_label = {'all'};
handles.zone{1}.t{1}.data = t;
handles.zone{1}.x{1}.data = x;
handles.zone{1}.y{1}.data = y;

handles.zone{1}.r{1}.data = r;
handles.zone{1}.theta{1}.data = theta;
handles.zone{1}.vx{1}.data = vx;
handles.zone{1}.vy{1}.data = vy;
handles.zone{1}.v{1}.data = v;
handles.zone{1}.vtheta{1}.data = vtheta;
handles.zone{1}.tau{1}.data = 'N/A';
handles.zone{1}.kappa{1}.data = 'N/A';
    


disp('Trajectory has been smoothed.')

update_listbox(handles, {'x_Full Arena_all','y_Full Arena_all','r_Full Arena_all','theta_Full Arena_all',...
    'vx_Full Arena_all','vy_Full Arena_all','v_Full Arena_all','vtheta_Full Arena_all'});
guidata(gcbo,handles);

function n_Callback(hObject, eventdata, handles)
% hObject    handle to n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n as text
%        str2double(get(hObject,'String')) returns contents of n as a double


% --- Executes during object creation, after setting all properties.
function n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dn_Callback(hObject, eventdata, handles)
% hObject    handle to dn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dn as text
%        str2double(get(hObject,'String')) returns contents of dn as a double


% --- Executes during object creation, after setting all properties.
function dn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in view_traj.
function view_traj_Callback(hObject, eventdata, handles)
% hObject    handle to view_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

t = handles.zone{1}.t{1}.data;
x = handles.zone{1}.x{1}.data;
y = handles.zone{1}.y{1}.data;
tmp = getfield(handles, 'zone');

if (isfield(tmp{1},'r'))
    r = handles.zone{1}.r{1}.data;
    theta = handles.zone{1}.theta{1}.data;
    vx = handles.zone{1}.vx{1}.data;
    vy = handles.zone{1}.vy{1}.data;
    v = handles.zone{1}.v{1}.data;
    vtheta =handles.zone{1}.vtheta{1}.data; 
end

str = get(handles.variables, 'String');
val = get(handles.variables, 'Value');

% Set current data to the selected data set.
      if strcmp(str(val),'(x,y)')
          figure
          plot(x,y)
          xlabel('x-position')
          ylabel('y-position')
          title('Trajectory')
          if (~isfield(tmp{1},'r'))
              disp('You can view your trajectory, but you MUST smooth your data before proceeding with any calculations!')
          end
      else
          
          var1 = 't';
          if exist(str{val})
              figure
              var2 = str{val};
              eval(['plot(' var1 ','  var2 ')'])
              xlabel('Time')
              ylabel(var2)
              xlim([0 t(end)])
          else
              errordlg('Can"t view this variable unless data is smoothed','Non-existant Variable')
             
              return;
          end
      end

% --- Executes on selection change in variables.
function variables_Callback(hObject, eventdata, handles)
% hObject    handle to variables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns variables contents as cell array
%        contents{get(hObject,'Value')} returns selected item from variables

set(hObject,'String',{'(x,y)','x','y','r','theta','vx','vy','v','vtheta'})

guidata(gcbo,handles);

% --- Executes during object creation, after setting all properties.
function variables_CreateFcn(hObject, eventdata, handles)
% hObject    handle to variables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String',{'(x,y)','x','y','r','theta','vx','vy','v','vtheta'})

guidata(gcbo,handles);
% --- Executes on button press in segment_space.
function segment_space_Callback(hObject, eventdata, handles)
% hObject    handle to segment_space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NumZones = str2double(get(handles.num_zones,'String')); 
if NumZones < 2
    errordlg('Must have at least 2 zones','Too few zones')
    return;
end

zone_indx = 1:NumZones;
list1 = {'x','y','r','theta','vx','vy','v','vtheta','tau','kappa'};
zone_names = handles.zone_names;
handles.zone_names = {};
handles.zone_names{1} = 'Full Arena';

t = handles.zone{1}.t{1}.data;
dt = t(2)-t(1);
x = handles.zone{1}.x{1}.data;
y = handles.zone{1}.y{1}.data;
r = handles.zone{1}.r{1}.data;
theta = handles.zone{1}.theta{1}.data;
vx = handles.zone{1}.vx{1}.data;
vy = handles.zone{1}.vy{1}.data;
v = handles.zone{1}.v{1}.data;
vtheta =handles.zone{1}.vtheta{1}.data;

disp('Segmenting Space.  Input parameters and wait...')

thresh(1) = 0;
for i=1:NumZones-1
    prompt = {['Threshold between Zones ' num2str(zone_indx(i)), ' and ' num2str(zone_indx(i+1)) '.'],...
        'Names for these zones:', ' '};
    dlg_title = 'Please Input Spatial Zone Data';
    num_lines = 1;
    def = {'7.3','CZ','RZ'};
    answer{i} = inputdlg(prompt,dlg_title,num_lines,def);
    thresh(i+1) = str2double(answer{i}(1));
end

zone_names = [zone_names answer{1}(2)];
for i=1:NumZones-2
    nameA = answer{i}(2);
    nameB = answer{i}(3);
    nameC = answer{i+1}(2);
    if (~strcmp(nameB,nameC))
        errordlg(['Zone ' num2str(i+1) ' name must be same as previous Zone ' num2str(i+1) ' name! Rename!'])
    else
        zone_names = [zone_names nameB];
    end
end    
zone_names = [zone_names answer{NumZones-1}(3)];

for i=1:length(zone_names)
    for  j=1:length(list1)
        if  strcmp(zone_names{i},'Full Arena') & (strcmp(list1{j},'kappa') | strcmp(list1{j},'tau'))
            temp_list{j,i} = [];
        elseif strcmp(list1{j},'kappa')
            temp_list{j,i} = [];
        else
        temp_list{j,i} = [list1{j} '_' zone_names{i} '_all'];
        end
    end
end

temp3 = {};
for i=1:NumZones+1
    temp3 = [temp3 temp_list{:,i}];
end

handles.zone_names = zone_names;

%Now using thresholds, divvy up space
for i = 2:length(thresh)
    
        % Can only do radial zones in this version of FAnalyze
        indxA =  find(r >= thresh(i-1) & r < thresh(i));
        indxB =  find(r >= thresh(i));
        
        handles.zone{i}.zone_label = zone_names(i);
        handles.zone{i}.seg_label = {'all'};
        handles.zone{i}.t{1}.data = t(indxA);
        handles.zone{i}.x{1}.data = x(indxA);
        handles.zone{i}.y{1}.data = y(indxA);
        handles.zone{i}.r{1}.data = r(indxA);
        handles.zone{i}.theta{1}.data = theta(indxA);
        handles.zone{i}.vx{1}.data = vx(indxA);
        handles.zone{i}.vy{1}.data = vy(indxA);
        handles.zone{i}.v{1}.data = v(indxA);
        handles.zone{i}.vtheta{1}.data = vtheta(indxA);
        handles.zone{i}.tau{1}.data = FindDuration(indxA);
        handles.zone{i}.kappa{1}.data = CalcCurvature(v, vtheta, indxA, dt);
        
        handles.zone{i+1}.zone_label = zone_names(i+1);
        handles.zone{i+1}.seg_label = {'all'};
        handles.zone{i+1}.t{1}.data = t(indxB);
        handles.zone{i+1}.x{1}.data = x(indxB);
        handles.zone{i+1}.y{1}.data = y(indxB);
        handles.zone{i+1}.r{1}.data = r(indxB);
        handles.zone{i+1}.theta{1}.data = theta(indxB);
        handles.zone{i+1}.vx{1}.data = vx(indxB);
        handles.zone{i+1}.vy{1}.data = vy(indxB);
        handles.zone{i+1}.v{1}.data = v(indxB);
        handles.zone{i+1}.vtheta{1}.data = vtheta(indxB);
        
        
        handles.zone{i+1}.tau{1}.data = FindDuration(indxB);
        handles.zone{i+1}.kappa{1}.data = CalcCurvature(v, vtheta, indxB, dt);
      
 
       
end

update_listbox(handles, temp3)
disp('Space has been segmented.')
guidata(gcbo,handles);

% --- Executes on button press in save_ws.
function save_ws_Callback(hObject, eventdata, handles)
% hObject    handle to save_ws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% %This whole function just defines variables and saves the workspace as a
% %.mat file.
[filename, pathname] = uiputfile('*.mat', 'Pick a MAT file to save to');

traj = handles.zone;
P = handles.P;
fps = handles.fps;
save(fullfile(pathname, filename),'traj','P')       
disp(['Data has been saved to ' fullfile(pathname, filename)])   

function update_listbox(handles, vars)
% this function updates the message center at the bottom of the GUI
% adapted from Mike Rieser's PControl GUI.

set(handles.ws_vars, 'String', vars);

%%%%%%%%%%%%%%%%
function varargout = get_var_names(handles)
list_entries = get(handles.ws_vars,'String');
index_selected = get(handles.ws_vars,'Value');

varargout = list_entries(index_selected);





function num_zones_Callback(hObject, eventdata, handles)
% hObject    handle to num_zones (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_zones as text
%        str2double(get(hObject,'String')) returns contents of num_zones as a double


% --- Executes during object creation, after setting all properties.
function num_zones_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_zones (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noise_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to noise_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise_thresh as text
%        str2double(get(hObject,'String')) returns contents of noise_thresh as a double




% --- Executes during object creation, after setting all properties.
function noise_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


