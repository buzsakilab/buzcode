function varargout = auto_classify(varargin)
% AUTO_CLASSIFY M-file for auto_classify.fig
%      AUTO_CLASSIFY, by itself, creates a new AUTO_CLASSIFY or raises the existing
%      singleton*.
%
%      H = AUTO_CLASSIFY returns the handle to a new AUTO_CLASSIFY or the handle to
%      the existing singleton*.
%
%      AUTO_CLASSIFY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTO_CLASSIFY.M with the given input arguments.
%
%      AUTO_CLASSIFY('Property','Value',...) creates a new AUTO_CLASSIFY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before auto_classify_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to auto_classify_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help auto_classify

% Last Modified by GUIDE v2.5 21-Jun-2006 00:34:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @auto_classify_OpeningFcn, ...
                   'gui_OutputFcn',  @auto_classify_OutputFcn, ...
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


% --- Executes just before auto_classify is made visible.
function auto_classify_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to auto_classify (see VARARGIN)

% input parameters

if isempty(varargin)
    handles.matrix2classify = [];
    handles.ncepstral = 5;
else
    handles.matrix2classify = varargin{1};
    handles.ncepstral = varargin{2};
end

handles.distancemeasure = 'sqEuclidean';

handles.cluster_method = 'kmeans';

set(handles.CepstralPopupMenu,'String',num2str([0:handles.ncepstral]'));
set(handles.CepstralPopupMenu,'Value',handles.ncepstral); %handles.ncepstral +1 );
% set(handles.CepstralPopupMenu,'Value',handles.ncepstral-1);

% Choose default command line output for auto_classify
% handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes auto_classify wait for user response (see UIRESUME)
% uiwait(handles.figure1);

uiwait(handles.figure1);
    
% --- Outputs from this function are returned to the command line.
function varargout = auto_classify_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

close;

function clusters = clusterdata(X,k,dmeasure)
% A wrapper function for clustering

try
    clusters = kmeans(X,k,'distance',dmeasure,'EmptyAction','singleton','replicates',30);
catch
    clusters = k;
end

% --- Executes on button press in AutoClassifyButton.
function AutoClassifyButton_Callback(hObject, eventdata, handles)
% hObject    handle to AutoClassifyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.AcceptButton,'Enable','off');
guidata(gcbo,handles);

cncepstral = get(handles.CepstralPopupMenu,'Value');

matrix2classify = handles.matrix2classify;
matrix2classify = madnormalize(matrix2classify,1); % normalize the columns
matrix2classify = matrix2classify(:,1:cncepstral+1); % include the cepstral coefficients that you want

if length(matrix2classify(1,:)) > 3
    matrix2classify = matrix2classify(:,[1,3:cncepstral+1]); % exclude the first cepestral coefficient
end

try
    segweight = str2num(get(handles.SegLenWeight,'String'));
catch
    segweight = 1;
    set(handles.SegLenWeight,'String','1');
end

% Weight the segment lengths by square root of the weight
matrix2classify(:,1) = matrix2classify(:,1)*sqrt(segweight);
 
% Other weightings can be added here but the cepestral coefficients have a
% natural exponential decline which acts a weighting

%Handle transformations to segment length

contents = get(handles.TransformPopupMenu,'String');
value = get(handles.TransformPopupMenu,'Value');

switch contents{value}
    case 'None'
        matrix2classify = matrix2classify;
    case 'Exclude'
        matrix2classify = matrix2classify(:,2:length(matrix2classify(1,:)));
    case 'Log'
        matrix2classify(:,1) = log(matrix2classify(:,1));
    case 'Exp'
        matrix2classify(:,1) = exp(matrix2classify(:,1));
end

if length(matrix2classify(1,:)) == 0 % matrix to classify has no information
    return
end

rangestate = get(handles.RangeSpecify,'Value');
diagnosticstate = get(handles.DiagnosticCheckbox,'Value');

if strcmp(handles.cluster_method,'hierarchical')
    handles.classification = cluster_hierarchical(matrix2classify);
    nclusters = length(unique(handles.classification)) - 1;
    set(handles.KClasses,'String',nclusters);
    
    if diagnosticstate
        pcaplot(matrix2classify,handles.classification);
    end
    
else if strcmp(handles.cluster_method, 'kmeans')

if not(rangestate)

    kclassstr = get(handles.KClasses,'String');

    if not(isempty(kclassstr))    
        try
            nclasses = str2num(kclassstr);
        catch
            return;
        end
    
        classification =  clusterdata(matrix2classify,nclasses,handles.distancemeasure);% this can be changed to a different clustering algorithm
    
        % clusterdata shows how a clustering algorithm can be hooked in where the number of
        % classes need to be assigned.
        
        if length(classification) == 1
            msgbox('Error in clustering. Try changing the number of target clusters or rerunning.');
        else
            if diagnosticstate
                figure();
                [silh h] = silhouette(matrix2classify,classification);
                pcaplot(matrix2classify,classification);
            end
        handles.classification = classification;
        end
    else
        return;
    end
else
    try
        minclust = str2num(get(handles.MinClusters,'String'));
    catch
        return;
    end
    try
        maxclust = str2num(get(handles.MaxClusters,'String'));
    catch
        return;
    end 
    
    resultsclust = cell((1+ maxclust)-minclust,1);
    
    hw = waitbar(0,'Clustering data set ');
    for i = 1:(1+maxclust)-minclust
        resultsclust{i} = clusterdata(matrix2classify,minclust + (i-1),handles.distancemeasure);
        if length(resultsclust{i}) == 1
            msgbox(['Clustering failed at ' num2str(resultsclust{i}) ' clusters. Try changing the range of clusters or rerunning.']);
            return;
        end
        waitbar(i/(1+maxclust-minclust));
    end
    close(hw);
    
    silhmean = zeros((1+ maxclust)-minclust,1);
    
    for i = 1:1+(maxclust)-minclust
        
        if diagnosticstate
            figure();
            [silh h] = silhouette(matrix2classify,resultsclust{i});
            pcaplot(matrix2classify,resultsclust{i});
        else
            silh = silhouette(matrix2classify,resultsclust{i});
        end
            
        silhmean(i) = mean(silh);
    end
    
    [x I] = max(silhmean);
    handles.classification = resultsclust{I(1)};
    
    set(handles.KClasses,'String',num2str(minclust + (I(1)-1)));
end
end
    end


set(handles.AcceptButton,'Enable','on');

guidata(gcbo,handles);

% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = 0;
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes on button press in AcceptButton.
function AcceptButton_Callback(hObject, eventdata, handles)
% hObject    handle to AcceptButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = handles.classification;
guidata(hObject,handles);
uiresume(handles.figure1);


function CepstralEdit_Callback(hObject, eventdata, handles)
% hObject    handle to CepstralEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CepstralEdit as text
%        str2double(get(hObject,'String')) returns contents of CepstralEdit as a double


% --- Executes during object creation, after setting all properties.
function CepstralEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CepstralEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function MinClusters_Callback(hObject, eventdata, handles)
% hObject    handle to MinClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinClusters as text
%        str2double(get(hObject,'String')) returns contents of MinClusters as a double


% --- Executes during object creation, after setting all properties.
function MinClusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function MaxClusters_Callback(hObject, eventdata, handles)
% hObject    handle to MaxClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxClusters as text
%        str2double(get(hObject,'String')) returns contents of MaxClusters as a double


% --- Executes during object creation, after setting all properties.
function MaxClusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in TransformPopupMenu.
function TransformPopupMenu_Callback(hObject, eventdata, handles)
% hObject    handle to TransformPopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns TransformPopupMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TransformPopupMenu

% value = get(hObject,'Value');
% 
% if value == 2
%     handles.matrix2classify(:,1) = log(handles.matrix2classify(:,1)); 
% else
%     handles.matrix2classify(:,1) = exp(handles.matrix2classify(:,1));
% end

% --- Executes during object creation, after setting all properties.
function TransformPopupMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TransformPopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function KClasses_Callback(hObject, eventdata, handles)
% hObject    handle to KClasses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of KClasses as text
%        str2double(get(hObject,'String')) returns contents of KClasses as a double


% --- Executes during object creation, after setting all properties.
function KClasses_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KClasses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in CepstralPopupMenu.
function CepstralPopupMenu_Callback(hObject, eventdata, handles)
% hObject    handle to CepstralPopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns CepstralPopupMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CepstralPopupMenu


% --- Executes during object creation, after setting all properties.
function CepstralPopupMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CepstralPopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --- Executes on button press in RangeSpecify.
function RangeSpecify_Callback(hObject, eventdata, handles)
% hObject    handle to RangeSpecify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RangeSpecify

state = get(hObject,'Value');

if state
    set(handles.KClasses,'Enable','off');
    set(handles.text3,'Visible','on');
    set(handles.MinClusters,'Visible','on');
    set(handles.text4,'Visible','on');
    set(handles.MaxClusters,'Visible','on');
else
    set(handles.KClasses,'Enable','on');
    set(handles.text3,'Visible','off');
    set(handles.MinClusters,'Visible','off');
    set(handles.text4,'Visible','off');
    set(handles.MaxClusters,'Visible','off');
end


% --- Executes on selection change in DistancePopupMenu.
function DistancePopupMenu_Callback(hObject, eventdata, handles)
% hObject    handle to DistancePopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns DistancePopupMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DistancePopupMenu

contents = get(hObject,'String');

dstring = contents{get(hObject,'Value')};

switch dstring
    case 'Squared Euclidean'
        handles.distancemeasure = 'sqEuclidean'
    case 'City Block (L1)'
        handles.distancemeasure = 'cityblock'
    case 'Cosine'
        handles.distancemeasure = 'cosine'
    case 'Correlation'
        handles.distancemeasure = 'correlation'
end

guidata(gcbo,handles);

% --- Executes during object creation, after setting all properties.
function DistancePopupMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DistancePopupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

handles.distancemeasure = 'sqEuclidean';

guidata(gcbo,handles);

% --- Executes on button press in DiagnosticCheckbox.
function DiagnosticCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to DiagnosticCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DiagnosticCheckbox

function pcaplot(matrix2classify,classification)
% Plots the first two principal components uses different shapes selected
% randomly and different colors selected randomly to show classes.

nclusters = length(unique(classification));

if not(isempty(find(classification == 0)))
    nclusters = nclusters -1;
end

colors = rand(nclusters,3);
mark = mod(1:nclusters,13)+1;

[coef, score] = princomp(matrix2classify);

figure();

for i = 1:length(classification)
    h = line(score(i,1),score(i,2));
    markers = set(h,'Marker');
    if classification(i) == 0 % For segments which are not classified
        set(h,'Marker',markers{4});   
        set(h,'Color',[0 0 0]);
    else
        set(h,'Marker',markers{mark(classification(i))});   
        set(h,'Color',colors(classification(i),:));
    end
end

xlabel('PCA 1');
ylabel('PCA 2');

function classifymatrix = madnormalize(classifymatrix, cols2normalize)
  for i = cols2normalize
    classifymatrix(:,i) = (classifymatrix(:,i) - median(classifymatrix(:,i))) / mad(classifymatrix(:,i));
  end
  
  
% This code is modified from http://phys.columbia.edu/~aylin/clustering

function c=cluster_aylin(p1,p2);
%enter the p1,p2 found from hneighbors.m and get the struct array c. 
%c(i).c will contain the indices of the objects that is in cluster i.
%c(1).c will be the cluster that has the maximum number of objects.
nn=size(p1,2);
p=std(p1);g=p(p2);p=p';
z=zeros(nn,1);zz=zeros(nn,1);
n=zeros(1,nn);
for i=1:nn
    n(i)=max(max(p1(:,p2(:,i))));
end
gg=[1:nn; max(g); n]';n=find(n<1);
[q1,q1]=sort(gg(n,2)');
g=gg(n(q1),1);
t=cputime;
j=1;m=1
	for i=1:length(n)
        ii=g(i);b=p2(:,ii);a=b;aa=a(find(gg(a,2)<1));
        if length([0; unique(z(aa))])<=2;a1=1;a2=0;
            while a1~=a2;a1=length(aa);
                a=unique(p2(:,a));a=a(find(gg(a,2)<1));
                a=unique(a(find(gg(a,2)<=mean(gg(aa,2))+std(gg(aa,2)))));
                a=a(find(ismember(a,aa)==0));
                if ~isempty(a);aa=[aa;a];
                    jj=aa(find(z(aa)));
                    if ~isempty(jj);;
                         u=unique(z(jj));
                         if length(u)==1;
                            zz(aa(find(~z(aa))))=m;z(aa)=u;m=m+1;
                         end
                         break;
                    end
                end;a2=length(aa);
            end;a=aa;
            jj=a(find(z(a)));
            if isempty(jj);
                z(a)=j;j=j+1;
                zz(a)=m;zz(b)=m-.1;zz(ii)=m-.2;m=m+1;
            end
        end
    end
    u=unique(z)
    u=u(find(u));v=length(u);vv=floor(1/v);if vv;vv=' is';else vv='s      are';end
    fprintf([int2str(v) ' cluster' vv ' found in ' num2str(cputime-t) 'sec\n\n'])
	q0=[];for i=1:v;qq=find(z==u(i));q0=[q0 length(qq)];end;
    [q1,q2]=sort(q0);%q2=q2(find(100*q1/nn>1));v=length(q2);
	c=[];for i=1:v;qq=find(z==u(q2(i)));[j,j]=sort(zz(qq));c(v-i+1).c=qq(j);end;
    if isempty(c);fprintf('\tNo cluster was found. \n \tscale the data (step size must be 1)\n');end 
    
    
    function [p1,p2]=hneighbors(e);
	% this function finds the neighbors of each object in 'e' within a unit hypercube
	% and returns the sorted object distances to q.q1 and their identities to q.q2 , 'e' is a 
	% matrix where i th row and j th column represents the j th component of the i th object. 
	s='find(';i='abs(e(:,%d)-e(i,%d))<1&';
    % assuming that the data is given scaled and the characteristic step size is 1, variable s 
	% keeps a script to find the objects that lie within a unit hypercube around the i th object.
    for j=1:size(e,2);
        s=[s sprintf(i,j,j)];
	end;s([end end+1])=');';
     % runs the script s for each of the objects and stores the sorted distances 
     % from the i th object in q(i).q1 and their indentities in q(i).q2
    nn=size(e,1);m=ceil(nn^(1/4));p1=ones(m,nn);p2=kron(ones(m,1),1:nn);
    for i=1:nn;
        j=eval(s);
        [q,qq]=sort(sqrt(sum((e(j,:)-kron(ones(length(j),1),e(i,:)))'.^2)));q=q(find(q<1));mm=length(q);
        qq=j(qq(1:mm))';mn=min([m mm]);
        p1(1:mn,i)=q(1:mn)';
        p2(1:mn,i)=qq(1:mn)';
    end
    
    
function classification = cluster_hierarchical(matrix2classify)
    [p1,p2] = hneighbors(matrix2classify);
    c = cluster_aylin(p1,p2);
    
    nclusters = size(c,2); 
    
    total_classified = 0
    
    nsegments = size(matrix2classify,1);
            
    classification = zeros(nsegments,1);
    
    for i = 1:nclusters
        for j = 1:length(c(i).c)
            classification(c(i).c(j)) = i;
        end
    end
;


% --- Executes on button press in HierarchicalRadioButton.
function HierarchicalRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to HierarchicalRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HierarchicalRadioButton

selected = get(hObject,'Value');

if selected
    handles.cluster_method = 'hierarchical';
end

guidata(gcbo,handles);


% --- Executes on button press in KmeansRadioButton.
function KmeansRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to KmeansRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of KmeansRadioButton

selected = get(hObject,'Value');

if selected
    handles.cluster_method = 'kmeans';
end

guidata(gcbo,handles);




function SegLenWeight_Callback(hObject, eventdata, handles)
% hObject    handle to SegLenWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SegLenWeight as text
%        str2double(get(hObject,'String')) returns contents of SegLenWeight as a double


% --- Executes during object creation, after setting all properties.
function SegLenWeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SegLenWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


