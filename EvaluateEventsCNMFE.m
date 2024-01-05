function varargout = EvaluateEventsCNMFE(varargin)
% EVALUATEEVENTSCNMFE MATLAB code for EvaluateEventsCNMFE.fig
%      EVALUATEEVENTSCNMFE, by itself, creates a new EVALUATEEVENTSCNMFE or raises the existing
%      singleton*.
%
%      H = EVALUATEEVENTSCNMFE returns the handle to a new EVALUATEEVENTSCNMFE or the handle to
%      the existing singleton*.
%
%      EVALUATEEVENTSCNMFE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVALUATEEVENTSCNMFE.M with the given input arguments.
%
%      EVALUATEEVENTSCNMFE('Property','Value',...) creates a new EVALUATEEVENTSCNMFE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EvaluateEventsCNMFE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EvaluateEventsCNMFE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EvaluateEventsCNMFE

% Last Modified by GUIDE v2.5 15-May-2020 13:40:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EvaluateEventsCNMFE_OpeningFcn, ...
                   'gui_OutputFcn',  @EvaluateEventsCNMFE_OutputFcn, ...
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


% --- Executes just before EvaluateEventsCNMFE is made visible.
function EvaluateEventsCNMFE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EvaluateEventsCNMFE (see VARARGIN)

% Choose default command line output for EvaluateEventsCNMFE
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EvaluateEventsCNMFE wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EvaluateEventsCNMFE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in acceptbutton.
function acceptbutton_Callback(hObject, eventdata, handles)
% hObject    handle to acceptbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.currlocs==length(handles.locs)
    msgbox('All events evaluated');
else
    handles.currlocs=handles.currlocs+1;
    clear handles.EventPlot;
    plot(handles.t,handles.currtrace,'Parent', handles.EventPlot);
    hold(handles.EventPlot,'on')
    plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(handles.currlocs)),'r*','Parent',handles.EventPlot);
    hold(handles.EventPlot,'off')
    currframes=imread(handles.fullnamevid,handles.locs(handles.currlocs));
    axes(handles.NeuronPlot)
    imagesc(currframes,[handles.lim1 handles.lim2]);
    hold(handles.NeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
    hold(handles.NeuronPlot,'off')
    set(handles.currenteventtext, 'String', handles.currlocs);
    
    axes(handles.NeuronZoomPlot);
    imagesc(currframes, [handles.lim1 handles.lim2]);
    colormap(gray);
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    th=0:pi/50:2*pi;
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    
    hold(handles.NeuronZoomPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
    hold(handles.NeuronZoomPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
end
end    
guidata(hObject,handles);


% --- Executes on button press in rejectbutton.
function rejectbutton_Callback(hObject, eventdata, handles)
% hObject    handle to rejectbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.currlocs==1
    firstloc=1;
else
    firstloc=handles.locs(handles.currlocs-1);
end
if handles.currlocs==length(handles.locs)
    lastloc=length(handles.traces.c(handles.currtracenum,:));
else
    lastloc=handles.locs(handles.currlocs+1);
end
while (handles.traces.c(handles.currtracenum,firstloc)>min(handles.traces.c(handles.currtracenum,firstloc:handles.locs(handles.currlocs))))&&firstloc<handles.locs(handles.currlocs)
    firstloc=firstloc+1;
end
while (handles.traces.c(handles.currtracenum,lastloc)>min(handles.traces.c(handles.currtracenum,handles.locs(handles.currlocs):lastloc)))&&lastloc>handles.locs(handles.currlocs)
    lastloc=lastloc-1;
end

for i=firstloc:lastloc
    handles.traces.c(handles.currtracenum,i)=0;
    handles.traces.s(handles.currtracenum,i)=0;
end

handles.currtrace=handles.traces.c(handles.currtracenum,:);
[handles.pks,handles.locs]=findpeaks(handles.traces.s(handles.currtracenum,:));
%clear handles.EventPlot;
plot(handles.t,handles.currtrace,'Parent', handles.EventPlot);
if handles.currlocs>length(handles.locs)
    msgbox('All events evaluated');
else
    hold(handles.EventPlot,'on')
    plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(handles.currlocs)),'r*','Parent',handles.EventPlot);
    hold(handles.EventPlot,'off')
    currframes=imread(handles.fullnamevid,handles.locs(handles.currlocs));
    axes(handles.NeuronPlot)
imagesc(currframes,[handles.lim1 handles.lim2]);
    hold(handles.NeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
    hold(handles.NeuronPlot,'off')
    set(handles.currenteventtext, 'String', handles.currlocs);
    
    axes(handles.NeuronZoomPlot);
    imagesc(currframes, [handles.lim1 handles.lim2]);
    colormap(gray);
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    th=0:pi/50:2*pi;
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
 
    hold(handles.NeuronZoomPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
    hold(handles.NeuronZoomPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
end
end
guidata(hObject,handles);


% --- Executes on button press in loaddatabutton.
function loaddatabutton_Callback(hObject, eventdata, handles)
% hObject    handle to loaddatabutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, foldername] = uigetfile('*.mat', 'Select a mat file');
fullname = fullfile(foldername, filename);
temp = load(fullname);
fn = fieldnames(temp);
handles.traces = temp.(fn{1});

handles.currtrace = handles.traces.c(1,:);
handles.currtracenum = 1;
handles.t=[1:length(handles.currtrace)]./20;

set(handles.currenttracenumbertext, 'String', handles.currtracenum);
set(handles.currentfiletext, 'String', fullname);

[handles.pks,handles.locs]=findpeaks(handles.traces.s(1,:));
handles.currlocs=1;
%axes(handles.EventPlot);
plot(handles.t,handles.currtrace,'Parent',handles.EventPlot);
hold(handles.EventPlot,'on')
plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(1)),'r*','Parent',handles.EventPlot);
hold(handles.EventPlot,'off')
set(handles.currenteventtext, 'String', handles.currlocs);
set(handles.totaleventtext, 'String', length(handles.locs));

xlabel('Time (sec)');
ylabel('Signal');
xticks([0 300 600 900 1200 1500 1800]);
xlim([0 1800]);

guidata(hObject,handles);

% --- Executes on button press in nexteventbutton.
function nexteventbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nexteventbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.currlocs==length(handles.locs)
    msgbox('You are viewing the last event');
else
    handles.currlocs=handles.currlocs+1;
    handles.currtrace=handles.traces.c(handles.currtracenum,:);
    [handles.pks,handles.locs]=findpeaks(handles.traces.s(handles.currtracenum,:));
    clear handles.EventPlot;
    plot(handles.t,handles.currtrace,'Parent', handles.EventPlot);
    hold(handles.EventPlot,'on')
    plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(handles.currlocs)),'r*','Parent',handles.EventPlot);
    hold(handles.EventPlot,'off')
    currframes=imread(handles.fullnamevid,handles.locs(handles.currlocs));
    axes(handles.NeuronPlot)
    imagesc(currframes,[handles.lim1 handles.lim2]);
    hold(handles.NeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
    hold(handles.NeuronPlot,'off')
    set(handles.currenteventtext, 'String', handles.currlocs);
    
    axes(handles.NeuronZoomPlot);
    imagesc(currframes, [handles.lim1 handles.lim2]);
    colormap(gray);
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
        th=0:pi/50:2*pi;
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);

    hold(handles.NeuronZoomPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
    hold(handles.NeuronZoomPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
    end
end
guidata(hObject,handles);

% --- Executes on button press in previouseventbutton.
function previouseventbutton_Callback(hObject, eventdata, handles)
% hObject    handle to previouseventbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.currlocs==1
    msgbox('You are viewing the last event');
else
    handles.currlocs=handles.currlocs-1;
    handles.currtrace=handles.traces.c(handles.currtracenum,:);
    [handles.pks,handles.locs]=findpeaks(handles.traces.s(handles.currtracenum,:));
    clear handles.EventPlot;
    plot(handles.t,handles.currtrace,'Parent', handles.EventPlot);
    hold(handles.EventPlot,'on')
    plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(handles.currlocs)),'r*','Parent',handles.EventPlot);
    hold(handles.EventPlot,'off') 
    currframes=imread(handles.fullnamevid,handles.locs(handles.currlocs));
    axes(handles.NeuronPlot)
    imagesc(currframes,[handles.lim1 handles.lim2]);
    hold(handles.NeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
    hold(handles.NeuronPlot,'off')
    set(handles.currenteventtext, 'String', handles.currlocs);
    
    axes(handles.NeuronZoomPlot);
    imagesc(currframes, [handles.lim1 handles.lim2]);
    colormap(gray);
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    th=0:pi/50:2*pi;
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    
    hold(handles.NeuronZoomPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
    hold(handles.NeuronZoomPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
    end
end
guidata(hObject,handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function gototext_Callback(hObject, eventdata, handles)
% hObject    handle to gototext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gototext as text
%        str2double(get(hObject,'String')) returns contents of gototext as a double
gotonum = str2double(get(hObject,'String'));

if handles.currlocs>length(handles.locs)||handles.currlocs<1
    msgbox('Outside of bounds');
else
    handles.currlocs=gotonum;
    plot(handles.t,handles.currtrace,'Parent', handles.EventPlot);
    hold(handles.EventPlot,'on')
    plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(handles.currlocs)),'r*','Parent',handles.EventPlot);
    hold(handles.EventPlot,'off')
    currframes=imread(handles.fullnamevid,handles.locs(handles.currlocs));
    axes(handles.NeuronPlot)
    imagesc(currframes,[handles.lim1 handles.lim2]);
    hold(handles.NeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
    hold(handles.NeuronPlot,'off')
    set(handles.currenteventtext, 'String', handles.currlocs);
    
    axes(handles.NeuronZoomPlot);
    imagesc(currframes, [handles.lim1 handles.lim2]);
    colormap(gray);
    th=0:pi/50:2*pi;
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    
    hold(handles.NeuronZoomPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
    hold(handles.NeuronZoomPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
end
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function gototext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gototext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in alldonebutton.
function alldonebutton_Callback(hObject, eventdata, handles)
% hObject    handle to alldonebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','modneuron',handles.traces)

% --- Executes on button press in nextneuronbutton.
function nextneuronbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nextneuronbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.currtracenum == length(handles.traces.c(:,1)))
    msgbox('You are viewing the last trace');
else 
    c=[randi(100)/100 randi(100)/100 randi(100)/100];
    hold(handles.AllNeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'MarkerEdgeColor',c,'Parent', handles.AllNeuronPlot);
    hold(handles.AllNeuronPlot,'off')
    handles.currtracenum = handles.currtracenum+1;
    handles.currtrace = handles.traces.c(handles.currtracenum,:);
    handles.currlocs=1;
    [handles.pks,handles.locs]=findpeaks(handles.traces.s(handles.currtracenum,:));
    plot(handles.t,handles.currtrace,'Parent', handles.EventPlot);
    hold(handles.EventPlot,'on')
    plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(handles.currlocs)),'r*','Parent',handles.EventPlot);
    hold(handles.EventPlot,'off')
    xlabel('Time (sec)');
    ylabel('Signal');
    xticks([0 300 600 900 1200 1500 1800]);
    xlim([0 1800]);
    set(handles.currenttracenumbertext, 'String', handles.currtracenum);
    currframes=imread(handles.fullnamevid,handles.locs(handles.currlocs)); 
    axes(handles.NeuronPlot)
    imagesc(currframes,[handles.lim1 handles.lim2]);
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    th=0:pi/50:2*pi;
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    
    hold(handles.NeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
    hold(handles.NeuronPlot,'off')
    set(handles.currenteventtext, 'String', handles.currlocs);
    set(handles.totaleventtext, 'String', length(handles.locs));
    guidata(hObject, handles);
    
    axes(handles.NeuronZoomPlot);
    imagesc(currframes, [handles.lim1 handles.lim2]);
    colormap(gray);
    hold(handles.NeuronZoomPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
    hold(handles.NeuronZoomPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
    end

    n=1;
    while ismember(n,handles.locs)
    n=n+1;
    end
    cf=imread(handles.fullnamevid,n);
    sm=double(cf);
    track=1;
    addamt=round(length(handles.traces.c(1,:))/10);
    n=n+addamt;
        while n<length(handles.traces.c(1,:))
                cf=imread(handles.fullnamevid,n);
                if ~ismember(n,handles.locs)
                    sm=double(sm)+double(cf);
                    track=track+1;
                    n=n+addamt;
                else
                    n=n+addamt;
                end
        end
        mncf=double(sm/track);

        axes(handles.NeuronNegPlot);
    imagesc(mncf, [handles.lim1 handles.lim2]);
    th=0:pi/50:2*pi;
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    hold(handles.NeuronNegPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronNegPlot);
    hold(handles.NeuronNegPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
    end
    hold(handles.AllNeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.AllNeuronPlot);
    hold(handles.AllNeuronPlot,'off')
end


% --- Executes on button press in previousneuronbutton.
function previousneuronbutton_Callback(hObject, eventdata, handles)
% hObject    handle to previousneuronbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.currtracenum == 1)
    msgbox('You are viewing the first trace');
else 
    c=[randi(100)/100 randi(100)/100 randi(100)/100];
    hold(handles.AllNeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'MarkerEdgeColor',c,'Parent', handles.AllNeuronPlot);
    hold(handles.AllNeuronPlot,'off')
    handles.currtracenum = handles.currtracenum-1;
    handles.currtrace = handles.traces.c(handles.currtracenum,:);
    handles.currlocs=1;
    plot(handles.t,handles.currtrace,'Parent', handles.EventPlot);
    hold(handles.EventPlot,'on')
    plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(handles.currlocs)),'r*','Parent',handles.EventPlot);
    hold(handles.EventPlot,'off')
    xlabel('Time (sec)');
    ylabel('Signal');
    xticks([0 300 600 900 1200 1500 1800]);
    xlim([0 1800]);
    set(handles.currenttracenumbertext, 'String', handles.currtracenum);
    [handles.pks,handles.locs]=findpeaks(handles.traces.s(handles.currtracenum,:));
    currframes=imread(handles.fullnamevid,handles.locs(handles.currlocs));
    axes(handles.NeuronPlot)
    imagesc(currframes,[handles.lim1 handles.lim2]);
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    th=0:pi/50:2*pi;
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    
    handles.currn.yext=handles.currn.yext;
    hold(handles.NeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
    hold(handles.NeuronPlot,'off')
    set(handles.currenteventtext, 'String', handles.currlocs);
    set(handles.totaleventtext, 'String', length(handles.locs));
    guidata(hObject, handles);
    
    axes(handles.NeuronZoomPlot);
    imagesc(currframes, [handles.lim1 handles.lim2]);
    colormap(gray);
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    hold(handles.NeuronZoomPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
    hold(handles.NeuronZoomPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
    end

    n=1;
    while ismember(n,handles.locs)
     n=n+1;
    end
    cf=imread(handles.fullnamevid,n);
    sm=double(cf);
    track=1;
    addamt=round(length(handles.traces.c(1,:))/10);
    n=n+addamt;
        while n<length(handles.traces.c(1,:))
                cf=imread(handles.fullnamevid,n);
                if ~ismember(n,handles.locs)
                    sm=double(sm)+double(cf);
                    track=track+1;
                    n=n+addamt;
                else
                    n=n+addamt;
                end
        end
        mncf=double(sm/track);

        axes(handles.NeuronNegPlot);
    imagesc(mncf, [handles.lim1 handles.lim2]);
    th=0:pi/50:2*pi;
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    hold(handles.NeuronNegPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronNegPlot);
    hold(handles.NeuronNegPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
    end

    hold(handles.AllNeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.AllNeuronPlot);
    hold(handles.AllNeuronPlot,'off')
end

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function gotoneurontext_Callback(hObject, eventdata, handles)
% hObject    handle to gotoneurontext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gotoneurontext as text
%        str2double(get(hObject,'String')) returns contents of gotoneurontext as a double
gotonum = str2double(get(hObject,'String'));
if gotonum > length(handles.traces.c(:,1)) || gotonum <= 0
    msgbox('Record is out of bounds');
else 
    c=[randi(100)/100 randi(100)/100 randi(100)/100];
    hold(handles.AllNeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'MarkerEdgeColor',c,'Parent', handles.AllNeuronPlot);
    hold(handles.AllNeuronPlot,'off')
    
    handles.currtracenum = gotonum;
    handles.currtrace = handles.traces.c(handles.currtracenum,:);
    handles.currlocs=1;
    plot(handles.t,handles.currtrace,'Parent', handles.EventPlot);
    hold(handles.EventPlot,'on')
    plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(handles.currlocs)),'r*','Parent',handles.EventPlot);
    hold(handles.EventPlot,'off')
    xlabel('Time (sec)');
    ylabel('Signal');
    xticks([0 300 600 900 1200 1500 1800]);
    xlim([0 1800]);
    set(handles.currenttracenumbertext, 'String', handles.currtracenum);
    [handles.pks,handles.locs]=findpeaks(handles.traces.s(handles.currtracenum,:));
    currframes=imread(handles.fullnamevid,handles.locs(handles.currlocs));
    axes(handles.NeuronPlot)
    imagesc(currframes,[handles.lim1 handles.lim2]);
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    th=0:pi/50:2*pi;
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    
    hold(handles.NeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
    hold(handles.NeuronPlot,'off')
    set(handles.currenteventtext, 'String', handles.currlocs);
    set(handles.totaleventtext, 'String', length(handles.locs));
    guidata(hObject, handles);
    
    axes(handles.NeuronZoomPlot);
    imagesc(currframes, [handles.lim1 handles.lim2]);
    colormap(gray);
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    th=0:pi/50:2*pi;
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    
    hold(handles.NeuronZoomPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
    hold(handles.NeuronZoomPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
    end

    n=1;
    while ismember(n,handles.locs)
        n=n+1;
    end
    cf=imread(handles.fullnamevid,n);
    sm=double(cf);
    track=1;
    addamt=round(length(handles.traces.c(1,:))/10);
    n=n+addamt;
        while n<length(handles.traces.c(1,:))
                cf=imread(handles.fullnamevid,n);
                if ~ismember(n,handles.locs)
                    sm=double(sm)+double(cf);
                    track=track+1;
                    n=n+addamt;
                else
                    n=n+addamt;
                end
        end
        mncf=double(sm/track);

        axes(handles.NeuronNegPlot);
    imagesc(mncf, [handles.lim1 handles.lim2]);
    th=0:pi/50:2*pi;
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    hold(handles.NeuronNegPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronNegPlot);
    hold(handles.NeuronNegPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
    end
    hold(handles.AllNeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.AllNeuronPlot);
    hold(handles.AllNeuronPlot,'off')
end



% --- Executes during object creation, after setting all properties.
function gotoneurontext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gotoneurontext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadvideobutton.
function loadvideobutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadvideobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.traces.rejectnums=[];
[filename, foldername] = uigetfile('*.tiff', 'Select a tiff file');
handles.fullnamevid = fullfile(foldername, filename);
%for i=1:length(info)
%    handles.rawframes(:,:,i)=imread(filename,i);
%end
currframes=imread(handles.fullnamevid,handles.locs(1));
xr=[1 length(currframes(1,:))];
yr=[1 length(currframes(:,1))];
axes(handles.NeuronPlot)
imagesc(currframes,[handles.lim1 handles.lim2]);
colormap('gray');
th=0:pi/50:2*pi;
%if ~exist('handles.currn.yext','var')
%   handles.circles=1;
%end
handles.currcenter=handles.traces.center(handles.currtracenum,:);
handles.currn.yext=3*sin(th)+handles.currcenter(1);
handles.currn.xext=3*cos(th)+handles.currcenter(2);
hold(handles.NeuronPlot,'on')
scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
hold(handles.NeuronPlot,'off')
set(handles.currenttracenumbertext, 'String', handles.currtracenum);

axes(handles.AllNeuronPlot);
imagesc(handles.traces.Cn,'Parent', handles.AllNeuronPlot);
colormap('gray');
hold(handles.AllNeuronPlot,'on')
for i=1:length(handles.traces.c(:,1))
    handles.currcenter=handles.traces.center(i,:);
    handles.currneu.yext=3*sin(th)+handles.currcenter(1);
    handles.currneu.xext=3*cos(th)+handles.currcenter(2);
     c=[randi(100)/100 randi(100)/100 randi(100)/100];
    scatter(handles.currneu.xext,handles.currneu.yext,1,'MarkerEdgeColor',c,'Parent', handles.AllNeuronPlot);
end
scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.AllNeuronPlot);
hold(handles.AllNeuronPlot,'off')

axes(handles.NeuronZoomPlot);
handles.currcenter=handles.traces.center(handles.currtracenum,:);
imagesc(currframes, [handles.lim1 handles.lim2]);
colormap(gray);
hold(handles.NeuronZoomPlot,'on')
scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
hold(handles.NeuronZoomPlot,'off')
x0=handles.currcenter(2); y0=handles.currcenter(1);
if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
end


n=1;
while ismember(n,handles.locs)
    n=n+1;
end
cf=imread(handles.fullnamevid,n);
sm=double(cf);
track=1;
addamt=round(length(handles.traces.c(1,:))/10);
n=n+addamt;
        while n<length(handles.traces.c(1,:))
                cf=imread(handles.fullnamevid,n);
                if ~ismember(n,handles.locs)
                    sm=double(sm)+double(cf);
                    track=track+1;
                    n=n+addamt;
                else
                    n=n+addamt;
                end
        end
        mncf=double(sm/track);

        axes(handles.NeuronNegPlot);
imagesc(mncf, [handles.lim1 handles.lim2]);
th=0:pi/50:2*pi;
hold(handles.NeuronNegPlot,'on')
scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronNegPlot);
hold(handles.NeuronNegPlot,'off')
x0=handles.currcenter(2); y0=handles.currcenter(1);
if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
end

axes(handles.NeuronPlot);

guidata(hObject, handles);


% --- Executes on button press in rejectneuronbutton.
function rejectneuronbutton_Callback(hObject, eventdata, handles)
% hObject    handle to rejectneuronbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.traces.rejectnums=[handles.traces.rejectnums handles.currtracenum];
hold(handles.AllNeuronPlot,'on')
c=[randi(100)/100 randi(100)/100 randi(100)/100];
scatter(handles.currn.xext,handles.currn.yext,1,c,'Parent', handles.AllNeuronPlot);
hold(handles.AllNeuronPlot,'off')

handles.traces.c(handles.currtracenum,:)=[];
handles.traces.s(handles.currtracenum,:)=[];
handles.traces.center(handles.currtracenum,:)=[];
handles.currtracenum=handles.currtracenum-1;
    handles.currtrace = handles.traces.c(handles.currtracenum,:);
    handles.currlocs=1;
    plot(handles.t,handles.currtrace,'Parent', handles.EventPlot);
    hold(handles.EventPlot,'on')
    plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(handles.currlocs)),'r*','Parent',handles.EventPlot);
    hold(handles.EventPlot,'off')
    xlabel('Time (sec)');
    ylabel('Signal');
    xticks([0 300 600 900 1200 1500 1800]);
    xlim([0 1800]);
    set(handles.currenttracenumbertext, 'String', handles.currtracenum);
    [handles.pks,handles.locs]=findpeaks(handles.traces.s(handles.currtracenum,:));
    currframes=imread(handles.fullnamevid,handles.locs(handles.currlocs));
    axes(handles.NeuronPlot)
    imagesc(currframes,[handles.lim1 handles.lim2]);
    th=0:pi/50:2*pi;
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    hold(handles.NeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
    hold(handles.NeuronPlot,'off')
    set(handles.currenteventtext, 'String', handles.currlocs);
    set(handles.totaleventtext, 'String', length(handles.locs));
    
    axes(handles.NeuronZoomPlot);
    imagesc(currframes, [handles.lim1 handles.lim2]);
    colormap(gray);
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    hold(handles.NeuronZoomPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
    hold(handles.NeuronZoomPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
    end

    n=1;
    while ismember(n,handles.locs)
        n=n+1;
    end
    cf=imread(handles.fullnamevid,n);
    sm=double(cf);
    track=1;
    addamt=round(length(handles.traces.c(1,:))/10);
    n=n+addamt;
        while n<length(handles.traces.c(1,:))
                cf=imread(handles.fullnamevid,n);
                if ~ismember(n,handles.locs)
                    sm=double(sm)+double(cf);
                    track=track+1;
                    n=n+addamt;
                else
                    n=n+addamt;
                end
        end
        mncf=double(sm/track);

        axes(handles.NeuronNegPlot);
    imagesc(mncf, [handles.lim1 handles.lim2]);
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    hold(handles.NeuronNegPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronNegPlot);
    hold(handles.NeuronNegPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
    if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
    end
    
    hold(handles.AllNeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.AllNeuronPlot);
    hold(handles.AllNeuronPlot,'off')

    
    guidata(hObject, handles);



function Lim1Text_Callback(hObject, eventdata, handles)
% hObject    handle to Lim1Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lim1Text as text
%        str2double(get(hObject,'String')) returns contents of Lim1Text as a double
handles.lim1 = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Lim1Text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lim1Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lim2Text_Callback(hObject, eventdata, handles)
% hObject    handle to Lim2Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lim2Text as text
%        str2double(get(hObject,'String')) returns contents of Lim2Text as a double
handles.lim2= str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Lim2Text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lim2Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in statcheckbox.
function statcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to statcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of statcheckbox
handles.circles=get(hObject,'Value');
guidata(hObject, handles);



function nteststext_Callback(hObject, eventdata, handles)
% hObject    handle to nteststext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nteststext as text
%        str2double(get(hObject,'String')) returns contents of nteststext as a double
handles.testnum= str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nteststext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nteststext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in seemaxbutton.
function seemaxbutton_Callback(hObject, eventdata, handles)
% hObject    handle to seemaxbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.pks,handles.locs]=findpeaks(handles.traces.s(handles.currtracenum,:));

lng=1:length(handles.locs);
idx=lng(handles.locs==handles.locs(handles.pks==max(handles.pks)));

    handles.currlocs=idx;
    plot(handles.t,handles.currtrace,'Parent', handles.EventPlot);
    hold(handles.EventPlot,'on')
    plot(handles.locs(handles.currlocs)/20,handles.currtrace(handles.locs(handles.currlocs)),'r*','Parent',handles.EventPlot);
    hold(handles.EventPlot,'off')
    currframes=imread(handles.fullnamevid,handles.locs(handles.currlocs));
    axes(handles.NeuronPlot)
    imagesc(currframes,[handles.lim1 handles.lim2]);
    hold(handles.NeuronPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronPlot);
    hold(handles.NeuronPlot,'off')
    set(handles.currenteventtext, 'String', handles.currlocs);
    
    axes(handles.NeuronZoomPlot);
    imagesc(currframes, [handles.lim1 handles.lim2]);
    colormap(gray);
    th=0:pi/50:2*pi;
    handles.currcenter=handles.traces.center(handles.currtracenum,:);
    handles.currn.yext=3*sin(th)+handles.currcenter(1);
    handles.currn.xext=3*cos(th)+handles.currcenter(2);
    
    hold(handles.NeuronZoomPlot,'on')
    scatter(handles.currn.xext,handles.currn.yext,1,'y','Parent', handles.NeuronZoomPlot);
    hold(handles.NeuronZoomPlot,'off')
    x0=handles.currcenter(2); y0=handles.currcenter(1);
if ~isnan(x0)
        xlim(x0+[-12, 12]*2);
        ylim(y0+[-12, 12]*2);
end


guidata(hObject, handles);
