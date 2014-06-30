function varargout = gui(varargin)
%GUI M-file for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('Property','Value',...) creates a new GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUI('CALLBACK') and GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help gui
% Last Modified by GUIDE v2.5 29-Jun-2014 15:26:39
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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
% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
% Choose default command line output for gui
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;
% --- Executes on slider movement.
function ebn0_Callback(hObject, eventdata, handles)
% hObject    handle to ebn0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ebn0=get(hObject,'Value');
set(handles.s2,'String',round(get(hObject,'Value')));
guidata(hObject, handles); 
% Hints: get(hObject,'String') returns contents of ebn0 as text
%        str2double(get(hObject,'String')) returns contents of ebn0 as a double
% --- Executes during object creation, after setting all properties.
function ebn0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebn0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes on button press in b2.
function b2_Callback(hObject, eventdata, handles)
% hObject    handle to b2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.s3,'String',simBER);
% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%[ t, d_d, dane, bpsk_mod, szum, nosna_p, odfiltr, wnew, simBER ] = elo2( handles.ebn0 );
a = get(handles.dzikie,'String');


b=31;
register1=[1 1 1 1 1];
code1=zeros(1,b);
for i=1:b
   temp = mod(register1(2)+register1(5),2);
   code1(i) = 2*register1(5)-1;
   for j=5:-1:2
      register1(j)=register1(j-1);
   end
   register1(1) = temp;
end
m_sequence_1=code1;
m_sequence_1=m_sequence_1';
m_sequence_1 =m_sequence_1>0;


l=length(m_sequence_1);
czas=0:1:(l-1);
axes(handles.axes2)
stem(czas,m_sequence_1)
m_sequence_1=num2str(m_sequence_1);
l = length(a);
b = a -'0';
c=m_sequence_1-'0';
c=c';
k=1;
g=length(m_sequence_1);
for i=1:l
for j=1:g
    
    spread(1,k)=xor(b(1,i),c(1,j));
    k=k+1;
end
end
guidata(hObject, handles);


axes(handles.axes3)
l=length(spread);
z=0:1:(l-1);
stem(z,spread);
% % Modulating the hopped signal
% dsss_sig=[];
% z=0:1:(l-1);   
% fc=30;
% c1=cos(pi*fc*z);
% c2=cos(pi*fc*z+pi);
% for k=1:120
%     if spread(1,k)==0
%         dsss_sig=[dsss_sig c1];
%     else
%         dsss_sig=[dsss_sig c2];
%     end
% end
% axes(handles.axes6)
% plot(dsss_sig);
a = get(handles.dzikie,'String');
a=a-'0';
NRZ = 2*a-1;                                        % Kodowanie NRZ
T = 1;                                              % Czas trwania bitu (okres)
fc = 3/T;                                           % Czestotlowosc nosna
                                         % Czestotlowosc nosna
t = linspace(0,length(NRZ),length(NRZ)*100);      % probki, czasy
N = length(t);                                      % Liczba probek
Lpnb = N/length(NRZ);                              % Liczba probek na bit
                    % powtarzamy bity Lpnb razy
dod_dane2 = repmat(NRZ',1,Lpnb);
                                   % Transpozycja wierszy i kolumn
d_d2 = dod_dane2';
                                    % Konwersja do jednego wiersza (100x1, 100x0, itd..)
d_d2 = d_d2(:)'; 
spread2 = 2.*spread-1;

spread2 = repmat(NRZ',1,Lpnb);
                                   % Transpozycja wierszy i kolumn
spread2 = spread2';
                                    % Konwersja do jednego wiersza (100x1, 100x0, itd..)
spread2 = spread2(:)'; 
nosna = sin(pi*fc*t);                               % Fala nosna 
%bpsk_mod = d_d2.*nosna;
bpsk_mod = spread2.*nosna;
axes(handles.axes6)
plot(t,bpsk_mod);

szum = awgn(bpsk_mod,handles.ebn0);    
axes(handles.axes5);
plot(t,szum); 
axis([0 length(NRZ) -2.5 2.5])


for i=1:length(spread)
    odfiltr(i)=szum(i).*nosna(i);
end

odfiltr=szum.*nosna;

for i=1:length(odfiltr)
    if odfiltr(i)>0
        bpsk_demod(i)=1;
    else
        bpsk_demod(i)=0;
    end
end


% Bledy - zwraca ich ilosc
y = bpsk_demod;                  
w = real(y)>0;
z=zeros(1,length(NRZ));
j=z;
bit=z;
for q=0:length(NRZ)-1
    for b=1:15
        if w(100*q+6*b)==0
            z(q+1)=z(q+1)+1;
        else j(q+1)=j(q+1)+1;
        end     
    end
    
end
for q=0:length(NRZ)-1
    for b=1:15
        if z(q+1)>j(q+1)
            wnew(15*q+b)=0;
            bit(q+1)=0;
        else
            wnew(15*q+b)=1;
            bit(q+1)=1;
        end
    end
end
% nErr = size(find([d_d- w]),2);
% simBER = nErr/length(d_d)
nErr = size(find([a- bit]),2);
simBER = nErr/length(a);
% axes(handles.axes1);
% plot(t,d_d); 
% axis([0 length(dane) -0.5 1.5])
% axes(handles.axes2);
% plot(t,bpsk_mod,'.'); 
% axis([0 length(dane) -1.5 1.5])
% axes(handles.axes3);
% plot(t,szum); 
% axis([0 length(dane) -2.5 2.5])
% axes(handles.axes6);
% plot(t, nosna_p)
% axis([0 length(dane) -2.5 2.5])
% axes(handles.axes5);
% plot(t, odfiltr)
% axis([0 length(dane) -1.5 1.5])
axes(handles.axes4);
%t = linspace(0,length(NRZ),length(NRZ)*15);  

% axis([0 length(dane) -1.5 1.5])

a = get(handles.dzikie,'String');
l=length(a);
z=0:1:(l-1);
%bit = a - '0';
axes(handles.axes4);
stem(z,bit); 
axis([0 length(NRZ)-1 0 1])
set(handles.s3,'String',simBER);


% stem(t,wnew); 
% axis([0 length(NRZ) -0.5 1.5])
% set(handles.s3,'String',simBER);

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BBER()
function dzikie_Callback(hObject, eventdata, handles)
% hObject    handle to dzikie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of dzikie as text
%        str2double(get(hObject,'String')) returns contents of dzikie as a double
a = get(handles.dzikie,'String');
guidata(hObject, handles);
l = length(a);
t=0:1:(l-1);
b = a-'0';
 
axes(handles.axes1)
stem(t,b)
axis([0 l-1 0 1])
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dzikie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dzikie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
