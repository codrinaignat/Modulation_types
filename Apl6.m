function varargout = Apl6(varargin)
% Last Modified by GUIDE v2.5 04-May-2023 15:58:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Apl6_OpeningFcn, ...
                   'gui_OutputFcn',  @Apl6_OutputFcn, ...
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

function Apl6_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
ah = axes('unit', 'normalized', 'position', [0 0 1 1]); 
bg = imread('sky.jpg'); imagesc(bg);
set(ah,'handlevisibility','off','visible','off')
uistack(ah, 'bottom');

guidata(hObject, handles);

function varargout = Apl6_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function pushbutton1_Callback(hObject, eventdata, handles)
AppInfo();

function popupmenu2_Callback(hObject, eventdata, handles)
A_carrier = str2num(get(handles.edit8, 'string'));
F_carrier = str2num(get(handles.edit9, 'string'));
A_modulator = str2num(get(handles.edit10, 'string'));
F_modulator = str2num(get(handles.edit11, 'string'));
M_index = str2num(get(handles.edit12, 'string'));
Bits = get(handles.edit13, 'string');
Bit_stream = str2num(Bits);
durata = str2num(get(handles.edit14, 'string'));
Square_duty_cycle = str2num(get(handles.edit19, 'string'));
Fs_PCM = str2num(get(handles.edit17, 'string'));
Numar_biti_cuantizare = str2num(get(handles.edit18, 'string'));

Fs = 40000;
Ts = 1/Fs;
N = 5000;
t = 0: Ts :durata;

Mod_type = get(handles.popupmenu2, 'value')
switch Mod_type
    case 1 %Amplitude Modulation
        Semnal_purtator=A_carrier*sin(2*pi*F_carrier*t);
        Semnal_modulator=sin(2*pi*F_modulator*t);
        ModulatAM = (A_carrier*(1+M_index*Semnal_modulator)).*sin(2*pi*F_carrier*t);
        
        axes(handles.axes6);
        plot(t,Semnal_purtator);
        title('Carrier signal - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes8);
        plot(t, Semnal_modulator);
        title('Modulator signal - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes9);
        plot(t,ModulatAM);
        title('AM modulated signal - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes7);
        axa1 = linspace(-Fs/2, Fs/2, length(ModulatAM)); 
        AmpMod = abs(fftshift(fft(ModulatAM))); 
        plot(axa1, AmpMod);
        title('AM modulated signal - frequency domain');
        zoom on
        
    case 2 %Frequency Modulation
        Semnal_purtator=A_carrier*sin(2*pi*F_carrier*t);
        Semnal_modulator=sin(2*pi*F_modulator*t);
        ModulatFM = A_carrier*sin(2*pi*F_carrier*t + M_index.*Semnal_modulator);
        
        axes(handles.axes6);
        plot(t,Semnal_purtator);
        title('Carrier signal - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes8);
        plot(t, Semnal_modulator);
        title('Modulator signal - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes9);
        plot(t,ModulatFM);
        title('FM modulated signal - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes7);
        axa1 = linspace(-Fs/2, Fs/2, length(ModulatFM)); 
        FreqMod = abs(fftshift(fft(ModulatFM))); 
        plot(axa1, FreqMod);
        title('FM modulated signal - frequency domain');
        zoom on     
        
    case 3 % Phase Modulation
        Semnal_purtator=A_carrier*sin(2*pi*F_carrier*t);
        Semnal_modulator = cos(2*pi*F_carrier*t + M_index*(2*pi*F_modulator*t)); 
        ModulatPM=A_carrier*sin(2*pi*F_carrier*t + cos(2*pi*F_carrier*t + M_index*(2*pi*F_modulator*t)));
        
        axes(handles.axes6);
        plot(t,Semnal_purtator);
        title('Carrier signal 1 - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes8);
        plot(t, Semnal_modulator);
        title('Carrier signal 2- time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes9);
        plot(t,ModulatPM);
        title('PM modulated signal - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes7);
        axa1 = linspace(-Fs/2, Fs/2, length(ModulatPM)); 
        PhaseMod = abs(fftshift(fft(ModulatPM))); 
        plot(axa1, PhaseMod);
        title('PM modulated signal - frequency domain');
        zoom on  
    case 4 % PAM
        Semnal_modulator=A_modulator*sin(2*pi*F_modulator*t);
        Semnal_purtator=A_carrier*square(2*pi*F_carrier*t);
        ModulatPAM=Semnal_purtator.*Semnal_modulator;
        
        axes(handles.axes6);
        plot(t,Semnal_purtator);
        axis([0 durata 0 (A_carrier+0.5)]);
        title('Carrier signal - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes8);
        plot(t, Semnal_modulator);
        title('Modulator signal - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes9);
        plot(t,ModulatPAM);
        title('PAM modulated signal - time domain');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes7);
        axa1 = linspace(-Fs/2, Fs/2, length(ModulatPAM)); 
        PAMMod = abs(fftshift(fft(ModulatPAM))); 
        plot(axa1, PAMMod);
        title('PAM modulated signal - frequency domain');
        zoom on
    
    case 5 %PWM
        Semnal_dinte_fierastrau = A_carrier.*sawtooth(2*pi*F_carrier*t); 
        Semnal_modulator = A_modulator.*sin(2*pi*F_modulator*t);
        for i=1:length(Semnal_dinte_fierastrau)
            if (Semnal_modulator(i)>=Semnal_dinte_fierastrau(i))
                ModulatPWM(i) = 1;
            else
                ModulatPWM(i) = 0;
            end
        end
        
        axes(handles.axes6);
        plot(t,Semnal_dinte_fierastrau);
        xlabel('t');
        ylabel('Amplitude [V]');
        title('Comparator signal');
        zoom on
        
        axes(handles.axes8);
        plot(t, Semnal_modulator);
        xlabel('t');
        ylabel('Amplitude [V]');
        title('Informational signal');
        zoom on
        
        axes(handles.axes9);
        plot(t,ModulatPWM);
        xlabel('t');
        ylabel('Amplitude [V]');
        title('PWM resultant signal');
        axis([0 durata 0 1.1]);
        zoom on
        
        axes(handles.axes7);
        axa1 = linspace(-Fs/2, Fs/2, length(ModulatPWM)); 
        PWMMod = abs(fftshift(fft(ModulatPWM))); 
        plot(axa1, PWMMod);
        title('PWM signal spectrum');
        zoom on
        
    case 6 %PPM
        t_ppm = 1;
        n = [0:1/Fs:t_ppm];
        n = n(1:end - 1);
        Period = Fs/F_carrier;
        T_on = Period/Square_duty_cycle;

        Semnal_purtator1 = square(2*pi*F_carrier*n,Square_duty_cycle);
        Semnal_purtator1(find(Semnal_purtator1 < 0)) = 0;
        Semnal_modulator = A_carrier*sin(2*pi*F_modulator*n);
        Semnal_purtator2=A_carrier.*sawtooth(2*pi*F_carrier*n);%Carrier sawtooth

        ModulatPPM = zeros(1,length(Semnal_purtator1));

        Compare1 = find(Semnal_purtator2 > Semnal_modulator);
        Compare2 = diff(Compare1);
        Compare3 = find(Compare2 ~= 1);
        Temporary(1) = Compare1(1);
        Temporary(2:length(Compare3)+1) = Compare1(Compare3 + 1);

        for i = 1:length(Temporary)
            ModulatPPM(Temporary(i) : Temporary(i) + T_on - 1) = 1;
        end
        
        axes(handles.axes6);
        plot(n,Semnal_purtator1,'LineWidth',2);
        axis([0 durata 0 1.1]);
        xlabel('t');
        ylabel('Amplitude [V]');
        title('Carrier signal 1');
        zoom on
        
        axes(handles.axes8);
        plot(n,Semnal_purtator2,'LineWidth',2);
        xlabel('t');
        ylabel('Amplitude [V]');
        title('Carrier signal 2');
        zoom on
        
        axes(handles.axes9);
        plot(n,Semnal_modulator,'LineWidth',2);
        xlabel('t');
        ylabel('Amplitude [V]');
        title('Modulator signal - time domain');
        zoom on
        
        axes(handles.axes7);
        plot(n,ModulatPPM,'LineWidth',2);
        axis([0 durata 0 1.1]);
        xlabel('t');
        ylabel('Amplitude [V]');
        title('PPM modulated signal - time domain');
        zoom on
        
    case 7 %ASK
        numar_biti = length(Bit_stream);
        t_local = 0:0.01:numar_biti;
        x = 1:1:(numar_biti+1)*100;
        for i = 1:numar_biti
            for j = i:0.1:i+1
                Bit_array(x(i*100:(i+1)*100)) = Bit_stream(i);
            end
        end
        

        Bit_array = Bit_array(100:end);
        Semnal_purtator = sin(2*pi*F_carrier*t_local);
        ModulatASK = Bit_array.*Semnal_purtator;
        
        axes(handles.axes6);
        plot(t_local,Bit_array);
        xlabel('t');
        ylabel('Amplitude [V]');
        title('Bit Array');
        axis([0 numar_biti -0.1 1.1])
        zoom on
        
        axes(handles.axes8);
        plot(t_local,Semnal_purtator);
        xlabel('t');
        ylabel('Amplitude [V]');
        title('Carrier signal - time domain ');
        zoom on
        
        axes(handles.axes9);
        plot(t_local, ModulatASK);
        xlabel('t');
        ylabel('Amplitude [V]');
        title('ASK modulated signal - time domain');
        zoom on
        
        blank = imread('white.jpg');
        axes(handles.axes7);
        imshow(blank);

    case 8 %FSK
        Bit_period = 0.000001;
        bit=[];
        for n=1:1:length(Bit_stream)
            if Bit_stream(n)==1;
                Level = ones(1,100);
            else Bit_stream(n)==0;
                Level = zeros(1,100);
            end
            bit=[bit Level];

        end
        t_local1=Bit_period/100:Bit_period/100:100*length(Bit_stream)*(Bit_period/100);
        
        Bit_rate = 1/Bit_period;                                                         % bit rate
        F_carrier = Bit_rate*8;                           % carrier frequency for information as 1
        F_modulator=Bit_rate*2;                           % carrier frequency for information as 0
        t_local2=Bit_period/99:Bit_period/99:Bit_period;
        ss=length(t_local2);
        ModulatFSK=[];
        for (i=1:1:length(Bit_stream))
            if (Bit_stream(i)==1)
                Output = A_carrier*cos(2*pi*F_carrier*t_local2);
            else
                Output = A_carrier*cos(2*pi*F_modulator*t_local2);
            end
            ModulatFSK=[ModulatFSK Output];
        end
        t_local3=Bit_period/99:Bit_period/99:Bit_period*length(Bit_stream);
        
        axes(handles.axes6);
        plot(t_local2, A_carrier*cos(2*pi*F_carrier*t_local2));
        title('Signal corresponding to logical value 1');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes8);
        plot(t_local2, A_carrier*cos(2*pi*F_modulator*t_local2));
        title('Signal corresponding to logical value 0');
        xlabel('t');
        ylabel('Amplitude [V]');
        zoom on
        
        axes(handles.axes9);
        plot(t_local1,bit,'lineWidth',2.5);grid on;
        axis([ 0 Bit_period*length(Bit_stream) -.5 1.5]);
        ylabel('Amplitude [V]');
        xlabel('t');
        title('Binary information');
        zoom on
        
        axes(handles.axes7);
        plot(t_local3,ModulatFSK);
        ylabel('Amplitude [V]');
        xlabel('t');
        title('FSK modulated signal - time domain');
        zoom on
        
    case 9 %PSK-BPSK
        numar_biti = length(Bit_stream);
        Bit_period = 0.0001;   

        bit = []; 
        for n = 1:1:numar_biti    
            if Bit_stream(n) == 1;   
               Level = ones(1,100);
            else Bit_stream(n) == 0;
                Level = zeros(1,100);
            end
             bit = [bit Level];
        end

        t_local1=Bit_period/100:Bit_period/100:100*numar_biti*(Bit_period/100);   
        
        axes(handles.axes6);
        plot(t_local1,bit,'lineWidth',2.5);
        grid on;
        axis([0 Bit_period*numar_biti -0.5 1.5]);
        ylabel('Amplitude [V]');
        xlabel('t');
        title('Binary information');
        zoom on
        
        Bit_rate = 1/Bit_period;   
        F_carrier = Bit_rate*10;  
        Phasec1 = 0;     
        Phasec2 = pi;    
        t_local2 = Bit_period/100:Bit_period/100:Bit_period;                  

        ModulatBPSK = [];
        for (i = 1:1:numar_biti)
            if (Bit_stream(i)==1)
                Output = A_carrier*cos(2*pi*F_carrier*t_local2+Phasec1);   
            else
                Output = A_carrier*cos(2*pi*F_carrier*t_local2+Phasec2);   
            end
            ModulatBPSK=[ModulatBPSK Output];
        end

        t_local3=Bit_period/100:Bit_period/100:Bit_period*numar_biti;  
        
        axes(handles.axes8);
        plot(t_local3,ModulatBPSK);
        ylabel('Amplitude [V]');
        xlabel('t');
        title('BPSK modulated signal - time domain');
        zoom on
        
        blank = imread('white.jpg');
        axes(handles.axes9);
        imshow(blank);
        
        blank = imread('white.jpg');
        axes(handles.axes7);
        imshow(blank);

    case 10 %PCM
%         Fs_PCM = 40;
        nT = 0: 1/Fs_PCM :durata;
        
        Semnal_modulator = A_modulator*sin(2*pi*F_modulator*nT);
%         Numar_biti_cuantizare = 8; % 256 quantization levels
        Nivele_cuantizare = 2^Numar_biti_cuantizare;
        A_minim = -A_modulator  ;
        A_maxim = A_modulator;  
        factor_scalare = (A_maxim-A_minim)/Nivele_cuantizare;
        Semnal_prelucrat = Semnal_modulator / factor_scalare ;
        Semnal_prelucrat = floor(Semnal_prelucrat);
        Semnal_prelucrat = Semnal_prelucrat * factor_scalare;
        
        axes(handles.axes6);
        plot(nT, Semnal_modulator);
        xlabel('Time (sec)');
        ylabel('Amplitude');
        title('Sine wave - time domain');
        
        axes(handles.axes8);
        stem(nT, Semnal_modulator);
        xlabel('Time (sec)');
        ylabel('Amplitude');
        title('Sampled signal');
        
        axes(handles.axes9);
        plot(nT, Semnal_prelucrat);
        stairs(nT, Semnal_prelucrat);
        xlabel('Time (sec)');
        ylabel('Amplitude');
        title('Quantized signal');
        
        axes(handles.axes7);
        stem(nT,Semnal_prelucrat);
        xlabel('Time (sec)');
        ylabel('Amplitude');
        title('Sampled and quantized signal');
end

function popupmenu2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function axes3_CreateFcn(hObject, eventdata, handles)

function pushbutton9_Callback(hObject, eventdata, handles)
Help();

function edit2_Callback(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)

function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)

function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)

function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)

function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit8_Callback(hObject, eventdata, handles)

function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit9_Callback(hObject, eventdata, handles)

function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit10_Callback(hObject, eventdata, handles)

function edit10_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit11_Callback(hObject, eventdata, handles)

function edit11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit12_Callback(hObject, eventdata, handles)

function edit12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit13_Callback(hObject, eventdata, handles)

function edit13_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit14_Callback(hObject, eventdata, handles)

function edit14_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit15_Callback(hObject, eventdata, handles)

function edit15_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit16_Callback(hObject, eventdata, handles)

function edit16_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)

function edit17_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit18_Callback(hObject, ~, ~)

function edit18_CreateFcn(hObject, eventdata, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton10_Callback(hObject, ~, ~)
theory = imread('Modulation_types.jpg');
figure(1)
imshow(theory);

function edit19_Callback(hObject, eventdata, handles)

function edit19_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
