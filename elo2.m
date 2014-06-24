function [ t, d_d, dane, bpsk_mod, szum, nosna_p, odfiltr, ipHatn, simBER ] = elo2( Eb_N0_dB )
x = inputdlg('Wprowadz sekwencje zer i jedynek, oddzielajac je spacjami:',...
             'Dane', [1 100]);
dane = str2num(x{1,:});                             %  sekwencje bitow
%dane = round(rand(1,1200));  
SNR = Eb_N0_dB;                                     %10*log10(Eb_N0_dB)+10*log10(1/1); % multiple Eb/N0 values; log10(bitrate/Bandwidth)
NRZ = 2*dane-1;                                     % Kodowanie NRZ
T = 1;                                              % Czas trwania bitu (okres)
fc = 3/T;                                           % Czestotlowosc nosna
E = 5;                                              % Energia sygna?l
t = linspace(0,length(dane),length(dane)*100);      % probki, czasy
N = length(t);                                      % Liczba probek
Lpnb = N/length(dane);                              % Liczba probek na bit
dod_dane = repmat(dane',1,Lpnb);                    % powtarzamy bity Lpnb razy
dod_dane2 = repmat(NRZ',1,Lpnb);
d_d = dod_dane';                                    % Transpozycja wierszy i kolumn
d_d2 = dod_dane2';
d_d = d_d(:)';                                      % Konwersja do jednego wiersza (100x1, 100x0, itd..)
d_d2 = d_d2(:)';     
nosna = sin(pi*fc*t);                               % Fala nosna (nosna= sqrt(2*(E/T))*cos(2*pi*fc*t);)
bpsk_mod = d_d2.*nosna;                             % Modulowany przebieg
szum = awgn(bpsk_mod,SNR);                          % Szum

%nosna_p = sqrt(2*(E/T))*cos(2*pi*fc*t);
nosna_p = sin(pi*fc*t); 
bpsk_mod_p = d_d2.*nosna_p;
szum_p = awgn(bpsk_mod_p,SNR); 

kwadrat2 = szum_p.^2;                               % podnosze sygnal do kwadratu 

filtr=fir1(100,[29 31]/(length(dane)*100/2));       % tworzenie filtru
pofiltr=filter(filtr,1,kwadrat2);                   % filtrowanie sygnalu 
x = pofiltr;
%x=szum_p(1,1:length(szum_p)-127); %zmniejszam sygna? czasu
                                     %o odpowiednie op??nienie, wyliczone
                                     %jako stopie? wielomianu przy filtrze
                                     %(512) podzielony przez 2.
%x=pofiltr(1,128:length(szum_p)); %przesuwam sygna? zakodowany
                                     %o okre?lon? liczb? pr?bek, ma to na
                                     %celu wyeliminowanie b??du
                                    

nosna2=x./2;                                        %obnizam amplitude aby jej wartosc zgadzala sie z sygnalem poczatkowym
%nowanosna=0.5:1:length(nosna2);                    %bardziej 'g?sta' no?na, b?dzie potrzebna 
                                                    %przy interpolacji
%nosna2=interp1(1:length(nosna2),nosna2,nowanosna);  % interpoluj? no?n?, zwi?kszam liczb? pr?bek 
                                                    %do 2000 i nast?pnie wycinam 1000,
                                                    %dzi?ki temu zmniejszam cz?stotliwo??
                                                    %sygna?u o po?ow? 
%nosna2(1)=nosna2(2);     %problem z pierwszym elementem wektora, zapisuj? 
                     %do niego warto?? 2 elementu
%nosna2=nosna2(1,1:length(szum_p));
                   
% subplot(2,1,2);
% plot(t,nosna2)
% title('Odzyskana nosna ')





%{
[A, B] = butter(1,0.99,'low');
yout = filter(A,B,szum');
%}
%{
b=fidemo.fi_fir_coefficients;
zi = zeros(length(b),1);

%codegen_fi doc model example
%Initialize the output signal yout and the final conditions zf
Ty = numerictype(1,12,8);
yout = fi(zeros(size(szum)),'numerictype',Ty);
zf = zi;

% FIR filter code
for k=1:length(szum);
  % Update the states: z = [x(k);z(1:end-1)]  
  zf(:) = [szum(k);zf(1:end-1)];
  % Form the output: y(k) = b*z
  yout(k) = b*zf;
end

odfiltr=yout.*nosna;
%}
%nosna= sqrt(2*(E/T))*sin(2*pi*fc*t);

odfiltr=szum.*nosna;

for i=1:length(odfiltr)
    if odfiltr(i)>0
        bpsk_demod(i)=1;
    else
        bpsk_demod(i)=0;
    end
end

%%
%%%%%%%%%%%%%% Plotujemy sobie %%%%%%%%%%%%%%%%%%%%%

% figure(1)
% subplot(6,1,1);                 % Dana sekwencja
% plot(t,d_d); 
% axis([0 length(dane) -0.5 1.5])
% 
% % subplot(4,2,3);                 % Po NRZ
% % plot(t,d_d2); 
% % axis([0 length(dane) -1.5 1.5])
% 
% % subplot(4,2,5);                 % Nosna
% % plot(t,nosna); 
% % axis([0 length(dane) -1.5 1.5]) 
% 
% subplot(6,1,2);                 % Zmodulowany
% plot(t,bpsk_mod); 
% axis([0 length(dane) -3 3])
% 
% subplot(6,1,3);                 % Szumy
% plot(t,szum); 
% axis([0 length(dane) -2.5 2.5])
% 
% % subplot(6,1,4);                 % Po filtracji fir
% % plot(t,szum,'.'); 
% % axis([0 length(dane) -3 3])
% 
% subplot(6,1,5);                 % Demodulator PRACUJE
% plot(t,odfiltr);
% axis([0 length(dane) -3 3])
% 
% subplot(6,1,6);                 % Nasz sekwencja po demodulatorze
% plot(t,bpsk_demod); 
% axis([0 length(dane) -0.5 1.5])

%%
% Bledy - zwraca ich ilosc

y = bpsk_demod;                  
% odbiornik - dekodowanie hard decision
ipHat = real(y)>0;
z=zeros(1,length(dane));
j=z;

for a=0:length(dane)-1
    for b=1:15
        if ipHat(100*a+6*b)==0
            z(a+1)=z(a+1)+1;
        else j(a+1)=j(a+1)+1;
        end     
    end
    
end

for a=0:length(dane)-1
    for b=1:15
        if z(a+1)>j(a+1)
            ipHatn(15*a+b)=0;
        else
            ipHatn(15*a+b)=1;
        end
    end
end
z
j

nErr = size(find([d_d- ipHat]),2);
simBER = nErr/length(d_d)


end
