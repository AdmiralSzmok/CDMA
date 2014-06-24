function [ t, d_d, dane, bpsk_mod, szum, nosna_p, odfiltr, wnew, simBER ] = elo2( Eb_N0_dB )
x = inputdlg('Wprowadz sekwencje zer i jedynek, oddzielajac je spacjami:',...
             'Dane', [1 100]);
dane = str2num(x{1,:});                             %  sekwencje bitow
%dane = round(rand(1,1000000));  
SNR = Eb_N0_dB;                                     %10*log10(Eb_N0_dB)+10*log10(1/1); % multiple Eb/N0 values; log10(bitrate/Bandwidth)
NRZ = 2*dane-1;                                     % Kodowanie NRZ
T = 1;                                              % Czas trwania bitu (okres)
fc = 3/T;                                           % Czestotlowosc nosna
t = linspace(0,length(dane),length(dane)*100);      % probki, czasy
N = length(t);                                      % Liczba probek
Lpnb = N/length(dane);                              % Liczba probek na bit
dod_dane = repmat(dane',1,Lpnb);                    % powtarzamy bity Lpnb razy
dod_dane2 = repmat(NRZ',1,Lpnb);
d_d = dod_dane';                                    % Transpozycja wierszy i kolumn
d_d2 = dod_dane2';
d_d = d_d(:)';                                      % Konwersja do jednego wiersza (100x1, 100x0, itd..)
d_d2 = d_d2(:)';     
nosna = sin(pi*fc*t);                               % Fala nosna 
bpsk_mod = d_d2.*nosna;                             % Modulowany przebieg
szum = awgn(bpsk_mod,SNR);                          % Szum

nosna_p = sin(pi*fc*t); 
bpsk_mod_p = d_d2.*nosna_p;
szum_p = awgn(bpsk_mod_p,SNR); 

kwadrat2 = szum_p.^2;                               % podnosze sygnal do kwadratu 

filtr=fir1(30,0.99);                                % tworzenie filtru
pofiltr=filter(filtr,1,kwadrat2);                   % filtrowanie sygnalu 
x = pofiltr;
                                   

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
w = real(y)>0;
z=zeros(1,length(dane));
j=z;
bit=z;

for a=0:length(dane)-1
    for b=1:15
        if w(100*a+6*b)==0
            z(a+1)=z(a+1)+1;
        else j(a+1)=j(a+1)+1;
        end     
    end
    
end

for a=0:length(dane)-1
    for b=1:15
        if z(a+1)>j(a+1)
            wnew(15*a+b)=0;
            bit(a+1)=0;
        else
            wnew(15*a+b)=1;
            bit(a+1)=1;
        end
    end
end

% nErr = size(find([d_d- w]),2);
% simBER = nErr/length(d_d)
nErr = size(find([dane- bit]),2);
simBER = nErr/length(dane);


end
