function [simBER] = nowa(snr)
a = round(rand(1,10000));

b=40;
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

l=length(spread);
z=0:1:(l-1);
stem(z,spread);


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
nosna = sin(pi*fc*t);                               % Fala nosna 
bpsk_mod = d_d2.*nosna;    


szum = awgn(bpsk_mod,snr);    

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

nErr = size(find([a- bit]),2);
simBER = nErr/length(a)
