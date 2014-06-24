function [  ] = BBER(  )
    figure();
    SNR = -4:12;
    EbN0=SNR(1):0.001:SNR(length(SNR));
    error = (1/2)*erfc(sqrt(10.^(EbN0/10)));
    sError  =  importdata('matlab.mat');
    semilogy(SNR,sError,'*r',...
                   EbN0 ,error,'green');
    title('BER')
    xlabel('SNR [dB]');
    ylabel('BER')
    legend('BPSK Symulacja','BPSK Teoria')
   
end