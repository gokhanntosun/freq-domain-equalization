clear;close all;
%definitions---
h=[0.74 -0.514 0.37 0.216 0.062];%given channel model
N=1000;%block length
Ps=[1 3 5 10];%cyclic prefix count
frame_cnt=2000;%frame count
fer_lim=10000;%frame error limit
SNR_dB=0:2:20;%snr values in dbs
SNR=1./(2*(10.^(SNR_dB./10)));%snr values in bit energy
BER_FDE=zeros(length(Ps),length(SNR_dB));%define BER
FFT_Pt=N;%n-point fft
H=fft(h,FFT_Pt);%n-point fft of channel impulse response
C=1./H;%freq domain zero-forcing equalizer
for p=1:length(Ps)
    P=Ps(p);%pilot count
    for s=1:length(SNR_dB)
        var=SNR(s);%corresponding variance for snr
        fr=0;fer=0;%initialize for sim
        %begin monte carlo sim---
        while fr<frame_cnt && fer<fer_lim            
            %transmitter---
            x=randi([0 1],[1 N]);x(x==0)=-1;%generate block
            x_p=[x(:,end-P+1:end) x];%add x cyclic prefix
            %channel---
            noise=normrnd(0,sqrt(var),[1 N+P+length(h)-1]);%noise samples
            y=conv(h,x_p)+noise;%impose channel conditions
            %receiver---
            y=y(:,P+1:end-length(h)+1);%discard cyclic prefix
            Y=fft(y,FFT_Pt);%dft
            O=Y.*C;%equalize in freq domain
            o=ifft(O,FFT_Pt);%inverse dft
            o(o<0)=-1;o(o>=0)=1;%detection
            diff=nnz(x-o);%compute block error
            fer=fer+diff;%accumulate block errors
            fr=fr+1;%increment frame count
        end
        BER_FDE(p,s)=fer/(N*fr);%compute ber
        display([SNR_dB(s) fr fer]);%display sim parameters
    end
end
%plot results---
load 'CE_MMSE'
ss=get(0,'ScreenSize');
figure;
semilogy(SNR_dB,BER_FDE(1,:),'*-');hold on;
semilogy(SNR_dB,BER_FDE(2,:),'*-');hold on;
semilogy(SNR_dB,BER_FDE(3,:),'*-');hold on;
semilogy(SNR_dB,BER_FDE(4,:),'*-');hold on;
semilogy(SNR_dB,BER_MMSE(1:length(SNR_dB)),'k-^');hold on;
legend('P=1','P=3','P=5','P=10','MMSE','Location','SouthWest');
xlabel('SNR(dB)');ylabel('BER');
title('FDE | Zero-Forcing Eq.');
grid on;axis square;
set(gca,'FontSize',14);