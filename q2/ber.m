%%%***---------------------------求ber曲线----------------------------***%%%

tic
clc
clear
close all
SNRdB=[2:13];  %信噪比(dB)的范围
for j=1:length(SNRdB)
    snr(j)=10^(SNRdB(j)/10);
    err_rate(j)=qpsk_with_ber(snr(j)) ;
end
for i=1:length(SNRdB)
    snr=10^(SNRdB(i)/10);
    theo_pb(i)=(1/2)*erfc(sqrt(snr));
end
figure
semilogy(SNRdB,err_rate,'*');hold on
semilogy(SNRdB,theo_pb);grid;
title('QPSK的ber曲线');
xlabel('Eb/n0(dB)');ylabel('P_e');
legend('仿真误比特率','过莱斯信道+高斯噪声+迫零均衡误比特率');

toc