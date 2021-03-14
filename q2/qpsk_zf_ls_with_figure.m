%%%%%----------------------      题目三  qpsk3.m     ------------------%%%%%
%%%  没有乘法器那边的调制，是映射之后ifft求时域波形
clear
close all
tic
fs=240;%符号率：240MHz
dt=1/fs;
n=1024;
a=11;
m=n*2*(a-1);
symbol=60;%符号率：60MHz
SNR=30;%snr=30

%% 信号产生
signal1=randi([0 1],1,m);
signal=(signal1*2-1);
%% qpsk调制-->iq映射
ich=[];
qch=[];
for j=1:m/2
    ich(j)=signal(2*j-1);
    qch(j)=signal(2*j);
end

%% qpsk调制-->产生复信号
kmod=1./sqrt(2);
ich1=ich.*kmod;
qch1=qch.*kmod;
x=ich1+qch1.*1i;

%% 插导频
frank=frank(n);
x_re=reshape(x,[n,a-1]);
for j=1:a
    if(j==1)
        jiadaopin(:,j)=frank';
    else
        jiadaopin(:,j)=x_re(:,j-1);
    end
end
ch2=ifft(jiadaopin);%%转换为时域信号
%% 加CP
for j=1:a
        ch3(:,j)=[ch2(end-127:end,j);ch2(:,j)];
end

%% P/S
ch4=reshape(ch3,1,[]);

%% 上采样和脉冲成型滤波

% 脉冲形成：升余弦滤波器
r = 0.25; %滚降系数
delay = 8; %延时
IPOINT = 4; %内插倍数
h_rcos = rcosdesign(r, delay, IPOINT, 'sqrt'); %生成时域响应h
len_fir = length(h_rcos);
len_ch4 = length(ch4);

signal_txr = upfirdn(ch4, h_rcos, 4 ); 
signal_tx=signal_txr(2*delay-1:end-2*delay+1);
eyediagram(fft(signal_tx),2);
t=0:dt/4:(length(abs(signal_tx))-1)*dt/4;    %时间向量
figure;plot(t,abs(signal_tx));title('A点波形');xlabel('t/us');ylabel('幅值/V');
scatterplot(signal_tx);title('A点星座图');
signal_tx_psd=abs(fftshift(fft(signal_tx))).*conj(abs(fftshift(fft(signal_tx))));
w=(-length(signal_tx_psd)/2:1:length(signal_tx_psd)/2-1)*fs/4/length(signal_tx_psd); 
figure;plot(w,signal_tx_psd);title('A点功率谱');xlabel('频率/mHz');ylabel('幅值/V');
%% 过莱斯信道

%%%%%%%%莱斯信道
%%%%%%多径时延：[0 0.5us 1.5us]
%%%%%%多径功率：[0 -3dB -9 dB]
n_delay1=fix(0.5*symbol);
n_delay2=fix(1.5*symbol);
signal_1=signal_tx;
signal_2=10^(-3/20)*[zeros(1,n_delay1),signal_1(1:end-n_delay1)];
signal_3=10^(-9/20)*[zeros(1,n_delay2),signal_1(1:end-n_delay2)];
signal_rx_raly=signal_1+signal_2+signal_3;
% signal_rx_raly=signal_tx;
%% 加高斯噪声
signal_rx=awgn(signal_rx_raly,SNR - 10*log10(IPOINT),'measured');
% signal_rx=signal_tx;
%% 匹配滤波
% signal_irx=signal_itx;
% signal_qrx=signal_qtx;
signal_fira = upfirdn(signal_rx, h_rcos, 1); %对信号进行匹配滤波

signal_fir=signal_fira(2*delay+1:end-2*delay);

scatterplot(fft(signal_fir));title('B点星座图');
t=0:dt/4:(length(abs(signal_fir))-1)*dt/4;    %时间向量
figure;plot(t,abs(signal_fir));title('B点波形');xlabel('t/us');ylabel('幅值/V');
eyediagram(fft(signal_fir),2);
signal_fir_psd=abs(fftshift(fft(signal_tx))).*conj(abs(fftshift(fft(signal_tx))));
w=(-length(signal_fir_psd)/2:1:length(signal_fir_psd)/2-1)*fs/4/length(signal_fir_psd); 
figure;plot(w,signal_fir_psd);title('B点功率谱');xlabel('频率/mHz');ylabel('幅值/V');
%% 下采样

receice_s=downsample(signal_fir, 4);

%% S/P
receice_p=reshape(receice_s,[],a);

%% 去cp

receice_p_ncp=receice_p(129:end,:);

%% fft
fft_data=fft(receice_p_ncp);

%% 信道估计与均衡-->信道估计（频域ls估计）
frank_bn=fft_data(:,1);
%%%频域信道估计值：
H=frank_bn./frank';
%%%时域信道估计值
h=ifft(H);
figure,impz(h)%信道时域响应
figure,freqz(h)%信道频域响应
%% 信道估计与均衡-->频域迫零均衡
data_bal_f=fft_data./H;


%% P/S+去导频
data_bal_p=reshape(data_bal_f(:,2:end),1,[]);
scatterplot(data_bal_p);title('C点星座图');
eyediagram(data_bal_p,2);
t=0:dt/4:(length(abs(ifft(data_bal_p)))-1)*dt;    %时间向量
figure;plot(abs(ifft(data_bal_p)));title('C点波形');xlabel('t/us');ylabel('幅值/V');
data_bal_p_psd=abs(fftshift(data_bal_p)).*conj(abs(fftshift(data_bal_p)));
w=(-length(data_bal_p_psd)/2:1:length(data_bal_p_psd)/2-1)*fs/length(data_bal_p_psd); 
figure;plot(data_bal_p_psd);title('C点功率谱');xlabel('频率/mHz');ylabel('幅值/V');
%% 解调器

reciv_i=real(data_bal_p);
reciv_q=imag(data_bal_p);
receive_data=[];
for j=1:length(reciv_i)
    receive_data=[receive_data sign(reciv_i(j)) sign(reciv_q(j))]; 
end
receive_data=(receive_data+1)/2;
biterr(receive_data,signal1)


toc