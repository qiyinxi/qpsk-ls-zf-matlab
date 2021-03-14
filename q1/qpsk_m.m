clc
clear
close all
no_seq = 1000;%%如果出眼图的话这里最好设成1k，不然运算太慢
Rb=100;                                                 %信息速率，单位mbit/s
Tb=1/Rb;
Ts=2*Tb;
fc=2*Rb;                                                %载波速率
A=1;
constant=4;                                             %抽样常数
fs=constant*fc;                                         %抽样频率
ts=1/fs;                                                %抽样间隔
N=Tb/ts;                                                %一个bit内样点数
a=randi([0 1],1,no_seq);                                %bit源
No_sample=length(a)*N;                                  %样点总数

%% 调制前内插与A点波形显示
a1=a(:);
a2=a1*ones(1,N);
a3=a2';
a4=reshape(a3,1,[]);
t=[0:No_sample-1]*ts;

figure(1)
stem(a);axis([0 20 -0.2 1.2]);title('bit源（A）')

%% 产生B方式的QPSK信号并显示B点图形
mapping=[5*pi/4 3*pi/4 7*pi/4 pi/4];
gl=[-1-i -1+i 1-i 1+i];
qpsk=[];
qpsk_gl=[];
for k=1:2:no_seq
    index=0;
    for j=k:k+1
        index=2*index+a(j);
    end
    index=index+1;
    theta=mapping(index);
    thegl=gl(index);
    qpsk_gl=[qpsk_gl thegl];
    qpsk=[qpsk,A*cos(2*pi*fc*t((k-1)*N+1:(k+1)*N)+theta)];
end
scatterplot(qpsk_gl);title('调制好的星座图');
figure(3)
subplot(311)
plot(t*1000,qpsk);grid;axis([0 100 -1.2 1.2]);title('过调制信号(B)')
xlabel('t/ns');ylabel('幅度/V');

df=fs/2000;
f=[0:df:df*(2000-1)]-fs/2;
QPSK=fft(qpsk,2000)/fs;
QPSK0=abs(fftshift(QPSK));
figure(3)
subplot(312)
plot(f,QPSK0);grid;
axis([-400 400 0 0.2]);title('过调制信号频谱')
ylabel('频谱/V');xlabel('频率/mHz');

subplot(313)
qpsk0_psd=QPSK0.*conj(QPSK0);
plot(f,qpsk0_psd);grid;
axis([-400 400 0 0.1]);title('过调制信号功率谱')
ylabel('频谱/V');xlabel('频率/mHz');


%% 上采样与脉冲成型滤波
% 脉冲形成：升余弦滤波器
r = 0.25; %滚降系数
delay = 8; %延时
IPOINT = 4; %内插倍数
h_rcos = rcosdesign(r, delay, IPOINT, 'sqrt'); %生成时域响应h
len_fir = length(h_rcos);
len_qpsk = length(qpsk);
qpsk_unfir=qpsk;
qpsk1 = upfirdn(qpsk_unfir, h_rcos, 4 ); %对信号进行匹配滤波
qpsk=qpsk1(fix(len_fir/2+1):fix((len_fir/2)+len_qpsk*4));
% C点图形
figure(4)
subplot(211)
t1=[0:No_sample*4-1]*ts/4;
plot(t1*1000,qpsk);grid;axis([0 100 -1.2 1.2]);title('脉冲成型信号(C)')
xlabel('t/ns');ylabel('幅度/V');

df=fs/2000;
f=[0:df*4:df*(2000-1)*4]-fs*4/2;
QPSK=fft(qpsk,2000)/fs;
QPSK1=abs(fftshift(QPSK));
figure(4)
subplot(212)
plot(f,QPSK1);grid;
axis([-400 400 0 0.2]);title('脉冲成型信号频谱')
ylabel('频谱/V');xlabel('频率/mHz');
eyediagram(qpsk,6);title('脉冲成型（无噪声原始）眼图');

%% 误码率仿真
Eb=(1/2)*(A^2)*Tb;
snr_in_db1=-5:10;
for k=1:length(snr_in_db1)
%% 信号过高斯噪声信道
    snr=10^(snr_in_db1(k)/10);
    n0=Eb/snr;
    delta=sqrt(n0*(fs/2));
    noise=delta*randn(1,No_sample*4);
    rqpsk_unfir=qpsk+noise;
    if (k==10)
       eyediagram(rqpsk_unfir,6);title('ebn0=5情况下加噪声后眼图') 
    end
    if (k==15)
       eyediagram(rqpsk_unfir,6);title('ebn0=10情况下加噪声后眼图') 
    end
    
    len_qpsk = length(rqpsk_unfir);
    rqpsk_fir = upfirdn(rqpsk_unfir, h_rcos, 1); %对信号进行匹配滤波
    len_fir = fix(length(h_rcos)/4);
    len_rqpsk = length(rqpsk_unfir)/4;
    rqpsk_fir_f=rqpsk_fir(len_fir*2+1:length(rqpsk_unfir)+len_fir*2);
    if (k==10)
        
        figure
        subplot(311)
        plot(t1*1000,rqpsk_fir_f);grid;axis([0 40 -5 5]);title('ebn0=5情况下匹配滤波信号(D)')
        xlabel('t/ns');ylabel('幅度/V');
        
        subplot(312)
        f_rfir=[0:df*4:df*(2000-1)*4]-fs*4/2;
        QPSK_rfir=fft(rqpsk_fir_f,2000)/fs;
        QPSK1_rfir=abs(fftshift(QPSK_rfir));
        plot(f,QPSK1_rfir);grid;
        axis([-400 400 0 0.7]);title('ebn0=5情况下匹配滤波信号频谱')
        ylabel('频谱/V');xlabel('频率/mHz');
        
        subplot(313)
        qpsk1_rfir_psd=QPSK1_rfir.*conj(QPSK1_rfir);
        plot(f,qpsk1_rfir_psd);grid;
        axis([-400 400 0 0.1]);title('ebn0=5情况下匹配滤波信号功率谱')
        ylabel('频谱/V');xlabel('频率/mHz');
        
        
    end
    rqpsk=rqpsk_fir(17:4:end-16);%降采样
    if (k==10)
        figure
        plot(t*1000,rqpsk);axis([0 40 -5 5]);title('ebn0=5情况下下采样信号(E)')
        xlabel('t/ns');ylabel('幅度/V');
        eyediagram(rqpsk,6);title('ebn0=5情况下下采样信号眼图')
    end
    if (k==15)
        eyediagram(rqpsk,6);title('ebn0=10情况下下采样信号眼图')
    end
    
    
%     rqpsk=rqpsk_unfir;
%% 信号解调  
    Xc1_unfil=rqpsk.*cos(2*pi*fc*t);
    Xs1_unfil=-rqpsk.*sin(2*pi*fc*t);
    
%% 低通滤波
    Srrc=srrc;
    S_n=Srrc.Numerator;
    len_srrc=fix(length(S_n)/2);
    len_Xc1_unfil=length(Xc1_unfil);
    Xc1_fil=conv(S_n,Xc1_unfil);
    Xs1_fil=conv(S_n,Xs1_unfil);
    Xc1=Xc1_fil(len_srrc+1:len_Xc1_unfil+len_srrc);
    Xs1=Xs1_fil(len_srrc+1:len_Xc1_unfil+len_srrc);
    
    

%% 抽样判决    
    Xc2=reshape(Xc1,Ts/ts,[]);
    Xc=sum(Xc2)*ts;
%     Xc=16*Xc2(9,:)*ts;
    
    Xs2=reshape(Xs1,Ts/ts,[]);
    Xs=sum(Xs2)*ts;   
%     Xs=16*Xs2(9,:)*ts;%%单纯取一点抽样会导致误码率上升
    if(k==10)
        scatterplot(200*(Xc+1i*Xs));title('ebn0=5情况下接收端星座图');
    end
    if(k==15)
        scatterplot(200*(Xc+1i*Xs));title('ebn0=10情况下接收端星座图');
    end
    decis=[];
    for m=1:no_seq/2
        theta=mod(angle(Xc(m)+1i*Xs(m)),2*pi);
        if(theta>3*pi/2)
            decis=[decis 1 0];
        elseif(theta<pi/2)
            decis=[decis 1 1];
        elseif(theta<pi)
            decis=[decis 0 1];
        else
            decis=[decis 0 0];
        end
    end
    
%% 求ber
    biterror=0;
    for n=1:no_seq
        if(a(n)~=decis(n))
            biterror=biterror+1;
        end
    end
    pb(k)=biterror/no_seq;
end
snr_in_db2=-5:0.1:10;
for i=1:length(snr_in_db2)
    snr=10^(snr_in_db2(i)/10);
    theo_pb(i)=(1/2)*erfc(sqrt(snr));
end

figure
semilogy(snr_in_db1,pb,'*');hold on
semilogy(snr_in_db2,theo_pb);grid;
title('QPSK的ber曲线');
xlabel('Eb/n0(dB)');ylabel('P_e');
legend('仿真误比特率','理论误比特率');








