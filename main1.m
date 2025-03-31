clear all
close all
warning off

%% 实验设置
NMonte=500; % 蒙特卡洛试验次数
SnrdBMat=[-30:5:10]; % 信噪比矩阵
NSnr=length(SnrdBMat);
SnapshotsF=3; % 频域快拍数
DOAS=[-15.5,10.3]; % 信源方向
NSignals=length(DOAS); % 信源个数

%% 阵列设置
c=344; % 波速
dd=0.1; % ULA阵元间距
NEle=10; % ULA阵元数
LocArray=[0:dd:dd*(NEle-1)]'; % 阵列位置

%% 信号采样设置
Fs=8000; % 采样频率
Nfft=128; % FFT点数
Framelength=1; % 帧长512，0%overlap
Snapshots=Nfft*SnapshotsF*Framelength; % 时域快拍数

Fl=1200; % 感兴趣频率下限
Fh=2200; % 感兴趣频率上限
IFl=ceil(Fl/Fs*Nfft); % 感兴趣频点序号下限
IFh=floor(Fh/Fs*Nfft); % 感兴趣频点序号上限
NBand=IFh-IFl+1; % 感兴趣子带个数
FreqMat=[IFl:IFh]/Nfft*Fs; % 感兴趣频点

%% 统计结果暂存
Time_NarrowMUSIC=zeros(1,NSnr);
Time_IssmMUSIC=zeros(1,NSnr);
Time_IssmMUSIC2=zeros(1,NSnr);
Time_CssmMUSIC=zeros(1,NSnr);
Time_mTOPS=zeros(1,NSnr);
Time_SBL=zeros(1,NSnr);
Time_proposed=zeros(1,NSnr);

%% 蒙特卡洛实验
for iMonte=1:NMonte
    for iSnr=1:NSnr
%     iSnr=1;
    %% 生成阵列观测
        Snr=10^(SnrdBMat(iSnr)/10); % 信噪比（功率）
%         RandSeed=round(rand*(length(data1)-Snapshots)); % 随机选取部分音频
%         Signals=repmat(data1(RandSeed:RandSeed+Snapshots-1)',NSignals); % 选取的信源
%         Noise=randn(NEle,Snapshots);
        Signals=1/sqrt(2)*(randn(NSignals,Snapshots)+1j*randn(NSignals,Snapshots));
%         Signals=1/sqrt(2)*(randn(1,Snapshots)+1j*randn(1,Snapshots));
%         Signals=repmat(Signals,NSignals,1);
        Noise=1/sqrt(2)*(randn(NEle,Snapshots)+1j*randn(NEle,Snapshots));
        
        
        [bI,aI]=butter(8,[Fl/(Fs/2),Fh/(Fs/2) ]);
        for i=1:NSignals
            for ii=1:SnapshotsF
                temp1=filter(bI,aI,Signals(i,1+(ii-1)*Framelength*Nfft:ii*Framelength*Nfft));
                temp1=temp1/sqrt(trace(temp1*temp1')/Snapshots)*sqrt(Snr);
                temp2=fft(temp1,Nfft);
                for iii=1:NBand
                    SignalsF{iii}(i,ii)=temp2(IFl+iii-1);
                end
            end
        end
        for i=1:NEle
            for ii=1:SnapshotsF
        %         temp3=filter(bI,aI,Noise(i,1+(ii-1)*Framelength*Nfft:ii*Framelength*Nfft));
                temp4=fft(Noise(i,1+(ii-1)*Framelength*Nfft:ii*Framelength*Nfft),Nfft);
                for iii=1:NBand
                    NoiseF{iii}(i,ii)=temp4(IFl+iii-1);
                end
            end
        end
        
        
        
        
        
        
        
        % sound(Signals,Fs); % 试听
        % NoisySignals=awgn(Signals,Snr,'measured');
        % sound(NoisySignals,Fs); % 试听
       %% 真实导向矢量
        for i=1:NBand
            for ii=1:NSignals
                f=FreqMat(i); % 频率
                AReal{i}(:,ii)=exp(-1j*2*pi*f*LocArray'/c*sin(pi/180*DOAS(ii))); % 真实导向矢量
            end
        end
       %% 信源STFT
% %         PowerSignal=0;
%         for i=1:NSignals
%             for ii=1:SnapshotsF
%                 temp=fft(Signals(i,1+(ii-1)*Framelength*Nfft:ii*Framelength*Nfft),Nfft);
%                 % plot([1:Nfft]/Nfft*Fs,abs(temp)) % 频谱图
%                 for iii=1:NBand
%                     SignalsF{iii}(i,ii)=temp(iii+IFl-1);
% %                     PowerSignal=PowerSignal+abs(SignalsF{iii}(i,ii))^2/SnapshotsF/NSignals/NBand;
%                 end
%             end
%         end
%         for i=1:NBand
%             SignalsF{i}=SignalsF{i}/sqrt(PowerSignal)*sqrt(Snr); % 放大到指定的倍数
%         end
        for i=1:NBand
            SignalPowerSpectrum(i)=abs(trace(SignalsF{i}*SignalsF{i}'))/NSignals/SnapshotsF; % 感兴趣频带的信号功率谱
        end
        % plot(FreqMat,SignalPowerSpectrum)
       %% 噪声STFT
%         PowerNoise=0;
%         for i=1:NEle
%             for ii=1:SnapshotsF
%                 temp=fft(Noise(i,1+(ii-1)*Framelength*Nfft:ii*Framelength*Nfft),Nfft);
%                 for iii=1:NBand
%                     NoiseF{iii}(i,ii)=temp(iii+IFl-1);
% %                     PowerNoise=PowerNoise+abs(NoiseF{iii}(i,ii))^2/NEle/SnapshotsF/NBand;
%                 end
%             end
%         end
% %         for i=1:NBand
% %             NoiseF{i}=NoiseF{i}/sqrt(PowerNoise);
% %         end
        for i=1:NBand
            NoisePowerSpectrum(i)=abs(trace(NoiseF{i}*NoiseF{i}'))/NEle/SnapshotsF; % 感兴趣频带的噪声功率谱
        end
        
%         figure
%         plot(FreqMat,SignalPowerSpectrum)
%         hold on
%         plot(FreqMat,NoisePowerSpectrum)
       %% 阵列观测
        for i=1:NBand
            X{i}=AReal{i}*SignalsF{i};
            Xn{i}=X{i}+NoiseF{i};
        end
        

            
        
        
        
        %% DOA估计
        GridRes=2; % 搜索网格精度
        
        %% 窄带MUSIC
        tic
        PickBand=9; % 选择第一个子带
        [DOA_NarrowMUSIC,P_NarrowMUSIC]=NarrowMUSIC(Xn, LocArray, FreqMat, NSignals, c, PickBand, GridRes);
        Time_NarrowMUSIC(iSnr)=Time_NarrowMUSIC(iSnr)+toc/NMonte;
        
        %% IssmCBF
%         [DOA_IssmCBF,P_IssmCBF]=IssmCBF(Xn, LocArray, FreqMat, NSignals, c, GridRes);
        
        %% IssmMUSIC（算术平均）
        tic
        [DOA_IssmMUSIC,P_IssmMUSIC]=IssmMUSIC(Xn, LocArray, FreqMat, NSignals, c, GridRes);
        Time_IssmMUSIC(iSnr)=Time_IssmMUSIC(iSnr)+toc/NMonte;
        
        %% IssmMUSIC（几何平均）
        tic
        [DOA_IssmMUSIC2,P_IssmMUSIC2]=IssmMUSIC2(Xn, LocArray, FreqMat, NSignals, c, GridRes);
        Time_IssmMUSIC2(iSnr)=Time_IssmMUSIC2(iSnr)+toc/NMonte;
        
        %% CssmMUSIC（Rss）
        tic
        [DOA_CssmMUSIC,P_CssmMUSIC]=CssmMUSIC(Xn, LocArray, FreqMat, NSignals, c, GridRes);
        Time_CssmMUSIC(iSnr)=Time_CssmMUSIC(iSnr)+toc/NMonte;  
        %% mTOPS
        tic
        [DOA_mTOPS,P_mTOPS]=mTOPS(Xn, LocArray, FreqMat, NSignals, c, GridRes);
        Time_mTOPS(iSnr)=Time_mTOPS(iSnr)+toc/NMonte;  
        %% SBL
        tic
        [DOA_SBL,P_SBL]=SBL(Xn, LocArray, FreqMat, NSignals, c, GridRes);
        Time_SBL(iSnr)=Time_SBL(iSnr)+toc/NMonte;  
        %% Proposed
        tic
        [DOA_proposed,P_proposed]=proposed(Xn, LocArray, FreqMat, NSignals, c, GridRes);
        Time_proposed(iSnr)=Time_proposed(iSnr)+toc/NMonte;  
        
        %% 统计误差
        Err_NarrowMUSIC(iSnr,iMonte,:)=abs(DOAS-DOA_NarrowMUSIC);
        Err_IssmMUSIC(iSnr,iMonte,:)=abs(DOAS-DOA_IssmMUSIC);
        Err_IssmMUSIC2(iSnr,iMonte,:)=abs(DOAS-DOA_IssmMUSIC2);
        Err_CssmMUSIC(iSnr,iMonte,:)=abs(DOAS-DOA_CssmMUSIC);
        Err_mTOPS(iSnr,iMonte,:)=abs(DOAS-DOA_mTOPS);
        Err_SBL(iSnr,iMonte,:)=abs(DOAS-DOA_SBL);
        Err_proposed(iSnr,iMonte,:)=abs(DOAS-DOA_proposed);
        
%         MSE_NarrowMUSIC(iSnr)=MSE_NarrowMUSIC(iSnr)+(DOAS-DOA_NarrowMUSIC)*(DOAS-DOA_NarrowMUSIC)'/NSignals/NMonte;
%         MSE_IssmMUSIC(iSnr)=MSE_IssmMUSIC(iSnr)+(DOAS-DOA_IssmMUSIC)*(DOAS-DOA_IssmMUSIC)'/NSignals/NMonte;
%         MSE_IssmMUSIC2(iSnr)=MSE_IssmMUSIC2(iSnr)+(DOAS-DOA_IssmMUSIC2)*(DOAS-DOA_IssmMUSIC2)'/NSignals/NMonte;
%         MSE_CssmMUSIC(iSnr)=MSE_CssmMUSIC(iSnr)+(DOAS-DOA_CssmMUSIC)*(DOAS-DOA_CssmMUSIC)'/NSignals/NMonte;
%         MSE_MTOPS(iSnr)=MSE_MTOPS(iSnr)+(DOAS-DOA_mTOPS)*(DOAS-DOA_mTOPS)'/NSignals/NMonte;
%         MSE_SBL(iSnr)=MSE_SBL(iSnr)+(DOAS-DOA_SBL)*(DOAS-DOA_SBL)'/NSignals/NMonte;
%         MSE_proposed(iSnr)=MSE_proposed(iSnr)+(DOAS-DOA_proposed)*(DOAS-DOA_proposed)'/NSignals/NMonte;
        
        
        
    end
    disp(['进度：',num2str(iMonte/NMonte*100),'%'])
end


[~,NMonte,~]=size(Err_proposed);
%% 统计RMSE
for iSnr=1:NSnr
    RMSE_NarrowMUSIC(iSnr)=sqrt(sum(sum(Err_NarrowMUSIC(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_IssmMUSIC(iSnr)=sqrt(sum(sum(Err_IssmMUSIC(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_IssmMUSIC2(iSnr)=sqrt(sum(sum(Err_IssmMUSIC2(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_CssmMUSIC(iSnr)=sqrt(sum(sum(Err_CssmMUSIC(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_mTOPS(iSnr)=sqrt(sum(sum(Err_mTOPS(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_SBL(iSnr)=sqrt(sum(sum(Err_SBL(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_proposed(iSnr)=sqrt(sum(sum(Err_proposed(iSnr,:,:).^2))/NMonte/NSignals);
    
end

%% 统计成功率
th=2;
for iSnr=1:NSnr
    Success_NarrowMUSIC(iSnr)=sum(max(Err_NarrowMUSIC(iSnr,:,:),[],3)<=th)/NMonte;
    Success_IssmMUSIC(iSnr)=sum(max(Err_IssmMUSIC(iSnr,:,:),[],3)<=th)/NMonte;
%     Success_IssmMUSIC2(iSnr)=sum(max(Err_IssmMUSIC2(iSnr,:,:),[],3)<=th)/NMonte;
    Success_CssmMUSIC(iSnr)=sum(max(Err_CssmMUSIC(iSnr,:,:),[],3)<=th)/NMonte;
    Success_mTOPS(iSnr)=sum(max(Err_mTOPS(iSnr,:,:),[],3)<=th)/NMonte;
    Success_SBL(iSnr)=sum(max(Err_SBL(iSnr,:,:),[],3)<=th)/NMonte;
    Success_proposed(iSnr)=sum(max(Err_proposed(iSnr,:,:),[],3)<=th)/NMonte;
end

 %% 绘制RMSE图
figure
semilogy(SnrdBMat,RMSE_NarrowMUSIC,'-.')
hold on
semilogy(SnrdBMat,RMSE_IssmMUSIC,'--')
% plot(SnrdBMat,RMSE_CssmMUSIC,'-')
semilogy(SnrdBMat,RMSE_mTOPS,'-.')
semilogy(SnrdBMat,RMSE_SBL,'--')
semilogy(SnrdBMat,RMSE_proposed,'-')
grid on
legend('窄带MUSIC','ISSM MUSIC','mTOPS','宽带SBL','提出的方法')
% legend('窄带MUSIC','ISSM MUSIC','CSSM MUSIC','mTOPS','宽带SBL','提出的方法')
xlabel('信噪比（dB）')
ylabel('均方根误差（°）')
                    
%% 绘制成功率图
figure
plot(SnrdBMat,Success_NarrowMUSIC,'-.')
hold on
plot(SnrdBMat,Success_IssmMUSIC,'--')
% plot(SnrdBMat,Success_CssmMUSIC,'-')
plot(SnrdBMat,Success_mTOPS,'-.')
plot(SnrdBMat,Success_SBL,'--')
plot(SnrdBMat,Success_proposed,'-')
grid on
legend('窄带MUSIC','ISSM MUSIC','mTOPS','宽带SBL','提出的方法')
% legend('窄带MUSIC','ISSM MUSIC','CSSM MUSIC','mTOPS','宽带SBL','提出的方法')
xlabel('信噪比（dB）')
ylabel('估计成功率')

%% 绘制运行时间图
figure
semilogy(SnrdBMat,Time_NarrowMUSIC,'-.')
hold on
semilogy(SnrdBMat,Time_IssmMUSIC,'--')
% plot(SnrdBMat,Success_CssmMUSIC,'-')
semilogy(SnrdBMat,Time_mTOPS,'-.')
semilogy(SnrdBMat,Time_SBL,'--')
semilogy(SnrdBMat,Time_proposed,'-')
grid on
legend('窄带MUSIC','ISSM MUSIC','mTOPS','宽带SBL','提出的方法')
% legend('窄带MUSIC','ISSM MUSIC','CSSM MUSIC','mTOPS','宽带SBL','提出的方法')
xlabel('信噪比（dB）')
ylabel('平均运行时间（s）')       
