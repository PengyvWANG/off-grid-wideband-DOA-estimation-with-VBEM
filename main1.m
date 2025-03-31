clear all
close all
warning off

%% ʵ������
NMonte=500; % ���ؿ����������
SnrdBMat=[-30:5:10]; % ����Ⱦ���
NSnr=length(SnrdBMat);
SnapshotsF=3; % Ƶ�������
DOAS=[-15.5,10.3]; % ��Դ����
NSignals=length(DOAS); % ��Դ����

%% ��������
c=344; % ����
dd=0.1; % ULA��Ԫ���
NEle=10; % ULA��Ԫ��
LocArray=[0:dd:dd*(NEle-1)]'; % ����λ��

%% �źŲ�������
Fs=8000; % ����Ƶ��
Nfft=128; % FFT����
Framelength=1; % ֡��512��0%overlap
Snapshots=Nfft*SnapshotsF*Framelength; % ʱ�������

Fl=1200; % ����ȤƵ������
Fh=2200; % ����ȤƵ������
IFl=ceil(Fl/Fs*Nfft); % ����ȤƵ���������
IFh=floor(Fh/Fs*Nfft); % ����ȤƵ���������
NBand=IFh-IFl+1; % ����Ȥ�Ӵ�����
FreqMat=[IFl:IFh]/Nfft*Fs; % ����ȤƵ��

%% ͳ�ƽ���ݴ�
Time_NarrowMUSIC=zeros(1,NSnr);
Time_IssmMUSIC=zeros(1,NSnr);
Time_IssmMUSIC2=zeros(1,NSnr);
Time_CssmMUSIC=zeros(1,NSnr);
Time_mTOPS=zeros(1,NSnr);
Time_SBL=zeros(1,NSnr);
Time_proposed=zeros(1,NSnr);

%% ���ؿ���ʵ��
for iMonte=1:NMonte
    for iSnr=1:NSnr
%     iSnr=1;
    %% �������й۲�
        Snr=10^(SnrdBMat(iSnr)/10); % ����ȣ����ʣ�
%         RandSeed=round(rand*(length(data1)-Snapshots)); % ���ѡȡ������Ƶ
%         Signals=repmat(data1(RandSeed:RandSeed+Snapshots-1)',NSignals); % ѡȡ����Դ
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
        
        
        
        
        
        
        
        % sound(Signals,Fs); % ����
        % NoisySignals=awgn(Signals,Snr,'measured');
        % sound(NoisySignals,Fs); % ����
       %% ��ʵ����ʸ��
        for i=1:NBand
            for ii=1:NSignals
                f=FreqMat(i); % Ƶ��
                AReal{i}(:,ii)=exp(-1j*2*pi*f*LocArray'/c*sin(pi/180*DOAS(ii))); % ��ʵ����ʸ��
            end
        end
       %% ��ԴSTFT
% %         PowerSignal=0;
%         for i=1:NSignals
%             for ii=1:SnapshotsF
%                 temp=fft(Signals(i,1+(ii-1)*Framelength*Nfft:ii*Framelength*Nfft),Nfft);
%                 % plot([1:Nfft]/Nfft*Fs,abs(temp)) % Ƶ��ͼ
%                 for iii=1:NBand
%                     SignalsF{iii}(i,ii)=temp(iii+IFl-1);
% %                     PowerSignal=PowerSignal+abs(SignalsF{iii}(i,ii))^2/SnapshotsF/NSignals/NBand;
%                 end
%             end
%         end
%         for i=1:NBand
%             SignalsF{i}=SignalsF{i}/sqrt(PowerSignal)*sqrt(Snr); % �Ŵ�ָ���ı���
%         end
        for i=1:NBand
            SignalPowerSpectrum(i)=abs(trace(SignalsF{i}*SignalsF{i}'))/NSignals/SnapshotsF; % ����ȤƵ�����źŹ�����
        end
        % plot(FreqMat,SignalPowerSpectrum)
       %% ����STFT
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
            NoisePowerSpectrum(i)=abs(trace(NoiseF{i}*NoiseF{i}'))/NEle/SnapshotsF; % ����ȤƵ��������������
        end
        
%         figure
%         plot(FreqMat,SignalPowerSpectrum)
%         hold on
%         plot(FreqMat,NoisePowerSpectrum)
       %% ���й۲�
        for i=1:NBand
            X{i}=AReal{i}*SignalsF{i};
            Xn{i}=X{i}+NoiseF{i};
        end
        

            
        
        
        
        %% DOA����
        GridRes=2; % �������񾫶�
        
        %% խ��MUSIC
        tic
        PickBand=9; % ѡ���һ���Ӵ�
        [DOA_NarrowMUSIC,P_NarrowMUSIC]=NarrowMUSIC(Xn, LocArray, FreqMat, NSignals, c, PickBand, GridRes);
        Time_NarrowMUSIC(iSnr)=Time_NarrowMUSIC(iSnr)+toc/NMonte;
        
        %% IssmCBF
%         [DOA_IssmCBF,P_IssmCBF]=IssmCBF(Xn, LocArray, FreqMat, NSignals, c, GridRes);
        
        %% IssmMUSIC������ƽ����
        tic
        [DOA_IssmMUSIC,P_IssmMUSIC]=IssmMUSIC(Xn, LocArray, FreqMat, NSignals, c, GridRes);
        Time_IssmMUSIC(iSnr)=Time_IssmMUSIC(iSnr)+toc/NMonte;
        
        %% IssmMUSIC������ƽ����
        tic
        [DOA_IssmMUSIC2,P_IssmMUSIC2]=IssmMUSIC2(Xn, LocArray, FreqMat, NSignals, c, GridRes);
        Time_IssmMUSIC2(iSnr)=Time_IssmMUSIC2(iSnr)+toc/NMonte;
        
        %% CssmMUSIC��Rss��
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
        
        %% ͳ�����
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
    disp(['���ȣ�',num2str(iMonte/NMonte*100),'%'])
end


[~,NMonte,~]=size(Err_proposed);
%% ͳ��RMSE
for iSnr=1:NSnr
    RMSE_NarrowMUSIC(iSnr)=sqrt(sum(sum(Err_NarrowMUSIC(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_IssmMUSIC(iSnr)=sqrt(sum(sum(Err_IssmMUSIC(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_IssmMUSIC2(iSnr)=sqrt(sum(sum(Err_IssmMUSIC2(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_CssmMUSIC(iSnr)=sqrt(sum(sum(Err_CssmMUSIC(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_mTOPS(iSnr)=sqrt(sum(sum(Err_mTOPS(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_SBL(iSnr)=sqrt(sum(sum(Err_SBL(iSnr,:,:).^2))/NMonte/NSignals);
    RMSE_proposed(iSnr)=sqrt(sum(sum(Err_proposed(iSnr,:,:).^2))/NMonte/NSignals);
    
end

%% ͳ�Ƴɹ���
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

 %% ����RMSEͼ
figure
semilogy(SnrdBMat,RMSE_NarrowMUSIC,'-.')
hold on
semilogy(SnrdBMat,RMSE_IssmMUSIC,'--')
% plot(SnrdBMat,RMSE_CssmMUSIC,'-')
semilogy(SnrdBMat,RMSE_mTOPS,'-.')
semilogy(SnrdBMat,RMSE_SBL,'--')
semilogy(SnrdBMat,RMSE_proposed,'-')
grid on
legend('խ��MUSIC','ISSM MUSIC','mTOPS','���SBL','����ķ���')
% legend('խ��MUSIC','ISSM MUSIC','CSSM MUSIC','mTOPS','���SBL','����ķ���')
xlabel('����ȣ�dB��')
ylabel('���������㣩')
                    
%% ���Ƴɹ���ͼ
figure
plot(SnrdBMat,Success_NarrowMUSIC,'-.')
hold on
plot(SnrdBMat,Success_IssmMUSIC,'--')
% plot(SnrdBMat,Success_CssmMUSIC,'-')
plot(SnrdBMat,Success_mTOPS,'-.')
plot(SnrdBMat,Success_SBL,'--')
plot(SnrdBMat,Success_proposed,'-')
grid on
legend('խ��MUSIC','ISSM MUSIC','mTOPS','���SBL','����ķ���')
% legend('խ��MUSIC','ISSM MUSIC','CSSM MUSIC','mTOPS','���SBL','����ķ���')
xlabel('����ȣ�dB��')
ylabel('���Ƴɹ���')

%% ��������ʱ��ͼ
figure
semilogy(SnrdBMat,Time_NarrowMUSIC,'-.')
hold on
semilogy(SnrdBMat,Time_IssmMUSIC,'--')
% plot(SnrdBMat,Success_CssmMUSIC,'-')
semilogy(SnrdBMat,Time_mTOPS,'-.')
semilogy(SnrdBMat,Time_SBL,'--')
semilogy(SnrdBMat,Time_proposed,'-')
grid on
legend('խ��MUSIC','ISSM MUSIC','mTOPS','���SBL','����ķ���')
% legend('խ��MUSIC','ISSM MUSIC','CSSM MUSIC','mTOPS','���SBL','����ķ���')
xlabel('����ȣ�dB��')
ylabel('ƽ������ʱ�䣨s��')       
