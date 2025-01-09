function [DOA_proposed,P_proposed]=proposed(Xn, LocArray, FreqMat, NSignals, c, GridRes)
NBand=length(Xn);
SearchGrid=[-90:GridRes:90];
NGrid=length(SearchGrid);
SpatialSpectrum=ones(1,NGrid);

epsilon=1e-4;
[NEle,SnapshotsF]=size(Xn{1});

MaxIter=30;
P_proposed=zeros(NGrid,1);betaFinal=zeros(1,NGrid);
for iBand=1:NBand
%     Y{iBand}=xn{iBand}/max(max(abs(xn{iBand})));
    Y{iBand}=Xn{iBand};
    f(iBand)=FreqMat(iBand);
    A{iBand}=exp(-1j*2*pi*f(iBand)*LocArray/c*sin(SearchGrid*pi/180)); % 观测矩阵
    B{iBand}=exp(-1j*2*pi*f(iBand)*LocArray/c*sin(SearchGrid*pi/180)).*(-1j*2*pi*f(iBand)*LocArray/c).*cos(SearchGrid*pi/180)*pi/180;
%     alpha(:,iBand)=ones(1,NGrid)*1e-11; % hyperparameter
%     alpha0(iBand)=1e-1; % hyperparameter
% %     mu{iBand}=pinv(A{iBand})*Y{iBand};
    mu{iBand}=zeros(NEle,SnapshotsF);
end
iter=0;
z=ones(1,NGrid)*0.5;
w=ones(1,NGrid)*0.5;
% z=ones(1,NGrid)*(NEle-1)/NGrid;
% w=ones(1,NGrid)*(NEle-1)/NGrid;
beta=zeros(1,NGrid);

flag=1;


for iBand=1:NBand
    Rxn{iBand}=Xn{iBand}*Xn{iBand}'/SnapshotsF;
    for ii=1:NGrid
        alpha(ii,iBand)=NEle^2/real(A{iBand}(:,ii)'*Rxn{iBand}*A{iBand}(:,ii));
    end
end

for iBand=1:NBand
    alpha0(iBand)=1/real(trace(Rxn{iBand}/NEle));
end

% tic
% while iter<MaxIter && flag>epsilon && max(z)<1 && min(z)>0
while iter<MaxIter && flag>epsilon
    iter=iter+1;
    zOld=z;
    for iBand=1:NBand
        Phi{iBand}=A{iBand}+B{iBand}*diag(beta);
        sigma{iBand}=real(1./(alpha(:,iBand)+z'.*alpha0(iBand).*diag(Phi{iBand}'*Phi{iBand})))';
        mu{iBand}=z'*alpha0(iBand).*sigma{iBand}'.*Phi{iBand}'*Y{iBand};
    end

    for iGrid=1:NGrid
        tempz1(iGrid)=0;
        tempz0(iGrid)=0;
        for iBand=1:NBand
            tempz1(iGrid)=tempz1(iGrid)-alpha0(iBand)*(trace(Y{iBand}*Y{iBand}')-2*real(trace(Y{iBand}*(Phi{iBand}(:,iGrid)*mu{iBand}(iGrid,:))'))+...
                NEle*(mu{iBand}(iGrid,:)*mu{iBand}(iGrid,:)'+SnapshotsF*sigma{iBand}(iGrid)));
            tempz0(iGrid)=tempz0(iGrid)-alpha0(iBand)*trace(Y{iBand}*Y{iBand}');
        end
        z(iGrid)=w(iGrid)/((1-w(iGrid))*exp(tempz0(iGrid)-tempz1(iGrid))+w(iGrid));
%         z(iGrid)=w(iGrid)*exp(tempz1)/(w(iGrid)*exp(tempz1)+(1-w(iGrid))*exp(tempz0));
    end   
    z(z<1e-4)=1e-4;
%     z(z>1-1e-4)=1-1e-4;
%     z=z/sum(z);
    w=z;

    for iBand=1:NBand
        temptemp(:,iBand)=real(diag(mu{iBand}*mu{iBand}')/SnapshotsF+sigma{iBand}');
        alpha(:,iBand)=1./temptemp(:,iBand);
    end
    

    for iBand=1:NBand
        temp1(iBand)=NGrid*trace(Y{iBand}*Y{iBand}')-2*real(trace(Y{iBand}*(Phi{iBand}*diag(z)*mu{iBand})'))+trace(Phi{iBand}*diag(z)*mu{iBand}*...
            mu{iBand}'*Phi{iBand}')+trace(diag(z)*diag(sigma{iBand}));
        alpha0(iBand)=real(NEle*SnapshotsF*NGrid/temp1(iBand));
    end
%     flag=abs(max(z)-max(zOld))/abs(max(zOld));
    flag=norm(z-zOld)/norm(zOld);

    for iGrid=1:NGrid
        temp1=0;temp2=0;temp3=0;
        for iBand=1:NBand
%             for in=1:N
%                 temp1=temp1+real((Y{iBand}(:,in)-A{iBand}(:,iGrid))'*B{iBand}(:,iGrid)*mu{iBand}(iGrid,in));
%                 temp2=temp2+real(B{iBand}(:,iGrid)'*B{iBand}(:,iGrid)*(mu{iBand}(iGrid,in)'*mu{iBand}(iGrid,in)+M*sigma{iBand}(iGrid)));
            temp1=temp1+real(trace((Y{iBand}-A{iBand}(:,iGrid))*(B{iBand}(:,iGrid)*mu{iBand}(iGrid,:))'));
%             temp2=temp2+real(trace(A{iBand}(:,iGrid)'*B{iBand}(:,iGrid)*(mu{iBand}(iGrid,:))));
            temp2=temp2+real(B{iBand}(:,iGrid)'*B{iBand}(:,iGrid)*(mu{iBand}(iGrid,:)*mu{iBand}(iGrid,:)'+NEle*SnapshotsF*sigma{iBand}(iGrid)));
%             end
        end
        beta(iGrid)=temp1/temp2;
%         beta(iGrid)=(temp1-temp2)/temp3;
        if beta(iGrid)<-GridRes/2
            beta(iGrid)=-GridRes/2;
        elseif beta(iGrid)>GridRes/2
            beta(iGrid)=GridRes/2;
        end
    end
%     plot(Grid,w);
%     figure
%     plot(Grid,alpha.^-1)

end   
% toc
% pFinal=zOld;
% pFinal=z;
P_proposed=z'.*sum(alpha.^(-1),2);
% P_proposed=sum(alpha.^(-1),2);

% plot(SearchGrid,P_proposed)
% hold on

% beta=zeros(1,NGrid);


[x1,y1]=find(beta>90/(NGrid-1));
beta(x1,y1)=90/(NGrid-1);
[x1,y1]=find(beta<-90/(NGrid-1));
beta(x1,y1)=-90/(NGrid-1);
[peaks,index]=findpeaks(P_proposed);
[peaks,index2]=sort(peaks,'descend');
index=index(index2);
index=index(1:NSignals);
% figure
% plot(pFinal)
DOA_proposed=SearchGrid(index)+beta(index);
DOA_proposed=sort(DOA_proposed);
P_proposed=P_proposed/max(P_proposed);
% P_proposed=10*log(P_proposed)';


