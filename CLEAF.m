function [Output]=CLEAF(Input,Mode)
D_o=exp(-1/6);
kappa_o=1/2;
K_o=200;
Lambda_o=0.01;
rho_o=0.5;
if strcmp(Mode,'L')==1
    [Output.Ye,Output.Syst.ModelNumber,Output.Syst.prototypes,Output.Syst.center,Output.Syst.Support,Output.Syst.local_delta,Output.Syst.Global_mean,Output.Syst.Global_X,Output.Syst.A,Output.Syst.L]=EFStrain(Input.x,Input.y,D_o,kappa_o,K_o,Lambda_o,rho_o);
end
if strcmp(Mode,'T')==1
    [Output.Ye]=EFStest(Input.x,Input.Syst.ModelNumber,Input.Syst.prototypes,Input.Syst.center,Input.Syst.local_delta,Input.Syst.Global_mean,Input.Syst.Global_X,Input.Syst.A,Input.Syst.L,D_o,kappa_o,K_o);
end
Output.chunksize=K_o;
end
function [Ye,ModelNumber,prototypes,center,Support,local_delta,Global_mean,Global_X,A,L]=EFStrain(data0,Y0,threshold1,threshold2,chunksize,forgettingfactor,rho0)
CL=size(Y0,2);
[L,W]=size(data0);
Ye=zeros(L,CL);
omega=0;
seqck=[1:chunksize:L+1,L+1];
Lck=length(seqck)-1;
%%
datack0=data0(seqck(1):1:seqck(2)-1,:);
Yck0=Y0(seqck(1):1:seqck(2)-1,:);
Global_mean=sum(datack0,1);
Global_X=sum(datack0.^2,1);
Global_ymean=sum(Yck0,1);
Global_Y=sum(sum(Yck0.^2,2),1);
Global_Delta=(Global_X*(seqck(2)-1)-(Global_mean).^2)/(seqck(2)-1)^2;
Ystd=sum((Global_Y*(seqck(2)-1)-(Global_ymean).^2)./(seqck(2)-1)^2);
Xstd=sum(Global_Delta);
[prototypes,prototypeY,center,Local_X,centerY,Local_Y,Support,ModelNumber]=centerextraction(datack0,Xstd,Yck0,Ystd,rho0);
A=zeros(CL,W+1,ModelNumber);
C=repmat(eye(W+1)*omega,1,1,ModelNumber);
local_delta=abs(Local_X.*Support-center.^2)./Support.^2;
local_deltaY=abs(Local_Y.*Support-centerY.^2)./Support.^2;
[centerlambda,LocalDensity]=firingstrength(datack0,ModelNumber,prototypes,local_delta,Global_Delta,W);
sum_lambda=sum(LocalDensity,1)/(seqck(2)-1);
timeidx=ones(1,ModelNumber);
[centerlambda]=ActivatingRules(ModelNumber,centerlambda,LocalDensity,threshold2);
[A,C]=ConsequentUpdate(datack0,Yck0,centerlambda,A,C,ModelNumber);
for ii=2:1:Lck
    if seqck(ii+1)-1>seqck(ii)
        datack0=data0(seqck(ii):1:seqck(ii+1)-1,:);
        Yck0=Y0(seqck(ii):1:seqck(ii+1)-1,:);
        [centerlambda,LocalDensity]=firingstrength(datack0,ModelNumber,prototypes,local_delta,Global_Delta,W);
        Ye(seqck(ii):1:seqck(ii+1)-1,:)=OutputGeneration(datack0,A,centerlambda,LocalDensity,ModelNumber,CL,threshold2);
        Global_mean=Global_mean+sum(datack0,1);
        Global_X=Global_X+sum(datack0.^2,1);
        Global_ymean=Global_ymean+sum(Yck0,1);
        Global_Y=Global_Y+sum(Yck0.^2,1);
        Global_Delta=abs(Global_X*(seqck(ii+1)-1)-Global_mean.^2)/(seqck(ii+1)-1)^2;
        Global_DeltaY=abs(Global_Y*(seqck(ii+1)-1)-Global_ymean.^2)/(seqck(ii+1)-1)^2;
        Xstd=sum(Global_Delta);
        Ystd=sum(Global_DeltaY);
        [prototypes0,prototypeY0,center0,Local_X0,centerY0,Local_Y0,Support0,~]=centerextraction(datack0,Xstd,Yck0,Ystd,rho0);
        LocalDensity1=firingstrengthY(prototypes0,prototypeY0,ModelNumber,prototypes,prototypeY,local_delta,local_deltaY,Global_Delta,Global_DeltaY,W,rho0);
        [seq,idx1]=max(LocalDensity1,[],2);
        seqpp1=find(seq<threshold1);
        seqpp2=find(seq>=threshold1);
        for jj=seqpp2'
            idx=idx1(jj);
            Local_X(idx,:)=Local_X(idx,:)+Local_X0(jj,:);
            center(idx,:)=center(idx,:)+center0(jj,:);
            centerY(idx,:)=centerY(idx,:)+centerY0(jj,:);
            Local_Y(idx,:)=Local_Y(idx,:)+Local_Y0(jj,:);
            Support(idx)=Support(idx)+Support0(jj);
        end
        ModelNumber1=length(seqpp1);
        prototypes=[prototypes;prototypes0(seqpp1,:)];
        prototypeY=[prototypeY;prototypeY0(seqpp1,:)];
        Local_X=[Local_X;Local_X0(seqpp1,:)];
        Support=[Support;Support0(seqpp1)];
        center=[center;center0(seqpp1,:)];
        Local_Y=[Local_Y;Local_Y0(seqpp1,:)];
        centerY=[centerY;centerY0(seqpp1,:)];
        timeidx=[timeidx,ones(1,ModelNumber1)*ii];
        sum_lambda=[sum_lambda,zeros(1,ModelNumber1)];
        local_delta=abs(Local_X.*Support-center.^2)./Support.^2;
        local_deltaY=abs(Local_Y.*Support-centerY.^2)./Support.^2;
        A(:,:,ModelNumber+1:1:ModelNumber+ModelNumber1)=zeros(CL,W+1, ModelNumber1);
        C(:,:,ModelNumber+1:1:ModelNumber+ModelNumber1)=repmat(eye(W+1)*omega,1,1,ModelNumber1);
        ModelNumber=ModelNumber+ModelNumber1;
        [centerlambda,LocalDensity]=firingstrength(datack0,ModelNumber,prototypes,local_delta,Global_Delta,W);
        sum_lambda=(sum_lambda.*(ii-timeidx)+mean(LocalDensity,1))./(ii-timeidx+1);
        seq=find(sum_lambda>=forgettingfactor);
        ModelNumber0=length(seq);
        if ModelNumber0<ModelNumber
            center=center(seq,:);
            Local_X=Local_X(seq,:);
            centerY=centerY(seq,:);
            Local_Y=Local_Y(seq,:);
            timeidx=timeidx(seq);
            sum_lambda=sum_lambda(seq);
            Support=Support(seq);
            A=A(:,:,seq);
            C=C(:,:,seq);
            prototypes=prototypes(seq,:);
            prototypeY=prototypeY(seq,:);
            LocalDensity=LocalDensity(:,seq);
            ModelNumber=ModelNumber0;
            local_delta=local_delta(seq,:);
            local_deltaY=local_deltaY(seq,:);
        end
        [centerlambda]=ActivatingRules(ModelNumber,centerlambda,LocalDensity,threshold2);
        [A,C]=ConsequentUpdate(datack0,Yck0,centerlambda,A,C,ModelNumber);
    end
end
Global_mean=Global_mean/L;
Global_X=Global_X/L;
end
function [Ye]=EFStest(data1,ModelNumber,prototypes,center,local_delta,Global_mean,Global_X,A,L1,threshold1,threshold2,chunksize)
CL=size(A,1);
[L,W]=size(data1);
Ye=zeros(L,CL);
seqck=[1:chunksize:L+1,L+1];
Lck=length(seqck)-1;
Global_Delta1=abs(Global_X-Global_mean.^2);
for ii=1:1:Lck
    datack1=data1(seqck(ii):1:seqck(ii+1)-1,:);
    [centerlambda,LocalDensity]=firingstrength(datack1,ModelNumber,prototypes,local_delta,Global_Delta1,W);
    Ye(seqck(ii):1:seqck(ii+1)-1,:)=OutputGeneration(datack1,A,centerlambda,LocalDensity,ModelNumber,CL,threshold2);
end
end
function Ye=OutputGeneration(datain,A,centerlambda,LocalDensity,ModelNumber,CL,threshold2)
[L,W]=size(datain);
X=[ones(L,1),datain];
[centerlambda]=ActivatingRules(ModelNumber,centerlambda,LocalDensity,threshold2);
Ye=sum(pagemtimes(X,pagetranspose(A)).*reshape(centerlambda,L,1,ModelNumber),3);
end
%%
function [centerlambda]=ActivatingRules(ModelNumber,centerlambda,LocalDensity,threshold2)
centerlambda=LocalDensity;
for jj=1:1:size(centerlambda,1)
    LocalDensity1=LocalDensity(jj,:);
    [values,seq]=sort(LocalDensity1,'descend');
    values=sum(triu(repmat(values',1,ModelNumber)),1);
    a=find(values>=threshold2*sum(LocalDensity1));
    seq1=seq(1:1:a(1))';
    seq2=seq(a(1)+1:end)';
    centerlambda1=LocalDensity1;
    centerlambda1(seq1)=LocalDensity1(seq1)./sum(LocalDensity1(seq1));
    centerlambda1(seq2)=0;
    centerlambda(jj,:)=centerlambda1;
end
centerlambda(isnan(centerlambda))=0;
end
function [centerlambda,LocalDensity]=firingstrength(datain,ModelNumber,center,local_delta,Global_Delta,W)
[L,W]=size(datain);
Global_Delta1=sum((Global_Delta+local_delta)/2,2)';
datain1=pdist2(datain,center).^2;
LocalDensity=exp(-1*datain1./Global_Delta1);
LocalDensity(isnan(LocalDensity))=0;
centerlambda=LocalDensity./sum(LocalDensity,2);
centerlambda(isnan(LocalDensity))=0;
end
function [LocalDensity]=firingstrengthY(datain,dataout,ModelNumber,center,centerY,local_delta,local_deltaY,Global_Delta,Global_DeltaY,W,rho0)
LocalDensity=exp(-1*(rho0*pdist2(datain,center).^2./sum(Global_Delta+local_delta,2)'+(1-rho0)*pdist2(dataout,centerY).^2./sum(Global_DeltaY+local_deltaY,2)'));
end
function [prototype,prototypeY,center,LocalX,centerY,LocalY,support,LC]=centerextraction(datack0,Xstd,Yck0,Ystd,rho0)
[L,W]=size(datack0);
CL=size(Yck0,2);
if L==1
    prototype=datack0;
    prototypeY=Yck0;
    center=datack0;
    LocalX=datack0.^2;
    centerY=Yck0;
    LocalY=Yck0.^2;
    support=ones(size(datack0,1),1);
    LC=size(datack0,1);
else
    [datack0,ia]=unique(datack0,'rows');
    Yck0=Yck0(ia,:);
    distck=squareform(pdist(datack0).^2);
    distyck=squareform(pdist(Yck0).^2);
    distck1=exp(-(rho0*distck/Xstd+(1-rho0)*distyck/Ystd));
    distck_tempseq=ones(size(distck));
    distck_tempseq(distck>Xstd|distyck>Ystd)=0;
    distck_tempseq1=sum(distck_tempseq,2);
    distck_tempseq2=sum(distck1.*distck_tempseq,2)./distck_tempseq1;
    distck_tempseq3=distck_tempseq2.*distck_tempseq;
    seqpk=diag(distck_tempseq3)-max(distck_tempseq3,[],1)';
    seqpk=find(seqpk==0);
    prototype=datack0(seqpk,:);
    prototypeY=Yck0(seqpk,:);
    LC=length(seqpk);
    LocalX=zeros(LC,W);
    support=zeros(LC,1);
    center=zeros(LC,W);
    centerY=zeros(LC,CL);
    LocalY=zeros(LC,CL);
    distck=distck(:,seqpk);
    distyck=distyck(:,seqpk);
    [~,idx]=min((rho0*distck/Xstd+(1-rho0)*distyck/Ystd),[],2);
    for jj=1:1:LC
        tempseq=idx==jj;
        center(jj,:)=sum(datack0(tempseq,:),1);
        LocalX(jj,:)=sum(datack0(tempseq,:).^2,1);
        support(jj)=sum(tempseq);
        centerY(jj,:)=sum(Yck0(tempseq,:),1);
        LocalY(jj,:)=sum(Yck0(tempseq,:).^2,1);
    end
end
end
function [A,C]=ConsequentUpdate(datack0,Yck0,centerlambda,A,C,ModelNumber)
[L,W]=size(datack0);
X=[ones(L,1),datack0];
for jj=1:1:ModelNumber
   C(:,:,jj)=C(:,:,jj)+centerlambda(:,jj)'.*X'*X;
   A(:,:,jj)=A(:,:,jj)+centerlambda(:,jj)'.*((Yck0-X*A(:,:,jj)'))'*(X)*pinv(C(:,:,jj));
end
end