clear all
clc
close all
% name='LR'
% name='Hyperplan'
% name='SUSY'
for ii=1:1:1
    ii
    load (['D:\HPC\Kent\SAFE\Data2\CoverType\data_' num2str(ii) '.mat'])

    X=full(ind2vec(LTra1')');
    X(X==0)=-1;
    %% Training
    Input.x=DTra1;Input.y=X;
    tic
    [Output]=CLEF(Input,'L');
    tt(ii)=toc;
    %% Testing
    Input.x=DTes1;Input.Syst=Output.Syst;
    [Output]=CLEF(Input,'T');
    label_est=Output.Ye;
    [~,label_est]=max(label_est,[],2);
    Acc(ii)=sum(sum(confusionmat(LTes1,label_est).*(eye(length(unique(LTes1))))))/length(LTes1);
    Ye1{ii}=label_est;
end
mean(tt)
[mean(Acc),std(Acc)]
save(['CLEF_' name '_results.mat'],'tt','Acc','Ye1')