clear all
clc
close all
name{1}='NSLKDD';
name{2}='UNSWNB15';
name{3}='HIKARI';
for jj=1
for ii=1:1:1
    % ii
    load(['C:\Users\XwGu\Desktop\RRDEFNNADC\Data\' name{jj} '\data_' num2str(ii) '.mat'])
    % X=readtable('E:\Projects_published\DSTLProject\DSTL\Dataset\NSL-KDD\NSLKDD_trainsetoh.csv');
    % data=table2array(X(2:end,2:end));
    % 
    % LTra1=data(:,end)+1;
    % DTra1=data(:,1:1:end-1);
    % X=readtable('E:\Projects_published\DSTLProject\DSTL\Dataset\NSL-KDD\NSLKDD_testsetoh.csv');
    % data=table2array(X(2:end,2:end));
    % % seq=find(mean(abs(data),1)==0);
    % LTes1=data(:,end)+1;
    % DTes1=data(:,1:1:end-1);
    % seq=find(mean(abs([DTra1;DTes1]),1)==0);
    % DTra1=(DTra1-mean([DTra1;DTes1]))./std([DTra1;DTes1]);
    % DTes1=(DTes1-mean([DTra1;DTes1]))./std([DTra1;DTes1]);
    % DTra1(:,seq)=[];
    % DTes1(:,seq)=[];
    % seq=1:1:length(LTra1);%randperm(length(LTra1));
    X=full(ind2vec(LTra1')');
    % seq=1:1:100000;
    %% Training
    Input.x=DTra1;Input.y=X;
    tic
    [Output]=CLEF(Input,'L');
    MN(ii)=Output.Syst.ModelNumber;
    tt(ii)=toc;
    %% Testing
    Input.x=DTes1;Input.Syst=Output.Syst;
    [Output]=CLEF(Input,'T');
    label_est=Output.Ye;
    [~,label_est]=max(label_est,[],2);
    CM=confusionmat(LTes1,label_est);
    Acc(ii)=sum(sum(CM.*(eye(length(unique(LTes1))))))/length(LTes1);
    Ye1{ii}=label_est;
    TP=CM(1,1);
    TN=CM(2,2);
    FN=CM(1,2);
    FP=CM(2,1);
    Mcc(ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
end
% 
AA(jj)=1-mean(Acc);
[1-mean(Acc),std(Acc)]
[mean(Mcc),std(Mcc)]
[mean(MN),std(MN)]
[mean(tt),std(tt)]
save(['CLEF_' name{jj} '_results.mat'],'tt','Acc','Ye1','MN','Mcc')
end


% clear all
% clc
% close all
% name{1}='ELEC'
% name{2}='SKIN'
% name{3}='Weather'
% for jj=3
% for ii=1:1:10
%     % ii
%     load(['C:\Users\XwGu\Desktop\CLEF\Dataset\' name{jj} '\data_' num2str(ii) '.mat'])
%     %%
%     data=[DTra1;DTes1];
%     mind=min(data,[],1);
%     maxd=max(data,[],1);
%     DTra1=(DTra1-mind)./(maxd-mind);
%     DTes1=(DTes1-mind)./(maxd-mind);
%     %%
%     X=full(ind2vec(LTra1')');
%     % X(X==0)=-1;
%     %% Training
%     Input.x=DTra1;Input.y=X;
%     tic
%     [Output]=CLEF(Input,'L');
%     MN(ii)=Output.Syst.ModelNumber;
% 
%     tt(ii)=toc;
%     %% Testing
%     Input.x=DTes1;Input.Syst=Output.Syst;
%     [Output]=CLEF(Input,'T');
%     label_est=Output.Ye;
%     [~,label_est]=max(label_est,[],2);
%     Acc(ii)=sum(sum(confusionmat(LTes1,label_est).*(eye(length(unique(LTes1))))))/length(LTes1);
%     Ye1{ii}=label_est;
% end
% % 
% [1-mean(Acc),std(Acc)]
% [mean(MN),std(MN)]
% [mean(tt),std(tt)]
% % save(['CLEF_' name{jj} '_results.mat'],'tt','Acc','Ye1','MN')
% end
% % 
% clear all
% % clc
% close all
% name{1}='PR';
% name{2}='OR';
% name{3}='OD';
% name{4}='AB';
% name{5}='MG';
% name{6}='PW';
% name{7}='PB';
% name{8}='SP';
% name{9}='TE';
% name{10}='MA';
% for jj=[1:1:10]
% for ii=1:1:10
%     % ii
%     load(['D:\PapersPublished\SAFLE\EnsembleEFS\Data\' name{jj} '\data_' num2str(ii) '.mat'])
% 
%     X=full(ind2vec(LTra1')');
%     % X(X==0)=-1;
%     %% Training
%     Input.x=DTra1;Input.y=X;
%     tic
%     [Output]=CLEF(Input,'L');
%     MN(ii)=Output.Syst.ModelNumber;
% 
%     tt(ii)=toc;
%     %% Testing
%     Input.x=DTes1;Input.Syst=Output.Syst;
%     [Output]=CLEF(Input,'T');
%     label_est=Output.Ye;
%     [~,label_est]=max(label_est,[],2);
%     Acc(ii)=sum(sum(confusionmat(LTes1,label_est).*(eye(length(unique(LTes1))))))/length(LTes1);
%     Ye1{ii}=label_est;
% end
% % 
% AA(jj)=1-mean(Acc);
% [1-mean(Acc),std(Acc)]
% [mean(MN),std(MN)]
% [mean(tt),std(tt)]
% save(['CLEF_' name{jj} '_results.mat'],'tt','Acc','Ye1','MN')
% end