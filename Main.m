clear all
clc
close all
load(['demodata.mat'])
X=full(ind2vec(LTra1')');
%% Training
Input.x=DTra1;Input.y=X;
[Output]=CLEAF(Input,'L');
MN=Output.Syst.ModelNumber
%% Testing
Input.x=DTes1;Input.Syst=Output.Syst;
[Output]=CLEAF(Input,'T');
label_est=Output.Ye;
[~,label_est]=max(label_est,[],2);
Acc=sum(sum(confusionmat(LTes1,label_est).*(eye(length(unique(LTes1))))))/length(LTes1)