clc;clear all;close all;
addpath('util','solver','data');
randn('seed',0);
rand('seed',0);

ks = [4:4:40];

tts1 = []; tts2 = []; tts3 = [];
tts4 = []; tts5 = []; tts6 = [];
tts7 = []; tts8 = []; tts9 = [];

for iii=1:length(ks),
    
    k = ks(iii);
    [dat_x,dat_y] = SelectData(2);
    
Cov  = mycov(dat_x);
opt_A  = -Cov;
opt_C = eye(size(Cov,1));
HandleObj = @(x)ComputeMainObj(x,opt_A,opt_C,k);



[~,~,tt1] = SGkPCA(Cov,eye(size(Cov,1)),k,3);
[~,~,tt2] = SGkPCA(Cov,eye(size(Cov,1)),k,2);
[~,~,tt3] = SGkPCA(Cov,eye(size(Cov,1)),k,1);
[~,~,tt4] = SkPCA_l0(Cov,k);
[~,~,~,tt5] = TruncatedRwyleighFlow(Cov,eye(size(Cov,1)),k);
[~,~,~,tt6] = TPower(Cov,k);
[~,~,~,tt7] = cwPCA(Cov,k,'type','GCW');
[~,~,~,tt8] = DEC(opt_A,opt_C,k,[6;0],0);
[~,~,~,tt9] = DEC(opt_A,opt_C,k,[6;0],1);

tts1 = [tts1;tt1];
tts2 = [tts2;tt2];
tts3 = [tts3;tt3];
tts4 = [tts4;tt4];
tts5 = [tts5;tt5];
tts6 = [tts6;tt6];
tts7 = [tts7;tt7];
tts8 = [tts8;tt8];
tts9 = [tts9;tt9];

end

result{1}=tts1;
result{2}=tts2;
result{3}=tts3;
result{4}=tts4;
result{5}=tts5;
result{6}=tts6;
result{7}=tts7;
result{8}=tts8;
result{9}=tts9;
    

save(mfilename,'result');
