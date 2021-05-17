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
    [A,B] = getFisher_AC(dat_x,dat_y);
    opt_A = -A; opt_C = B;
    
    HandleObj = @(x)ComputeMainObj(x,opt_A,opt_C,k);
    [~,~,tt1] = SGkPCA(A,B,k,1);
    [~,~,tt2] = SGkPCA(A,B,k,2);
    [~,~,tt3] = SGkPCA(A,B,k,3);
    [~,~,~,tt4] = TruncatedRwyleighFlow(A,B,k);
    [~,~,~,tt5] = DEC(-A,B,k,[6;0],0);
    [~,~,~,tt6] = DEC(-A,B,k,[6;0],1);
    
    tts1 = [tts1;tt1]; tts2 = [tts2;tt2]; tts3 = [tts3;tt3]; 
    tts4 = [tts4;tt4]; tts5 = [tts5;tt5]; tts6 = [tts6;tt6];
    
end

result{1}=tts1;
result{2}=tts2;
result{3}=tts3;
result{4}=tts4;
result{5}=tts5;
result{6}=tts6;



save(mfilename,'result');
