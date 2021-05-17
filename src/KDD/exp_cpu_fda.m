clc;clear all;close all;
addpath('util','solver','data');
randn('seed',0);
rand('seed',0);

ds = [7;10];
ks = [15;15];



result = [];
for iii=1:length(ds),
    
    k = ks(iii);
    [dat_x,dat_y] = SelectData(ds(iii));
    
    [A,B] = getFisher_AC(dat_x,dat_y);
    opt_A = -A;
    opt_C = B;
    
    HandleObj = @(x)ComputeMainObj(x,opt_A,B,k);
    Alg1 = @(A,B,k)TruncatedRwyleighFlow(A,B,k);
    Alg2 = @(opt_A,opt_C,k)DEC(opt_A,opt_C,k,[6;0],0);
    Alg3 = @(opt_A,opt_C,k)DEC(opt_A,opt_C,k,[10;0],0);
    Alg4 = @(opt_A,opt_C,k)DEC(opt_A,opt_C,k,[0;6],0);
    Alg5 = @(opt_A,opt_C,k)DEC(opt_A,opt_C,k,[0;10],0);
    Alg6 = @(opt_A,opt_C,k)DEC(opt_A,opt_C,k,[4;4],0);
    Alg7 = @(opt_A,opt_C,k)DEC(opt_A,opt_C,k,[6;6],0);
    Alg8 = @(opt_A,opt_C,k)DEC(opt_A,opt_C,k,[8;8],0);
    
    [x1,his1,tt1] = Alg1(A,B,k);
    [x2,his2,tt2] = Alg2(opt_A,opt_C,k);
    [x3,his3,tt3] = Alg3(opt_A,opt_C,k);
    [x4,his4,tt4] = Alg4(opt_A,opt_C,k);
    [x5,his5,tt5] = Alg5(opt_A,opt_C,k);
    [x6,his6,tt6] = Alg6(opt_A,opt_C,k);
    [x7,his7,tt7] = Alg7(opt_A,opt_C,k);
    [x8,his8,tt8] = Alg8(opt_A,opt_C,k);
    
    %% color
    pcolor=[];
    
    pcolor.red      = [255,66,14]/256;    pcolor.red2     = [255,113,31]/256;   pcolor.red3     = [236,42,45]/256;
    pcolor.blue     = [0,149,201]/256;    pcolor.blue2    = [114,214,238]/256;  pcolor.blue3    = [0,183,240]/256;  pcolor.blue4    = [0,192,192]/256;   pcolor.blue5    = [65,146,228]/256;
    pcolor.green    = [0,170,77]/256;     pcolor.green2   = [81,157,28]/256;    pcolor.green3   = [36,178,76]/256;  pcolor.green4   = [192,192,0]/256;   pcolor.green5   = [76,181,60]/256;
    pcolor.purple   = [143,76,178]/256;   pcolor.purple2  = [125,66,210]/256; pcolor.purple3  = [0,0,255]/256;
    pcolor.crimson  = [192,52,148]/256;   pcolor.crimson2 = [245,147,202]/256; pcolor.crimson3  = [212,22,118]/256;
    pcolor.orange   = [254,181,89]/256;   pcolor.orange2  = [255,129,0]/256;    pcolor.orange3  = [219,130,1]/256;
    pcolor.gray     = [128,128,127]/256;
    pcolor.yellow   = [255,204,51]/256;   pcolor.yellow2  = [248,235,46]/256;
    pcolor.pink     = [255,130,160]/256; pcolor.pink2 = [255,0,255]/255;
    pcolor.def      = [256,256,230]/256;
    pcolor.darkred = [191/256 79/256 75/256];
    colors  = [];
    
    figure('color','w')
    %     plot([1:length(his1)],his1,'--*','LineWidth',5,'MarkerSize',22,'color', pcolor.crimson ); hold on;
    %     plot([1:length(his2)],his2,'-o','LineWidth',5,'MarkerSize',22,'color',  pcolor.blue ); hold on;
    %     plot([1:length(his3)],his3,'-.x','LineWidth',5,'MarkerSize',22,'color', pcolor.gray); hold on;
    %     plot([1:length(his4)],his4,'--^','LineWidth',5,'MarkerSize',22,'color', pcolor.darkred,'MarkerEdgeColor',pcolor.green3,'MarkerFaceColor', pcolor.yellow); hold on;
    %     plot([1:length(his5)],his5,'LineWidth',5,'MarkerSize',22,'color', pcolor.yellow); hold on;
    %     plot([1:length(his6)],his6,'LineWidth',5,'MarkerSize',22,'color', pcolor.purple,'MarkerEdgeColor',pcolor.pink,'MarkerFaceColor', pcolor.blue4); hold on;
    %     plot([1:length(his7)],his7,'-s','LineWidth',5,'MarkerSize',22,'color', pcolor.orange3); hold on;
    
    
    loglog(tt1,his1,'-','LineWidth',5,'MarkerSize',22,'color', pcolor.crimson3 ); hold on;
    loglog(tt2,his2,'-','LineWidth',5,'MarkerSize',22,'color',  pcolor.blue ); hold on;
    loglog(tt3,his3,'-','LineWidth',5,'MarkerSize',22,'color', pcolor.gray); hold on;
    loglog(tt4,his4,'-','LineWidth',5,'MarkerSize',22,'color', pcolor.darkred); hold on;
    loglog(tt5,his5,'-','LineWidth',5,'MarkerSize',22,'color', pcolor.yellow); hold on;
    loglog(tt6,his6,'-','LineWidth',5,'MarkerSize',22,'color', pcolor.purple); hold on;
    loglog(tt7,his7,'-','LineWidth',5,'MarkerSize',22,'color', pcolor.orange3); hold on;
    loglog(tt8,his8,'-','LineWidth',5,'MarkerSize',22,'color', pcolor.green); hold on;
    % loglog(tt9,his9,'-','LineWidth',5,'MarkerSize',22,'color', pcolor.red); hold on;
    % loglog(tt10,his10,'-','LineWidth',5,'MarkerSize',22,'color', pcolor.pink); hold on;
    
    hleg=legend('Power','CWM','TRF','R6S0','R10S0','R0S6','R0S10','R2S8','R8S2','R8S10') ;
    
    fprintf('%.5f\n',HandleObj(x1));
    fprintf('%.5f\n',HandleObj(x2));
    fprintf('%.5f\n',HandleObj(x3));
    fprintf('%.5f\n',HandleObj(x4));
    fprintf('%.5f\n',HandleObj(x5));
    fprintf('%.5f\n',HandleObj(x6));
    fprintf('%.5f\n',HandleObj(x7));
    fprintf('%.5f\n',HandleObj(x8));
    % fprintf('%.5f\n',HandleObj(x9));
    
    
    
    One = [];
    One.x1 = x1;One.x2 = x2;One.x3 = x3;One.x4 = x4;One.x5 = x5;One.x6 = x6;One.x7 = x7;One.x8 = x8;
    One.his1 = his1;One.his2 = his2;One.his3 = his3;One.his4 = his4;One.his5 = his5;One.his6 = his6;One.his7 = his7;One.his8 = his8;
    One.tt1 = tt1;One.tt2 = tt2;One.tt3 = tt3;One.tt4 = tt4;One.tt5 = tt5;One.tt6 = tt6;One.tt7 = tt7;One.tt8 = tt8;
    
    result{iii} = One;
    pause(1)
end


save(mfilename,'result');
