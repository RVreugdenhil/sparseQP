function test_NonlinearFractionalProgramming
addpath('util','solver','data');
clc;clear all;close all;

mexC;

for iter = 1:111111,
    
    randn('seed',iter)
    rand('seed',iter)
    [A,b,c,D,e,f]=GenData(12);
    
    
    n = length(b);
    x = 1*randn(n,1);
    HandleObj = @(x)computeObj(x,A,b,c,D,e,f);
    max_iter = 100;
    tic;
    [x1] = QuadFractionalProgrammingCOO(x,A,b,c,D,e,f,n,max_iter);
    toc;

%     tic;
%     [x2] = QuadFractionalProgrammingCOO_optimized(x,A,b,c,D,e,f,n,max_iter);
%     toc
%     norm(x2-x1)

    [x2] = QuadFractionalProgrammingGlobal(x,A,b,c,D,e,f);
    
    f1 = HandleObj(x1);
    f2 = HandleObj(x2);
    
    f1-f2
    
    if(f1+0.0001<f2),
        f1
        f2
        iter
        ddd
    end
    
    % plot([1:length(his0)],his0,[1:length(his1)],his1)
end





function [u] = solve_frac_problem1(O,g,delta,gamma,lb,ub)
% min_{u} (0.5u'Ou+g'u+0.5*delta) / (0.5u'u+0.5gamma)
n = length(g);
HandleObj = @(alpha)computeObjAlpha2(alpha,O,g,delta,gamma);
[alpha] = fmincon(HandleObj,lb,[],[],[],[],lb,ub,[],optimset('TolX',1e-5,'TolF',1e-5,'Display','off'));
[~,u]=HandleObj(alpha);



function [fobj,u] = computeObjAlpha2(alpha,O,g,delta,gamma)
% min_{u} (0.5u'Ou+g'u+0.5*delta) / (0.5u'u+0.5gamma)
% (0.5u'Ou+g'u+0.5*delta) - alpha (0.5u'u+0.5gamma)
% Ou + g - alpha u = 0
% (O - alpha I)u = -g
n = length(g);
u = inv(O-alpha*eye(n))*(-g);
fobj = (0.5*u'*O*u+g'*u+0.5*delta) / (0.5*u'*u+0.5*gamma);

function [x,his]=NonlinearFractionalProgrammingDB1(x,A,b,c,D,e,f)
% [lb,ub] = computeBound(A,b,c,D,e,f);
n = length(b);
HandleObj = @(alpha)computeObjAlpha(alpha,A,b,c,D,e,f);
[lb,ub] = computeBound(A,b,c,D,e,f);
% show(lb,ub,HandleObj)

warning off;
% [alpha] = gss(HandleObj,lb,ub,1e-12,10000,0);
% [alpha] = fminsearch(HandleObj,0,optimset('TolX',1e-12,'TolF',1e-12));
% [alpha] = fminsearch(HandleObj,0,optimset('TolX',1e-12,'TolF',1e-12));
[alpha] = fmincon(HandleObj,lb,[],[],[],[],lb,ub,[],optimset('TolX',1e-15,'TolF',1e-15,'Display','off'));
% 'by Matlab Line Search'
% alpha

% if(min(eig(A-alpha*D))<=eps)
%     A
%     alpha
%     eig(A-alpha*D)
%     ddddddd
% end



[~,x]=HandleObj(alpha);
his = [];

function [fobj,x] = computeObjAlpha(alpha,A,b,c,D,e,f)
O = A-alpha*D;
x = inv(O)*(alpha*e-b);
% min(eig(O))
% x
% '8888888888888888'
fobj = (0.5*x'*A*x+b'*x+c) / (0.5*x'*D*x+x'*e+f);

% [x3,his] = NonlinearFractionalProgramming3(x,-A,-b,-c,D,e,f,HandleObj);
% [fff0,ggg0]=HandleObj(x0);
% [fff1,ggg1]=HandleObj(x1);
% [fff2,ggg2]=HandleObj(x2);
% [fff3,ggg2]=HandleObj(x3);

% function [x,his] = NonlinearFractionalProgrammingDB(x,A,b,c,D,e,f)
% % min_x (0.5 x'Ax + b'x + c) / (0.5 x'Dx + x'e + f)
% % HandleObj = @(x)computeObj(x,A,b,c,D,e,f);
% alpha = 1;
% opt = @(x,alpha) 0.5*x'*A*x+b'*x+c - alpha*(0.5*x'*D*x+x'*e+f);
% his = [];
% for iter = 1:50
%     % min_{x} 0.5x'Ax+b'x+c - alpha (0.5x'Dx+x'e+f)
%     % Ax+b - alpha (Dx+e) = 0
%     % (A-alpha D)x = alpha e-b
%       x =  pinv(A-alpha*D)*(alpha*e-b);
%       fobj = (0.5*x'*A*x+b'*x+c) / (0.5*x'*D*x+x'*e+f);
%       his = [his;fobj];
%       optimality = opt(x,alpha);
%       fprintf('iter:%d, F:%.1e, fobj:%e, alpha:%f\n',iter,optimality,fobj,alpha);
% %       0.5*x'*D*x+x'*e+f
%       alpha = (0.5*x'*A*x+b'*x+c) / (0.5*x'*D*x+x'*e+f);
% end


