function [x,his]=NonlinearFractionalProgrammingGlobal2(x,A,b,c,D,e,f)
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

if(min(eig(A-alpha*D))<=eps)
    A
    alpha
    eig(A-alpha*D)
    ddddddd
end



[~,x]=HandleObj(alpha);
his = [];

function [L,U,delta,gamma,g,O] = computeBound(Q,p,w,R,c,v)
iR = inv(R);
iRsqrt = sqrtm(iR);
O = iRsqrt*Q*iRsqrt';
gamma = 2*v - norm(iRsqrt*c)^2;
if(gamma==0),gamma=eps;end
g = iRsqrt*p - iRsqrt*Q*iR*c;
delta = c'*iR*Q*iR*c - 2*c'*iR*p + 2*w;
Z = [O g/sqrt(gamma);g'/sqrt(gamma) delta/gamma];
Z = (Z+Z')/2; O = (O+O')/2;
L = min(eig(Z));
U = min(eig(O));


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


