function [x,his,tt,time_spend] = TruncatedRwyleighFlow(A,B,k)
% max_x x'Ax/x'Bx
max_iter = 500;
 
%   t=eigs(A,B,1)*1.2; A = A+t*B;
    t = 0;

n = size(A,1);

x = randn(n,1);    x  = x/norm(x);
x  = proj_l0(x,k);
eta = 0.5/eigs(B,1,'lm');
% eta = 1/eigs(B,1,'lm');

xs = [];
his = [];
tt = [];
t1 = clock();
tbegin = cputime;
last_fobj = -1e100;
for iter =1:max_iter
    
    Bx = B*x;
    Ax = A*x;
    xBx = x'*Bx;
    if(xBx<1e-3),break;end
    xAx = x'*Ax;
    rho = xAx/(xBx) -t ;
    curr_fobj = rho;
    if (abs((curr_fobj-last_fobj)/curr_fobj) < 1e-4)
        break;
    end
    last_fobj = curr_fobj;
%     xs{iter} = x;
    his(iter) = -rho ;
     tt = [tt;etime(clock(),t1)];
    if(~mod(iter,30))
%      fprintf('iter:%d, fobj:%f\n',iter,rho);
    end
%     C  = I + eta/rho * (A-rho*B);
%     Cx = C*x;
      if(abs(rho)<1e-8),break;end
    Cx  = x + eta/rho * (Ax-rho*Bx);
    normCx = norm(Cx);
    if(normCx==0), x = rndx(x); break;end
    x  = Cx/normCx;
    if(norm(x)==0), x = rndx(x);break; end
    x  = proj_l0(x,k);
    x  = x/norm(x);
 
end

if(~isLegal(x))
    x = randn(size(A,1),1);
end



x=proj_l0(x,k);
tend = cputime;
  time_spend = tend - tbegin;
  
% [~,index] = min(his);
% x = xs{index};

function [x] = rndx(x)
 x = x + randn(size(x))*1e-5;
function [legal] = isLegal(v)
legal = sum(any(imag(v(:))))==0 & sum(isnan(v(:)))==0 & sum(isinf(v(:)))==0;