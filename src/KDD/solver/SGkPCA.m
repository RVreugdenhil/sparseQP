function [x,his,time_spend] = SGkPCA(A,B,k,method)
% This program sovles the following optimization problem:
% max_x x'Ax/x'Bx, s.t. ||x||_0<=k

max_iter = 1;
for iter = 1:max_iter,
    fprintf('%d ',iter);
    [x,flag_succ,his,time_spend] = SGkPCA_one_run(A,B,k,method);
    if(1==flag_succ),break;end    
end

if(~isLegal(x))
    x = randn(size(A,1),1);
end

if(iter==max_iter),
    warning('Max Iter!');
end
x = proj_l0(x,k);



function [x,flag_succ,his,time_spend] = SGkPCA_one_run(A,B,k,method)

tbegin = cputime;
n = size(A,1);
x0 = randn(n,1);
x0 = x0/sqrt(x0'*B*x0);

lambda_min = 0.0001;
lambda_max = 10000;

mid = 0.5*(lambda_min+lambda_max);
tttt = 1;
[fun_value,x,his] = getCurrentK(A,B,x0,mid,method);

for iter = 1:30,
    if fun_value > k
        lambda_min = mid;
    elseif (fun_value<k)
        lambda_max = mid;
    elseif (fun_value==k)
%         output x;
        break;
    end
    mid = 0.5*(lambda_min+lambda_max);
    tttt = tttt + 1;
    [fun_value,x,his] = getCurrentK(A,B,x0,mid,method);
    fprintf('max:%f, min:%f, mid:%f\n',lambda_max,lambda_min,mid);
    if( abs(lambda_max-lambda_min) < 0.01*lambda_max),break;end
end
tend = cputime;

if(nnz(x)==k)
    flag_succ = 1;
  
else
    flag_succ = 0;
end
  time_spend = tend - tbegin;
  
function [k_cur,x,his] = getCurrentK(A,B,x0,rho,method)

if(method==1)
[x,his] = SGEP_solve(A,B,x0,rho,'Lp',0.1);
elseif(method==2)
[x,his] = SGEP_solve(A,B,x0,rho,'log',0.1);
elseif(method==3)
[x,his] = SGEP_solve(A,B,x0,rho,'exp',0.1);
end
x(abs(x)<1e-4)=0;
k_cur = nnz(x);
