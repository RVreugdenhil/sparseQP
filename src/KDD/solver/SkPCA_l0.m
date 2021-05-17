function [x,his,time_spend] = SkPCA_l0(A,k)
for iter = 1:100,
    [x,flag_succ,his,time_spend] = SkPCA_l0_one_run(A,k);
    if(1==flag_succ),break;end    
end
if(iter==100),
    warning('Max Iter!');
end
if(~isLegal(x))
    x = randn(size(A,1),1);
end
x=proj_l0(x,k);


function [x,flag_succ,his,time_spend] = SkPCA_l0_one_run(A,k)

tbegin = cputime;
% [x0,~] = eigs(A,1);
n = size(A,1);
x0 = randn(n,1);
x0 = x0/sqrt(x0'*x0);



lambda_min = 1e-5;
lambda_max = 1e5;
b = ones(n,1);
mid = 0.5*(lambda_min+lambda_max);
tttt = 1;
[fun_value,x,his] = getCurrentK(A,b,x0,mid);

for iter = 1:100,
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
    [fun_value,x,his] = getCurrentK(A,b,x0,mid);
      fprintf('max:%f, min:%f, mid:%f\n',lambda_max,lambda_min,mid);
    if( abs(lambda_max-lambda_min) < 0.01*lambda_max),break;end
end

tend = cputime;


if(nnz(x)==k)
    flag_succ = 1;
      
else
    flag_succ = 0;
end
x = full(x);

  time_spend = tend - tbegin;
function [k_cur,x,his] = getCurrentK(A,b,x0,rho_approx)
[x,his] = SPCA_L0(A,b,x0,rho_approx);
x(abs(x)<1e-4)=0;
k_cur = nnz(x);
