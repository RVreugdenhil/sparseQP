
function [x,his] = SkPCA_lp(A,k)
for iter = 1:100,
    [x,flag_succ,his] = SkPCA_lp_one_run(A,k);
    if(1==flag_succ),break;end    
end
if(iter==100),
    warning('Max Iter!');
end
x=proj_l0(x,k);


function [x,flag_succ,his] = SkPCA_lp_one_run(A,k)
[x0,~] = eigs(A,1);
n = size(A,1);
lambda_min = 1e-5;
lambda_max = 1e5;
b = ones(n,1);
mid = 0.5*(lambda_min+lambda_max);
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
    [fun_value,x,his] = getCurrentK(A,b,x0,mid);
      fprintf('max:%f, min:%f, mid:%f\n',lambda_max,lambda_min,mid);
    if( abs(lambda_max-lambda_min) < 0.01*lambda_max),break;end
end
if(nnz(x)==k)
    flag_succ = 1;
else
    flag_succ = 0;
end

function [k_cur,x,his] = getCurrentK(A,b,x0,rho_approx)
[x,his] = SPCA_approx(A,b,x0,rho_approx,'Lp',1);
x(abs(x)<1e-4)=0;
k_cur = nnz(x);
