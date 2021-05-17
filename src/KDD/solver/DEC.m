function [x,his,tt,time_spend] = DEC(A,C,k,workset,solver)
% This program solves the following optimization problem using decomposition method:
% min_x x'Ax / x'Cx, ||x||_0 = k

% Input parameters:
% workset: [i,j]. selection i coordinate using the random strategy and j coordinates using the swapping strategy
% solver: 0(bisection method) or 1 (coordinate descent method)


tbegin = cputime;
n = size(A,1);
if(isempty(C)),C = eye(n);end
fprintf('\n****************************************************************************************\n');
fprintf('dim:%d, sparsity:%d, working set: R%d-G%d, solver: %d\n',n,k,workset(1),workset(2),solver);
warning off;
A = (A+A')/2;
C = (C+C')/2;
% x = proj_l0(x,k);

x = randn(n,1);
% x = x/norm(x);
% [val,idx]=sort(diag(A)./diag(C), 'ascend');
% x = zeros(n,1);
% x(idx(1:k)) = randn(k,1);
% x = ones(n,1);
% x(1:k)=1;
%  [x,~]=eigs(A,C,1);
x = proj_l0(x,k);

his = [];
HandleObj = @(x) ComputeMainObj(x,A,C,k);
last_k = 50;
changes = ones(last_k,1);
fobj_old = HandleObj(x);

accuracy = 1e-5;
theta = 1e-5;

selection = 0;
if(selection==0)
worksetmethod = @(x)find_work_set_most_violate_pair(x,A,C,workset(2));
elseif(selection==1)
    worksetmethod = @(x)find_work_set_most_violate_coordinate(x,A,C,workset(2));
end

tt = [];
t1 = clock();
for iter = 1:10000
    iB1 = worksetmethod(x);
%       iB2 = find_work_set_random(x,workset(1));
    iB2 = randperm(n,workset(1))';

    iB = [iB1;iB2];
    iB = unique(iB);
    iN = setdiff([1:n],iB);
    
    xB = x(iB);
    xN = x(iN);
    sparse_ws = nnz(xB);
    size_ws = length(iB);
    
    sub_A = A(iB,iB);
    sub_b = A(iB,iN)*xN;
    sub_c = 0.5*xN'*A(iN,iN)*xN;
    
    sub_D = C(iB,iB);
    sub_e = C(iB,iN)*xN;
    sub_f = 0.5*xN'*C(iN,iN)*xN;
    
    %                       0.5*xB'*sub_A*xB + xB'*sub_b + sub_c
    % min_{||xB||_0=k1}  ------------------------------------------
    %                       0.5*xB'*sub_D*xB + xB'*sub_e + sub_f
    
    x(iB) = solve_subproblem_k(xB,sub_A+theta*eye(size_ws),sub_b-theta*xB,sub_c+0.5*theta*norm(xB,2)^2,sub_D,sub_e,sub_f,sparse_ws,solver);
    fobj = HandleObj(x);
    his(iter) = fobj;
    tt = [tt;etime(clock(),t1)];
    rel_change = abs((fobj - fobj_old)/max(1,fobj_old));
    changes = [changes(2:end);rel_change];
    fobj_old = fobj;
    mchanges = mean(changes);
     if(~mod(iter,20))
        fprintf('iter: %d, fobj: %f, changes: %f, nnzx: %d, nnz xB: %d, size working set: %d\n',iter,fobj,mchanges,nnz(x),nnz(xB),length(iB));
    end
    if(mchanges<accuracy),break;end
    
end

x = proj_l0(x,k);
x = reoptimizer(x,A,C);
x = proj_l0(x,k);
tend = cputime;
  time_spend = tend - tbegin;


function [fobj] = ComputeMainObj(x,A,C,k)
% min_{x} (x'*A*x) / (x'*C*x)
x = proj_l0(x,k);
fobj = (x'*A*x) / (x'*C*x);

function [x_best] = solve_subproblem_k(x0,A,b,c,D,e,f,k,solver)
% min_x (0.5 x'Ax + b'x + c) / (0.5 x'Dx + x'e + f), s.t. ||x||_0 = k
n = length(b); 
list = combs([1:n],k); 
% list = getlist_random(n,k); 

list_n = size(list,1);
HandleObj = @(x)(0.5*x'*A*x + b'*x + c) / (0.5*x'*D*x + x'*e + f);
x_best = x0;
f_best = HandleObj(x0);
his = [];
solver = 1;

if(solver==0)
HandleObjSolver = @(x,Q,p,w,R,c,v)QuadFractionalProgrammingGlobal(x,Q,p,w,R,c,v);
else
 HandleObjSolver = @(x,Q,p,w,R,c,v)QuadFractionalProgrammingCOO_optimized(x,Q,p,w,R,c,v,size(x,1),10);
% HandleObjSolver = @(x,Q,p,w,R,c,v)QuadFractionalProgrammingCOO(x,Q,p,w,R,c,v,size(x,1),30);
end

for iter = 1:list_n
    cur = list(iter,:);
    x = zeros(n,1);
    [x_sub] = HandleObjSolver (x0(cur),A(cur,cur),b(cur),c,D(cur,cur),e(cur),f);

    x(cur) = x_sub;
    f_curr = HandleObj(x);
    if(f_curr<f_best),
        f_best=f_curr;
        x_best = x;
    end
end



