function test_gspca
randn('seed',0);
rand('seed',0);
clc;clear all;close all;


n = 10;
A = randn(n,n);
A = -(A*A');
b = 0*rand(n,1);
c = 0;
D = randn(n);
D = D*D';
e = 0*randn(n,1);
f = 0;

HandleObj = @(x)computeObj(x,A,b,c,D,e,f);


while 1,
    x = randn(n,1);
    y = randn(n,1);
    lambda = rand(1);
    left = HandleObj(lambda*x+(1-lambda)*y);
    right = min(HandleObj(x),HandleObj(y));
    fprintf('%f %f\n',left-right,lambda);
end


function [x] = GS_PCA(A,B,k)
% min_x 0.5 x'Ax / 0.5 x'Bx, ||x||_0 = k
% We assume that B is strict positive definite
% We perform A <= A + theta B such that A + theta B is PSD
A = (A+A')/2;
B = (B+B')/2;
A_original = A;
B_original = B;

A = A + B * 1.001*abs(min(eig(A,B)));
n = size(A,1);
x = randn(n,1);
x = proj_l0(x,k);
num = 5;

for iter = 1:100,
    fobj = (0.5*x'*A_original*x)/(0.5*x'*B_original*x);
    fprintf('iter:%d, fobj:%f\n',iter,fobj);
    
    seq = randperm(n);
    iB = seq(1:num);
    iN = seq((num+1):n);
    k1 = nnz(x(iB));
    if(k1==0),continue;end
    
    sub_A = A(iB,iB);
    sub_b = A(iB,iN)*x(iN);
    sub_c = 0.5*x(iN)'*A(iN,iN)*x(iN);
    
    sub_D = B(iB,iB);
    sub_e = B(iB,iN)*x(iN);
    sub_f = 0.5*x(iN)'*B(iN,iN)*x(iN);
    
    x(iN) = solve_subproblem_k(sub_A,sub_b,sub_c,sub_D,sub_e,sub_f,k1);
end

nnz(x)

function [x] = solve_subproblem_k(A,b,c,D,e,f,k)
% min_x (0.5 x'Ax + b'x + c) / (0.5 x'Dx + x'e + f), s.t. ||x||_0=k
n = length(b);
[list] = combs([1:n],k);
nn = size(list,1);

xs = [];
fobjs = [];
for iter = 1:nn,
    cur = list(iter,:);
    x = zeros(n,1);
    x(cur) = solve_subproblem(x(cur),A(cur,cur),b(cur),c,D(cur,cur),e(cur),f);
    fobjs(iter) = (0.5*x'*A*x + b'*x + c) / (0.5*x'*D*x + x'*e + f);
    xs{iter} = x;
    if(norm(x)>10000),fobjs(iter)=1e100;end
end
[~,ind]=min(fobjs);
x  = xs{ind};

 

function [x,his] = solve_subproblem(x,A,b,c,D,e,f)
% min_x (0.5 x'Ax + b'x + c) / (0.5 x'Dx + x'e + f)
HandleObj = @(x)computeObj(x,A,b,c,D,e,f);
his=[];
for iter = 1:1000,
     [fobj_cur,grad_cur]=HandleObj(x);
     his = [his;fobj_cur];
      fprintf('iter: %d, fobj: %f, grad: %e\n',iter,fobj_cur,norm(grad_cur));
       if(iter==1)
            direction= - grad_cur ;
       else
%                 beta= mdot(grad_cur,grad_cur)/ mdot(grad_old,grad_old);                     % FR
           beta=  mdot(grad_cur,grad_cur-grad_old) / mdot(grad_old,grad_old);            % PRP
            beta(beta<0)=0;
            direction= - grad_cur + beta *   direction;
       end
       
      if(direction'*grad_cur>0)
           direction= - grad_cur ;
      end
        
    grad_old = grad_cur;
 
       % min_{s}  (0.5 a1 s^2 + a2 s + a3)  /  (0.5 a4 s^2 + a5 s + a6)
       a1 = direction'*A*direction;
       a2 = x'*A*direction + b'*direction;
       a3 = 0.5*x'*A*x + b'*x + c;
       a4 = direction'*D*direction;
       a5 = x'*D*direction + e'*direction;
       a6 = 0.5*x'*D*x + e'*x + f;
       
       % min_{s}  (0.5 a1 s^2 + a2 s + a3)  /  (0.5 a4 s^2 + a5 s + a6)
%       HandleLineSearch = @(s) (0.5*a1*s*s + a2*s + a3)  /  (0.5*a4*s*s + a5*s + a6);
%       s =  fminsearch(HandleLineSearch,0);

% grad_s =  [(a1s+a2)*(0.5 a4 s^2 + a5 s + a6)- (0.5 a1 s^2 + a2 s + a3)*(a4 s+a5) ]  /  (0.5 a4 s^2 + a5 s + a6)^2 =0

% (a1s+a2)*(0.5 a4 s^2 + a5 s + a6)- (0.5 a1 s^2 + a2 s + a3)*(a4 s+a5)=0

    c2 = 0.5*a1*a5-0.5*a2*a4;
    c1 = a1*a6-a3*a4;
    c0 = a2*a6-a3*a5;
    
    % To find roots to the following cubic equation where a3, a2, a1, and a0 are real
    % c2 x^2 + c1*x + a0 =0
    steps =roots([c2 c1 c0]);
    
    steps_candidates=[];
    for i=1:length(steps),
        if(isreal(steps(i)) && steps(i)> 0)
            steps_candidates=[steps_candidates;steps(i)];
        end;
    end;
    
    len=length(steps_candidates);
    if(len==0)
        % no good step, OK!
%         fprintf('!\n');
        break;
    elseif(len==1)
        % only one solution, the stepsize is optimal!
        s=steps_candidates(1);
    elseif(len>1)
        % multiple solutions, select the best one!
%         fprintf('*'); 
        find_min = inf;
        find_ind = -1;
        for i=1:length(steps_candidates),
            x_test = x + steps_candidates(i)*direction;
            find_fobj = HandleObj(x_test);
            if(find_fobj<find_min),
                find_min=find_fobj;
                find_ind = i;
            end
        end
        s = steps_candidates(find_ind);
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% line search end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
       x = x + s*direction;
       
end

%   plot(his)
%   ddd



function [fobj,grad] = computeObj(x,A,b,c,D,e,f)
% min_x (0.5x'Ax + b'x + c) / (0.5x'Dx + x'e + f)
up = (0.5*x'*A*x + b'*x + c);
down = (0.5*x'*D*x + x'*e + f);
up_g = A*x+b;
down_g = D*x+e;
fobj =  up / down;
grad = (up_g * down - up * down_g) / (down*down);


function [x] = subproblem(A,B,c,d,e)
% min_x 0.5x'Ax + c'x s.t. 0.5 x'Bx + d'x = e
% max_{theta} min_x 0.5x'Ax + c'x + theta * (0.5 x'Bx + d'x - e)
% A*x+c + theta *(B*x+d) = 0
% (A+theta *B) x = -c - theta * d
% x = inv(A+theta*B)*(-c-theta*d);
% HandleObj = @(theta) -0.5*x'Ax - c'x - theta * (0.5 x'Bx + d'x - e)
% HandleObj = @(theta) 0.5*x'(-A-thetaB)x - c'x - theta *d'x + theta *e
% HandleObj = @(theta) 0.5*x'(A+thetaB)x + theta *e
HandleObj = @(theta) 0.5* (c+theta*d)'*(inv(A+theta*B)*(c+theta*d)) + theta *e;
options.TolFun=1e-14;
options.TolX=1e-14;
theta = fminsearch(HandleObj,0,options);
x = -inv(A+theta*B)*(c+theta*d);

% 0.5*x'*B*x + d'*x - e
% 0.5*x'*A*x + c'*x



