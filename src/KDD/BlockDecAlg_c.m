function [x,his] = BlockDecAlg_c(x,A,b,k,workingset,flaginit)
% min_x 0.5||Ax-b||_2^2, s.t. ||x||_0 = k
% x: n x 1
% A: m x n
% b: m x 1
% min_x 0.5x'Qx + p'x, s.t. ||x||_0 = k

if(flaginit)
    [x] = proximal_gradient_l0c(x,A,b,k);
end

% accuracy
max_iter = 10000000;
last_k = 200;
stop_accuracy = 1e-5;

% x = proximal_gradient_l0c(x,Q,p,k);
const = 0.5*b'*b;
n = size(A,2);
x = proj_l0(x,k);

his = [];
theta = 1e-5;

changes = ones(last_k,1);

Q = A'*A;
p = -A'*b;
grad = Q*x+p;

HandleObj = @(x)0.5*x'*Q*x + p'*x;%0.5*norm(A*x-b,'fro')^2;
% fobj_old = HandleObj(x)+const;
fobj_old = (grad+p)'*x/2 + const;

% if(sum(workingset)>20),
%     subproblem_solver = @(a,b,c,d)quad_l0c_global_cplex_matlab(a,b,c,d);
% else
%     subproblem_solver = @(a,b,c,d)quad_l0c_global(a,b,c,d);
% end
n_ws_real=sum(workingset);
list = GenC(n_ws_real);
% fprintf('Block Decomposition Algorithm Start..\n'); 
t1 = clock();
for iter = 1:max_iter
    
    B1 = randperm(n,workingset(1));
    B2 = find_index_approximate_matrix(x,Q,p,workingset(2),B1);
    B = [B1(:);B2(:)];
    n_ws = length(B);
    k_B = nnz(x(B));
    if(n_ws~=n_ws_real),error('number of the working set is not correct!');end
    if(k_B==0),continue;end
    %     [fobj] = HandleObj(x)+const;
    [fobj] = (grad+p)'*x/2+const;
    
    rel_change = abs((fobj - fobj_old)/max(1,fobj_old));
    changes = [changes(2:end);rel_change];
    fobj_old = fobj;
 
    if(~mod(iter,50))
        fprintf('iter:%d, fobj:%f, sparsity:%d, work_set:%d, nnz_B: %d\n',iter,fobj,nnz(x),n_ws,k_B); his = [his;fobj];
   end
    
    H1 = Q(B,B) + eye(n_ws)*theta;
    x_B_old = x(B);
%     x_B_new = quad_l0c_global(H1,grad(B)-H1*x_B_old,k_B);
      x_B_new = quad_l0c_global(H1,grad(B)-H1*x_B_old,k_B,list{k_B});
    
    x(B) = x_B_new;
    % reconstruct the gradient
    grad = grad + Q(:,B)*(x_B_new-x_B_old);
    if(mean(changes)<stop_accuracy),break;end
end
t2 = clock();
tt = etime(t2,t1);
%  fprintf('time spent(s): %d\n',tt);

function [I] = find_index_slow(x,Q,p,num)
n = length(x);
I1 = find(x==0);
I2 = find(x~=0);

aa = []; bb = []; fobj=[];
for i=1:length(I1)
    for j=1:length(I2)
        first = I1(i);
        second = I2(j);
        aa = [aa;first];
        bb = [bb;second];
        x1 = x;
        x1(second) = 0;
        
        % minimize of the index for first
        B = first;
        N = setdiff([1:n],B);
        %  \min_{x1(B)}   0.5*(x1(N))'*Q(N,N)*x1(N) +     0.5*(x1(B))'*Q(B,B)*x1(B)  +     (x1(B))'*Q(B,N)*x1(N) + x1(N)'*p(N) + x1(B)'*p(B)
        %  \min_{x1(B)}   0.5*(x1(B))'*Q(B,B)*x1(B)  +  (x1(B))'*Q(B,N)*x1(N) + x1(B)'*p(B)
        x1(B) = -(Q(B,N)*x1(N)+p(B)) / Q(B,B);
        fobj = [fobj;0.5*x1'*Q*x1+p'*x1];
    end
end

[~,ind]=sort(fobj,'ascend');
aa = aa(ind); bb = bb(ind);
aa=find_unique_k(aa,num/2);
bb=find_unique_k(bb,num/2);
I = unique([bb;aa]);

function [index] = find_index(x,Q,p,num)
n = length(x);
I1 = find(x==0);
I2 = setdiff([1:n],I1);

Qx = Q*x;
xQx = x'*Qx;
ptx = p'*x;
aa = []; bb = []; fobj=[];
for i=1:length(I1)
    for j=1:length(I2)
        first = I1(i);
        second = I2(j);
        cof1 = x(first);
        cof2 = x(second);
        aa = [aa;first];
        bb = [bb;second];
        
        x(second) = 0;
        diff = cof2-x(second);
        Qx = Qx - Q(:,second)*(cof2-x(second));
        xQx = xQx - (2*Qx(second)*diff + diff*diff * Q(second,second));
        ptx = ptx - p(second) * diff;
        
        % minimize of the index for first
        B = first;
        %         N = setdiff([1:n],B);
        %  \min_{x1(B)}   0.5*(x1(N))'*Q(N,N)*x1(N) +     0.5*(x1(B))'*Q(B,B)*x1(B)  +     (x1(B))'*Q(B,N)*x1(N) + x1(N)'*p(N) + x1(B)'*p(B)
        %  \min_{x1(B)}   0.5*(x1(B))'*Q(B,B)*x1(B)  +  (x1(B))'*Q(B,N)*x1(N) + x1(B)'*p(B)
        
        
        %         x(B) = -(Q(B,N)*x(N)+p(B)) / Q(B,B);
        x(B) = -(Qx(B) - Q(B,B)*x(B)+p(B)) / Q(B,B);
        Qx = Qx - Q(:,B)*(cof1-x(B));
        diff = cof1-x(B);
        xQx = xQx - (2*Qx(B)*diff + diff*diff * Q(B,B));
        ptx = ptx - p(B) * diff;
        
        fobj = [fobj;0.5*xQx+ptx];
        
        % recover the original x
        ooo = x(second);
        x(second) = cof2;
        diff = ooo-cof2;
        Qx = Qx - Q(:,second)*diff;
        xQx = xQx - (2*Qx(second)*diff + diff*diff * Q(second,second));
        ptx = ptx - p(second) * diff;
        
        ooo = x(first);
        x(first)  = cof1;
        diff = ooo-cof1;
        Qx = Qx - Q(:,first)*diff;
        xQx = xQx - (2*Qx(first)*diff + diff*diff * Q(first,first));
        ptx = ptx - p(first) * diff;
        
        
    end
end


[~,ind]=sort(fobj,'ascend');
aa = aa(ind); bb = bb(ind);
aa=find_unique_k(aa,num/2);
bb=find_unique_k(bb,num/2);
index = unique([bb;aa]);



function [index] = find_index_approximate_matrix(X,Q,P,num,deleteSet)
% find the index based on the zero-ordered information
% min_X 0.5 ||AX-B||_2^2, s.t. ||X||_0 <=k
% min_X 0.5 trace(X'A'AX) - <AX,B>, s.t. ||X||_0 <=k
% min_X 0.5X'QX+P'X
% Q = A'A
% P = -A'B
% A: m x d
% X: d x r
% B: m x r

[Z] = find(X==0); [S] = find(X~=0);
Z = setdiff(Z,deleteSet(:)); 
S = setdiff(S,deleteSet(:)); 
diagQ = diag(Q);
O = Q*X;
[d,r] = size(X);

% E = zeros(d,r);
% i = 3;j = 3;
% E(i,j)=1;
% Q
% trace(E'*Q*E)
% % Q*E
% Q(3,3)
%
% ddd

%  fobjs1 = zeros(length(Z),1);
% for k=1:length(Z)
%     % we change zero to nonzero
%     i = Z(k);
%     [row,col]=co2xy(i,d);
%     % min_{alpha} 0.5 (X+alpha Eij)' Q (X+alpha Eij) + P' (X+alpha Eij)
%     % min_{alpha} 0.5 (alpha Eij)' Q (alpha Eij) + X' Q (alpha Eij) + P' (alpha Eij)
%     % grad_alpha = Q(i,i)alpha + O(ij) + P(ij)=0
%     alpha = -(P(i)+O(i)) / diagQ(row);
%     fobjs1(k) = 0.5*alpha*alpha*diagQ(row) + alpha*O(i) + alpha*P(i);
% end

rows = mod(Z-1,d)+1;
alpha = -(P(Z)+O(Z))./diagQ(rows);
fobjs1 = 0.5*alpha.*alpha.*diagQ(rows) + alpha.*O(Z) + alpha.*P(Z);


% for k=1:length(S)
%     % we change nonzero to zero
%     j = S(k);
%     [row,col]=co2xy(j,d);
%     % 0.5 (X-X(ij)Eij)' Q (X-X(ij)Eij) + P' (X-X(ij)Eij)
%     % 0.5 (X(ij)Eij)' Q (X(ij)Eij) - X' Q (X(ij)Eij) - P'(X(ij)Eij)
%     fobjs2(k) = 0.5*diagQ(row)*X(j)^2 -O(j)*X(j) -X(j)*P(j);
% end
rows = mod(S-1,d)+1;
fobjs2 = 0.5*diagQ(rows).*X(S).*X(S) - O(S).*X(S) - X(S).*P(S);
[~,ind1]=sort(fobjs1(:),'ascend');
[~,ind2]=sort(fobjs2(:),'ascend');
ind1_sort = Z(ind1);  
ind2_sort = S(ind2);  
n1 = min(length(ind2_sort),round(num/2));
index = [ind2_sort(1:n1);ind1_sort(1:num-n1)];



function [x_best] = quad_l0c_global(Q,p,k,list)
% This program solves the following l0 constrained problem
% min_x 0.5x'Qx + p'x, s.t. ||x||_0 = k
% Q: n x n
% p: n x 1
% k: 1 x 1

% Note:
% (1) Global optimum can be garanteed.
% (2) This program is practical only if n <= 30.

n = length(p);
% list = gen_list(n,k);
% list = getlist_random(n,k);
%  list = combs([1:n],k);


fs = [];
for iter = 1:size(list,1)
    S = list(iter,:);
    z = -Q(S,S)\p(S);
    fs(iter) = 0.5*z'*Q(S,S)*z + z'*p(S);
end

[~,index]=min(fs);
x_best = zeros(n,1);
S = list(index,:);
z = -Q(S,S)\p(S);
z(z==0) = 1e-14;
x_best(S) = z;


