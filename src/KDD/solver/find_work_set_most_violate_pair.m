function [index] = find_work_set_most_violate_pair(x,A,B,num)
% min_x x'Ax / x'Bx

if(num==0),index=[];return;end
n = size(A,1);
[Z] = find(x==0);
[S] = find(x~=0);

Ax = A*x;
Bx = B*x;
xAx = x'*Ax;
xBx = x'*Bx;
tt = xAx / xBx;
diagA = diag(A);
diagB = diag(B);
% decs=[];aa =[];bb =[];
% for i=1:length(Z)
%     for j=1:length(S)
%         ci = Z(i);
%         cj = S(j);
%         aa = [aa;ci];
%         bb = [bb;cj];
% 
%         % we set xj to zero exactly and set xi to nonzero
%         v = x; 
%         v(cj) = 0;
%         % min_t (v+t ei)' A (v+t ei) / (v+t ei)' B(v+t ei)
%         % min_t (v'Av + 2v'Aei t + ei'Aei tt) / (v'Bv + 2v'Bei t + ei'Bei tt)
%         Av = A*v;
%         Bv = B*v;
%         [~,ft] = quadfrac2(v'*Av,2*Av(ci),diagA(ci),v'*Bv,2*Bv(ci),diagB(ci));
%         dec  = ft - tt;
%         decs = [decs;dec];
% 
%     end
% end

aa = repmat(Z,1,length(S))';aa = aa(:);
bb = repmat(S,1,length(Z)); bb = bb(:);
u0 = xAx + 2*Ax(bb).*(-x(bb)) + diagA(bb).*(-x(bb)).*(-x(bb));
u1 = 2*(Ax(aa) + (-x(bb)).*A((bb-1)*n+aa));
u2 = diagA(aa);
d0 = xBx + 2*Bx(bb).*(-x(bb)) + diagB(bb).*(-x(bb)).*(-x(bb));
d1 = 2*(Bx(aa) + (-x(bb)).*B((bb-1)*n+aa));
d2 = diagB(aa);
[~,fs] = quadfrac2(u0,u1,u2,d0,d1,d2);
decs = fs - tt;

if(~isLegal(decs)),error('Wrong!');end
[~,ind]=sort(decs,'ascend');
aa = aa(ind); bb = bb(ind);
aa = unique(aa,'stable');
bb = unique(bb,'stable');
n2 = min(length(bb),floor(num/2));
n1 = min(length(aa),num-n2);
index = [aa(1:n1);bb(1:n2)];
