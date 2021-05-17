function [index] = find_work_set_most_violate_pair2(x,A,B,num)
% min_x x'Ax / x'Bx

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
% %         s = -x(cj);
%
%         % min_t 0.5(x+t ei+s ej)' A (x+t ei+s ej) / 0.5(x+t ei+s ej)' B(x+t ei+s ej)
%         %         y = x + s*ej;
% %         Ay_ci = Ax(ci) + s*A(ci,cj);
% %         yAy = xAx + 2*Ax(cj)*s + diagA(cj)*s*s;
% %         By_ci = Bx(ci) + s*B(ci,cj);
% %         yBy = xBx + 2*Bx(cj)*s + diagB(cj)*s*s;
%
%         % min_t (y+t ei)' A (y+t ei) / (y+t ei)' B(y+t ei)
%         % min_t (yAy + 2 y'Aei t + tt) / (yBy + 2 y'Bei t + tt)
% %         t = quadfrac(yAy,2*Ay_ci,diagA(ci),yBy,2*By_ci,diagB(ci));
% %         dec  = (yAy + 2*Ay_ci*t + t*t*diagA(ci)) / (yBy + 2*By_ci*t + t*t*diagB(ci)) - tt;
% %         decs = [decs;dec];
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
% index = [aa(1:num/2);bb(1:num/2)];
[index] = find_unique_k(aa,bb,num);

