function [index] = find_work_set_most_violate_coordinate(x,A,B,num)
% min_x x'Ax / x'Bx
[Z] = find(x==0);
[S] = find(x~=0);
Ax = A*x;
Bx = B*x;
xAx = x'*Ax;
xBx = x'*Bx;
diagA = diag(A);
diagB = diag(B);

% fobjs1 = zeros(length(Z),1);
% for k=1:length(Z)
%     % we change zero to nonzero
%     i = Z(k);
%     % min_t (x+t ei)' A (x+t ei) /  (x+t ei)' B (x+t ei)
%     % min_t (x'Ax+ 2 x'Aei t + tt A(i,i)) / (x'Bx+ 2 x'Bei t +  tt B(i,i))
%     [~,f] = quadfrac2(xAx,2*Ax(i),A(i,i),xBx,2*Bx(i),B(i,i));
%     fobjs1(k) = f;
% end

n1 = length(Z);
u0 = ones(n1,1)*xAx; u1 = 2*Ax(Z); u2 = diagA(Z);
d0 = ones(n1,1)*xBx; d1 = 2*Bx(Z); d2 = diagB(Z);
 
[~,fobjs1] = quadfrac2(u0,u1,u2,d0,d1,d2);

% fobjs2 = zeros(length(S),1);
% for k=1:length(S)
%     % we change nonzero to zero
%     j = S(k);
%     % min_t 0.5(x+t ei)' A (x+t ei) / 0.5 (x+t ei)' B (x+t ei)
%     t = -x(j);
%     fobjs2(k) = (xAx+2*Ax(j)*t+t*t*A(j,j)) / (xBx+2*Bx(j)*t+t*t*B(j,j)) ;
% end

t = -x(S);
fobjs2 = (xAx+2*Ax(S).*t + t.*t.*diagA(S)) ./ (xBx+2*Bx(S).*t+t.*t.*diagB(S)) ;
[~,ind1]=sort(fobjs1(:),'ascend');
[~,ind2]=sort(fobjs2(:),'ascend');
ind1_sort = Z(ind1);  ind1_sort = ind1_sort(:);
ind2_sort = S(ind2);  ind2_sort = ind2_sort(:);
n1 = min(length(ind2_sort),round(num/2));
index = [ind2_sort(1:n1);ind1_sort(1:num-n1)];
