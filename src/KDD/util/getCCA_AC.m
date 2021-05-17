function [A,B] = getCCA_AC(X,y)
% X: m x d
% y: m x 1

inx1 = find( y == 1);
inx2 = find( y == -1);
[C] = mycov(X');

XY = C(inx1,inx2);
XX = C(inx1,inx1);
YY = C(inx2,inx2);


[n1,n2]=size(XY);
A = [zeros(n1) XY;
     XY'  zeros(n2)];
B = [XX zeros(n1,n2);
      zeros(n2,n1) YY];
A = full(A);
B = full(B);
B = B + eye(size(B,1))*0.1;

function [C] = mycov(x)
m = size(x,1);
% xc = bsxfun(@minus,x,sum(x,1)/m);   C = (xc'*xc) / (m-1); 
xc = x; C = (xc'*xc) / m; 