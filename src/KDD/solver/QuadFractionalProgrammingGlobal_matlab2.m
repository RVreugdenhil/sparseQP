function [x] = NonlinearFractionalProgrammingGB(x,A,b,c,D,e,f)
Handle = @(alpha)computeObj1(alpha,A,b,c,D,e,f);
% xx = [-100:1:100];
% for ii=1:length(xx),
% [ff,x] = Handle(xx(ii));
% end
[alpha]=fminsearch(Handle,0);
[fobj,x]=Handle(alpha);
alpha


function [fobj,x] = computeObj1(alpha,A,b,c,D,e,f)
% F(alpha) = 0.5x'Ax+b'x+c + alpha (0.5x'Dx+e'x+f)
% grad = A*x+b+alpha (D*x+e)=0
% (A+alpha D) x = -alpha e -b
x = pinv(A+alpha*D)*(-alpha*e-b);
fobj = 0.5*x'*A*x+b'*x+c + alpha*(0.5*x'*D*x+e'*x+f);
fobj = -fobj;
% y = alpha*e-b;
% H = inv(A-alpha*D);
% fobj = (c - alpha * f - 0.5 * y'*H*y)^2;

