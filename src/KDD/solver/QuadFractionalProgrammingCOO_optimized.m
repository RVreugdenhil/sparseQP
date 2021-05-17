function [x,his] = NonlinearFractionalProgrammingCOO_optimized(x,A,b,c,Q,r,s,n,max_iter)
% min_x (0.5 x'Ax + b'x + c) / (0.5 x'Qx + x'r + s)


% x = randn(n,1);
HandleObj = @(x)(0.5*x'*A*x + b'*x + c) / (0.5*x'*Q*x + x'*r + s);
Ax = A*x; xAx = x'*Ax; bx = b'*x;
Qx = Q*x; xQx = x'*Qx; rx = r'*x;
diagA = diag(A);
diagQ = diag(Q);
his = [];
for iter = 1:max_iter
    up = 0.5*xAx + bx + c;
    down = 0.5*xQx + rx + s;
    up_g = Ax+b;
    down_g = Qx+r;
    fobj =  up / down;
       
%     grad = (up_g * down - up * down_g) / (down*down);
    grad = (up_g - fobj * down_g) /down;

    his(iter) = fobj;
    [~,i] = max(abs(grad));
   % i = iter;

%     if(abs(grad(i))<1e-5),break;end
%      i = randperm(n,1);
%      fprintf('iter:%d, fobj:%f, grad:%f, index:%d, down:%f\n',iter,fobj,grad(i),i,down);
    
   
    %              (0.5 (x+alpha e)'A(x+alpha e) + b'(x+alpha e) + c)
    % min_{alpha} --------------------------------------------------
    %              (0.5 (x+alpha e)'Q(x+alpha e) + r'(x+alpha e) + s)
    
    %              (0.5x'Ax + x'Ae alpha + 0.5 e'Ae alpha^2 + b'x+alpha b'e + c)
    % min_{alpha} --------------------------------------------------------------
    %              (0.5x'Qx + x'Qe alpha + 0.5 e'Qe alpha^2 + r'x+alpha r'e + s)
    
	u0 = 0.5*xAx + bx + c;
    u1 = Ax(i) + b(i);
	u2 = 0.5*diagA(i);
	d0 = 0.5*xQx + rx + s;
	d1 = Qx(i) + r(i);
	d2 = 0.5*diagQ(i);

    alpha =  quadfrac2(u0,u1,u2,d0,d1,d2);
    x(i) = x(i) + alpha;

    % x => x + alpha e
    % reconstruct Ax, xAx, bx, Qx, xQx, rx
    xAx = xAx + 2*alpha*Ax(i)+alpha*alpha*diagA(i);
    Ax = Ax + alpha*A(:,i);
    bx = bx + alpha*b(i);
    xQx = xQx + 2*alpha*Qx(i)+alpha*alpha*diagQ(i);
    Qx = Qx + alpha*Q(:,i);
    rx = rx + alpha*r(i);

%   fprintf('%f %f %f %f %f %f %f| xAx:%f %d\n',u0,u1,u2,d0,d1,d2,alpha, xAx,i);

end
x(x==0)=eps;

function [x,his] = NonlinearFractionalProgrammingCOO(x,A,b,c,Q,r,s)
% min_x (0.5 x'Ax + b'x + c) / (0.5 x'Qx + x'r + s)

n = length(x);
HandleObj = @(x)(0.5*x'*A*x + b'*x + c) / (0.5*x'*Q*x + x'*r + s);

for iter = 1:20,
%     i = randperm(n,1);
    [fff,ggg]=computeObj(x,A,b,c,Q,r,s);
      his(iter) = fff;
    [~,i] = max(abs(ggg));
    e = zeros(n,1);
    e(i) = 1;
    %              (0.5 (x+alpha e)'A(x+alpha e) + b'(x+alpha e) + c)
    % min_{alpha} --------------------------------------------------
    %              (0.5 (x+alpha e)'Q(x+alpha e) + r'(x+alpha e) + s)
    
    %              (0.5x'Ax + x'Ae alpha + 0.5 e'Ae alpha^2 + b'x+alpha b'e + c)
    % min_{alpha} --------------------------------------------------------------
    %              (0.5x'Qx + x'Qe alpha + 0.5 e'Qe alpha^2 + r'x+alpha r'e + s)
   alpha = quadfrac2(0.5*x'*A*x + b'*x + c, x'*A*e + b'*e, 0.5*e'*A*e,0.5*x'*Q*x + r'*x + s, x'*Q*e + r'*e, 0.5*e'*Q*e);
    x(i) = x(i) + alpha;
%     fobj = HandleObj(x);
    %     if(~mod(iter,100))
%          fprintf('iter:%d, fobj:%f, ggg:%f\n',iter,fobj,norm(ggg));
    %     end
      
end
x(x==0)=eps;
