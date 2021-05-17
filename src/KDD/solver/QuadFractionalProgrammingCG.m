function [x,his] = NonlinearFractionalProgrammingCG(x,A,b,c,D,e,f)
% min_x (0.5 x'Ax + b'x + c) / (0.5 x'Dx + x'e + f)
% x = x/norm(x);
n = length(x);
x = randn(n,1);
x0 = x;
HandleObj = @(x)computeObj(x,A,b,c,D,e,f);
% x = randn(size(x));
% maxIter=1000;
% [x]= smoothopt(x0,HandleObj,maxIter) ;
% return;

if(norm(x)<1e-5),x = randn(size(x));end
max_iter = 30;
his=[];
method = 1;
for iter = 1:max_iter
    [fobj_cur,grad_cur]=HandleObj(x);
    
    if(norm(grad_cur)<1e-8 || ~isLegal(fobj_cur) || ~isLegal(grad_cur)),break;end
    
    his = [his;fobj_cur];
%     fprintf('iter: %d, fobj: %f, grad: %e\n',iter,fobj_cur,norm(grad_cur));
    if(method==1)
        if(iter==1)
            direction= - grad_cur ;
        else
            %          beta= mdot(grad_cur,grad_cur)/ mdot(grad_old,grad_old);                     % FR
            beta=  mdot(grad_cur,grad_cur-grad_old) / mdot(grad_old,grad_old);            % PRP
            beta(beta<0)=0;
            direction= - grad_cur + beta * direction;
        end
        grad_old = grad_cur;
    elseif(method==2)
        if(iter==1)
            direction= - grad_cur ;
            old_dirs = zeros(n,1,0);
            old_stps = zeros(n,1,0);
            Hdiag = 1;
        else
            [old_dirs,old_stps,Hdiag] = lbfgsUpdate(grad_cur-grad_old,x-x_old,old_dirs,old_stps,Hdiag);
            direction = lbfgs(-grad_cur,old_dirs,old_stps,Hdiag);
        end;
        grad_old = grad_cur;
        x_old = x;
    end
    
    
    % min_s (0.5 (x+sd)'A(x+sd) + b'(x+sd) + c) / (0.5 (x+sd)'D(x+sd) + (x+sd)'e + f)
    % min_s (0.5 d'Ad ss + x'Ad s + 0.5 x'Ax + b'x + s b'd + c) / (0.5 d'Dd ss + x'Dd s + 0.5 x'Dx + e'x + s e'd + f)
    % min_s  (u0 + u1 s + u2 s^2)  /  (d0 + d1s + s2 s^2)
    u0 = 0.5*x'*A*x + b'*x + c;
    u1 = x'*A*direction + b'*direction;
    u2 = 0.5*direction'*A*direction;
    
    d0 = 0.5*x'*D*x + e'*x + f;
    d1 = x'*D*direction + e'*direction;
    d2 = 0.5*direction'*D*direction;
    
    [s] = quadfrac(u0,u1,u2,d0,d1,d2) ;
    
    if(~isLegal(s)),break;end
    
    x = x + s*direction;
    
    
end
