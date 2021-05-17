function [x]=NonlinearFractionalProgrammingGlobal(x,Q,p,w,R,c,v)
% The program solves the following optimziation globally:
% min_{x} (0.5x'Qx+p'x+w) / (0.5x'Rx+c'x+v)
% x = sqrtm(inv(R)) * u - inv(R)c
% min_{u} (0.5u'Ou+g'u+0.5*delta) / (0.5u'u+0.5gamma)
[lb,ub,delta,gamma,g,O] = computeBound(Q,p,w,R,c,v);
iR = inv(R);
iRs = sqrtm(iR);
% save hhh Q p w R c v
u = solve_frac_problem(O,g,delta,gamma,lb,ub);
x = iRs * u - iR*c;

function [L,U,delta,gamma,g,O] = computeBound(Q,p,w,R,c,v)
iR = inv(R);
iRsqrt = sqrtm(iR);
O = iRsqrt*Q*iRsqrt';
gamma = 2*v - norm(iRsqrt*c)^2;
if(gamma==0),gamma=eps;end
g = iRsqrt*p - iRsqrt*Q*iR*c;
delta = c'*iR*Q*iR*c - 2*c'*iR*p + 2*w;
Z = [O g/sqrt(gamma);g'/sqrt(gamma) delta/gamma];
Z = (Z+Z')/2; O = (O+O')/2;
L = min(eig(Z));
U = min(eig(O));

function [u] = solve_frac_problem(O,g,delta,gamma,lb,ub)
% min_{u} (0.5u'Ou+g'u+0.5*delta) / (0.5u'u+0.5gamma)
% 0 = F(alpha) = min_{u} (0.5u'Ou+g'u+0.5*delta) - alpha (0.5u'u+0.5gamma)
% Ou + g - alpha u = 0
% (O - alpha I)u = -g
% u = - inv(O-alpha I) g

% 0 = (0.5u'(O-alphaI)u+g'u+0.5*delta) - alpha 0.5gamma
% 0 = -0.5g'inv(O-alphaI)g + 0.5*delta - alpha 0.5gamma
n = length(g);
I = eye(n);
% F = @(alpha)-0.5*g'*inv(O-alpha*I)*g + 0.5*delta - alpha*0.5*gamma;
[U,D]=eig((O+O')/2); d = diag(D);
a = U'*g;
F = @ (alpha)-0.5*sum((a.*a)./(d-alpha))+ 0.5*delta- alpha*0.5*gamma;

% Find alpha in [lb,ub] such that F(alpha) = 0. Note that F is a decreasing function
% [u] = solve_frac_problem1(O,g,delta,gamma,lb,ub);

alpha = bisection(F,lb,ub);

% if(abs(F(alpha))>0.1),
%     F(alpha)
%     lb
%     ub
%     dddddd
% end

u = - inv(O-alpha*I)*g;


% '*****************************************'
function [ret] = bisection(f,a,b)
% Find alpha in [lb,ub] such that f(alpha) = 0. Note that f is a decreasing function

eps1 = 1e-5;
fa = f(a);
fb = f(b);
ret = a;
if( (fa<0&&fb<0) || (fa>0 && fb>0) || (fa<0 && fb>0) || (abs(a-b)<eps1) || fa<eps1),
    if(abs(fa)<abs(fb)),ret = a;else ret = b;end
    return;
elseif(fa>0 && fb<0)
    
    if(~isLegal(fb))

        % Never happen
        for iiii=1:100
        b = b - eps1;
        fb = f(b);
        if(isLegal(fb)),break;end
        dddddd
        end
    end
    
    for iter = 1:100
        if(abs(fa)<eps1 || abs(fb)<eps1 || abs(a-b)<eps1)
           if(abs(fa)<abs(fb)),
               ret = a;
           else
               ret = b;
           end
        end
        c = (a+b)/2;
        fc = f(c);
        if(fc>0)
            a = c;
            fa = fc;
        elseif(fc<0)
            b = c;
            fb = fc;
        end
        fprintf('--------iter:%d, (a,b,c): %.2f %.2f %.2f, (fa fb fc): %.2e %.2e %.2e\n',iter,a,b,c,fa,fb,fc);
    end
    
end








