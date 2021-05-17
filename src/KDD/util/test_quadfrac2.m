function test
for iter= 4:4,
    iter
    randn('seed',iter);
    u0 = randn(1);
    u1 = randn(1);
    u2 = randn(1);
    
    d0 = randn(1);
    d1 = randn(1);
    d2 = randn(1);
    
    x = quadfrac2(u0,u1,u2,d0,d1,d2)
end


function [s,fs] = quadfrac2(u0,u1,u2,d0,d1,d2)
% min_{s}  (u0 + u1 s + u2 s^2)  /  (d0 + d1s + d2 s^2)
Handle = @(s) (u2.*s.*s + u1.*s + u0)  ./  (d2.*s.*s + d1.*s + d0);
fprintf('%f %f %f %f %f %f\n',u0,u1,u2,d0,d1,d2);

% grad_s =  (2 u2 s+u1)*(d2 s^2 + d1 s + d0) - (u2 s^2 + u1 s + u0)*(2 d2 s + d1) =0
% (2 u2 s + u1)*(d2 s^2 + d1 s + d0) - (2 d2 s + d1)*(u2 s^2 + u1 s + u0) =0
c2 = u2.*d1 - d2.*u1;
c1 = 2.*u2.*d0 - 2.*d2.*u0;
c0 = u1.*d0 - d1.*u0;

% c2 x^2 + c1 x + c0 = 0
dd = sqrt(c1.*c1 - 4*c2.*c0);
x1 = (-c1+dd)./(2*c2);
x2 = (-c1-dd)./(2*c2);

x1
x2
ind = find(c2==0);
x1(ind) = -c0(ind)./c1(ind);
x2(ind) = x1(ind);

nanindex = find(isnan(x1) | isinf(x1));
x1(nanindex)=0; x2(nanindex)=0;


xx = [x1 x2]';

fs = [Handle(x1) Handle(x2)];

[fs,ind] = min(fs');

if(~isempty(find(ind~=1)))
    s = xx(ind);
(u2.*s.*s + u1.*s + u0)
    ddd
end
s = xx([0:2:(2*length(u0)-1)]+ind)';



