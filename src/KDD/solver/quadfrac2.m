function [s,fs] = quadfrac2(u0,u1,u2,d0,d1,d2)
% min_{s}  (u0 + u1 s + u2 s^2)  /  (d0 + d1s + d2 s^2)
Handle = @(s) (u2.*s.*s + u1.*s + u0)  ./  (d2.*s.*s + d1.*s + d0);
% fprintf('%f %f %f %f %f %f\n',u0,u1,u2,d0,d1,d2);

% global Example
% for out = 1:length(u0),
%     
%     u22 = u2(out);
%     u11 = u1(out);
%     u00 = u0(out);
%     d22 = d2(out);
%     d11 = d1(out);
%     d00 = d0(out);
%     
%     xs = [-200:1:200];
%     fs1 = [];
%     for iter =1:length(xs),
%         Handle1 = @(s) (u22.*s.*s + u11.*s + u00)  ./  (d22.*s.*s + d11.*s + d00);
%         fs1(iter) = Handle1(xs(iter));
%     end
% 
% fs1
%     Example = [Example;fs1];
% end
% save ('Example','Example')
% 
% ddd

% grad_s =  (2 u2 s+u1)*(d2 s^2 + d1 s + d0) - (u2 s^2 + u1 s + u0)*(2 d2 s + d1) =0
% (2 u2 s + u1)*(d2 s^2 + d1 s + d0) - (2 d2 s + d1)*(u2 s^2 + u1 s + u0) =0
c2 = u2.*d1 - d2.*u1;
c1 = 2.*u2.*d0 - 2.*d2.*u0;
c0 = u1.*d0 - d1.*u0;

% c2 x^2 + c1 x + c0 = 0
dd = sqrt(c1.*c1 - 4*c2.*c0);
dd = real(dd);
x1 = (-c1+dd)./(2*c2);
x2 = (-c1-dd)./(2*c2);
ind = find(c2==0);
x1(ind) = -c0(ind)./c1(ind);
x2(ind) = x1(ind);
nanindex = find(isnan(x1) | isinf(x1));
x1(nanindex)=0; x2(nanindex)=0;
xx = [x1 x2]';
fs = [Handle(x1) Handle(x2)];
[fs,ind] = min(fs');
s = xx([0:2:(2*length(u0)-1)]+ind)';
fs = fs(:);

