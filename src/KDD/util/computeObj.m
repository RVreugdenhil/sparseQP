function [fobj,grad] = computeObj(x,A,b,c,D,e,f)
% min_x (0.5x'Ax + b'x + c) / (0.5x'Dx + x'e + f)
up = (0.5*x'*A*x + b'*x + c);
down = (0.5*x'*D*x + x'*e + f);
up_g = A*x+b;
down_g = D*x+e;
fobj =  up / down;
grad = (up_g * down - up * down_g) / (down*down);
