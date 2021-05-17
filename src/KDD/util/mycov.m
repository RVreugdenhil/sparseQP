function [C] = mycov(x)
m = size(x,1);
% xc = bsxfun(@minus,x,sum(x,1)/m);   C = (xc'*xc) / (m-1); 
xc = x; C = (xc'*xc) / m; 