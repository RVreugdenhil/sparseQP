function [c] = find_k(a,k)
c = [];
for i=1:length(a)
    c = [c;a(i)];
    c = unique(c);
    if(length(c)==k),break;end
end
c = c(:);
