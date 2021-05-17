function [b] = find_unique_k(a,k)
b = [];
for i=1:length(a)
    b = [b;a(i)];
    b = unique(b);
    if(length(b)>=k),break;end
end

b = b(:);