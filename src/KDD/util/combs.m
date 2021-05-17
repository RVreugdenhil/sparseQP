function P = combs(v,m)
% combs computes all possible combinations.
v = v(:)';
n = length(v);
if n == m
    P = v;
elseif m == 1
    P = v';
else
    P = [];
    if m < n && m > 1
        for k = 1:n-m+1
            Q = combs(v(k+1:n),m-1);
            P = [P; [v(ones(size(Q,1),1),k) Q]];
        end
    end
end