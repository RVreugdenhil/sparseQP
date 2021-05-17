function [list] = getlist_random(n,k)
times = 1000;
list = zeros(times,k);
for i=1:times,
    seq = randperm(n);
    seq = seq(1:k)';
    list(i,:) = seq;
end

