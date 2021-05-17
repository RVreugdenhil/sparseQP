function [B] = find_work_set_random(x,num)
[Z]=find(x==0);
[S]=find(x~=0);
lenZ = length(Z);
lenS = length(S);
n1 = min(lenZ,floor(num/2));
n2 = min(lenS,num-n1);
iZ = randperm(lenZ,n1);
iS = randperm(lenS,n2);
B = [Z(iZ);S(iS)];
