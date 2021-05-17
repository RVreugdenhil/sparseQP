function [x,y]=SelectData(mynumber)
switch mynumber
    case 1
        load('a1a');
        [x,y] = get_random_subset(x,y) ;
    case 2
        load('w1a');
        [x,y] = get_random_subset(x,y) ;
    case 3
        load('w2a');
        [x,y] = get_random_subset(x,y) ;
    case 4
        load('madelon');
        [x,y] = get_random_subset(x,y) ;
    case 5
        load('mushroom');
        [x,y] = get_random_subset(x,y) ;
    case 6
        [x,y] = rnddata(300,100);
    case 7
        [x,y] = rnddata(300,500);
    case 8
        [x,y] = rnddata(300,1500);
    case 9
        [x,y] = rnddata(300,1800);
    case 10
        [x,y] = rnddata(300,2000);


end



function [x,y] = rnddata(n,d)
randn('seed',0);
x = randn(n,d);
y = ones(n,1);
y(1:(n/2)) = -1;

   

function [x,y] = get_random_subset(x,y) 
randn('seed',0);
rand('seed',0);
[n,d] = size(x);
max_iter = 2000;
if(n>max_iter),
    [seq] = randperm(n,max_iter);
    x = x(seq,:); y = y(seq);
end
x = full(x); y = full(y);
