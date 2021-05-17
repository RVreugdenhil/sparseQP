function [A] = generate_A(m,n,type)
switch type
    case 1
        A = randn(m,n);
    case 2
        A = randn(m,n);
        seq = randperm(m*n,round(0.02*m*n));
        A(seq) = A(seq)*100;
end
