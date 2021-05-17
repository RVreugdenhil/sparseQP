function [A,B] = getFisher_AC(data_x,data_y,mean_y)
% data_x: m x d
% data_y: m x 1

% inx1 = find( y == 1);
% inx2 = find( y == 2);
% n1 = length(inx1);
% n2 = length(inx2);
% X = X';
% m1 = mean(X(:,inx1),2);
% m2 = mean(X(:,inx2),2);
% 
% S1 = (X(:,inx1))*(X(:,inx1))';
% S2 = (X(:,inx2))*(X(:,inx2))';
% 
% A = (m1-m2)*(m1-m2)';
% B = S1 + S2;
% 
% 
% 
% 
% A = full(A);
% B = full(B);
% C = A + B;



[num_data,dim] = size(data_x);
data_x = data_x';

inx1 = find( data_y >= mean_y);
inx2 = find( data_y < mean_y);
n1 = length(inx1);
n2 = length(inx2);


m1 = mean(data_x(:,inx1),2);
m2 = mean(data_x(:,inx2),2);

S1 = (data_x(:,inx1)-m1*ones(1,n1))*(data_x(:,inx1)-m1*ones(1,n1))';
S2 = (data_x(:,inx2)-m2*ones(1,n2))*(data_x(:,inx2)-m2*ones(1,n2))';
B = S1 + S2;
B = B / num_data;

A = (m1 - m2)*(m1 - m2)';


% min_eig = min(abs(eig(A)));
% sigma = max(1,abs(min_eig));
sigma = 0.5;
B = B + sigma * eye(size(A,1));


A = full(A);
B = full(B);
