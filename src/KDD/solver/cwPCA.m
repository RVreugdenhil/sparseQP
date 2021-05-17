function [x_out,fun_all,tt,time_spend] = cwPCA(A,s,varargin)
%This function implements a coordinate-wise algorithm for solving the 
% sparse principal component analysis (PCA) problem.
%
% Based on the paper
% Amir Beck and Yakov Vaisbourd, "The Sparse Principal Component Analysis
%                           Problem: Optimality Conditions and Algorithms"
% -----------------------------------------------------------------------
% Copyright (2016): Amir Beck and Yakov Vaisbourd
% 
% cwPCA is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
%% INPUT
%
% A ............................. The covariance or the data matrix.
%
% s ............................. The sparsity level.
%
%
% OPTIONAL INPUTS   
% This arguments should be provided by the user as a - name, value - pairs.
%
% 'initial_support', vector ... This vector contains the initial
%                               support. Default: (1,2,...,s)';
%
% 'mat_type', 'CMat' ....... Covariance\correlation matrix (default).
% 'mat_type', 'DMat' ....... Data matrix.
%
% 'type', 'GCW' ............ Greedy CW (GCW) (default).
% 'type', 'PCW' ............ Partial CW (PCW).
%
% 'max_iter', scalar ....... Maximum number of iterations. If applied
%                            then there is no guarantee that the solution 
%                            is CW maximal. (default: 0 - Suppressed)).
% 
% 'display', true .......... Display the report on the iterations. 
% 'display', false ......... Silence the report on the iterations (default).
%
%
%% OUTPUT
% 
% X_out .............................. Solution to the sparse PCA problem.
%
% fun_all ............................ Array containing the sequence of
%                                      function values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                   




if nargin < 1
    error 'Wrong Input Parameters: Covariance\Data matrix is not defined.'
end
tbegin = cputime;
% Assign default values
display = false;
[~,n] = size(A);
max_iter = 0;
type = 'GCW';
Support = 1:s;
mat_type = 'CMat';

% Update variables according to input parameters
if (rem(length(varargin),2)==1)
    error('Optional parameters should always provided in pairs');
else
    for i = 1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'initial_support';  Support       = varargin{i+1};             
            case 'display';          display        = varargin{i+1};    
            case 'max_iter';         max_iter      = varargin{i+1};  
            case 'type';             type          = varargin{i+1}; 
            case 'mat_type';         mat_type      = varargin{i+1};             
        otherwise
                error(['Unrecognized option: ''', varargin{i}, '''']);
        end
    end
end

% Input validation: Optional parameters
if ~strcmp(type,'GCW') && ~strcmp(type,'PCW')
    error 'Undefined method type'
end

if numel(Support)<1 || numel(Support)>s
    error 'The size of the initial support must be between 1 to s.'
end

if sum(mod(Support,1)>0) || max(Support)>n || min(Support)<0
    str = ['The support should contain integer indices with values',...
            'between 1 and n.'];
	error(str);
end

if (s >= n || s<1 || mod(s,1)>0)
    error 'Sparsity Varaible should be an integer between 1 and n.'
end

if max_iter < 0 || mod(max_iter,1)>0
    error 'max_iter must be nonnegative integer.'
end

if ~strcmp(mat_type,'CMat') && ~strcmp(mat_type,'DMat')
    error 'Unrecognized matrix type.'
end

if ~islogical(display)
    error 'The display parameter should be logical.'
end

if display
    fprintf('\n Starting cwPCA method with type: %s',type);
    fprintf(' \n');
    fprintf('----------------------------------------\n');
    fprintf('|   k  |  out  |  in  |       val      |\n');
    fprintf('----------------------------------------\n');
end

% Initialization
iter=0;
 
fun_all = [];
 
tt = [];
t1 = clock();

if strcmp(mat_type,'CMat')
    diagA = diag(A);
    oracle = @(S) eigs(A(S,S),1);
else
    A = A./sqrt(n-1);   
    diagA = sum(A.^2)';
    oracle = @(S) eigs(A(:,S)'*A(:,S),1);
end

% Compute SO point that correspond to the initial support
if numel(Support)>0    
    [x_s,f_val] = oracle(Support);
    temp = x_s==0;
    if sum(temp)>0        
        Support(temp) = [];
        x_s(temp) = [];
    end
    OSupport = setdiff(1:n,Support);
else
    f_val = 0;
    OSupport = 1:n;
end

% Main loop
flag = true;
while flag && (max_iter==0 || iter < max_iter)
    iter = iter+1;
    
    % Phase one: The support is smaller then the sparsity level
    
    if numel(Support)<s
        f_max = f_val;
        j_max = 0;
        for j=1:numel(OSupport)            
            [v,val] = oracle([Support,OSupport(j)]);
            if val>f_max
                j_max = j;
                v_max = v;
                f_max = val;
            end
        end
        if j_max>0    
            index_in = OSupport(j_max);
            x_s = v_max;
            Support = [Support,index_in];               %#ok
            f_val = f_max;
            OSupport(j_max) = [];
            temp = x_s==0;
            if sum(temp)>0
                OSupport = [OSupport,Support(temp)];    %#ok
                Support(temp) = [];
                x_s(temp) = [];            
            end
            if display
                fprintf('|%5.0f |  NaN  | %5.0f | %5.10f  |\n',iter,...
                                            index_in,f_val); 
            end
            continue;
        end
    end      
    
                 
    % Phase two:
    
    % Compute the gradient
    if strcmp(mat_type,'CMat')            
        g = 2*A(:,Support)*x_s;
    else
        temp = A(:,Support)*x_s;
        g = 2*A'*temp;            
    end
    
    % Update parameters
    f_dif_max = 0;               
    i_max = 0;
    j_max = 0;
    if strcmp(type,'GCW')        
        checkSet = Support;        
    else
        temp = sortrows([abs(x_s),Support',(1:numel(Support))'],1);        
        checkSet = temp(:,2)';
        ordered_perm = temp(:,3)';
    end
    
    % Find swap
    for k=1:numel(checkSet)
        i = checkSet(k);    
        if strcmp(type,'GCW') 
            i_s = k;
        else
            i_s = ordered_perm(k);
        end
        if strcmp(mat_type,'CMat')
            DDi = A(:,i);
            h = g-2*x_s(i_s)*A(:,i);
        else
            DDi = A'*A(:,i);
            h = g-2*x_s(i_s)*DDi;            
        end
        val = (diagA(OSupport)+diagA(i)-2*sign(g(i))*...
              sign(h(OSupport)).*DDi(OSupport))*x_s(i_s)^2 + ...
              sign(h(OSupport)).*g(OSupport)*abs(x_s(i_s)) - g(i)*x_s(i_s);
               
        [val,j] = max(val);
        if val > f_dif_max
            if strcmp(type,'GCW')
                i_max = k;
                j_max = j;
                f_dif_max = val;
            else
                i_max = ordered_perm(k);
                j_max = j;
                break;
            end
        end
    end
    
    % Perform a swap and find an SO point
    
    if i_max==0
        flag = false;
    else
        index_out = Support(i_max);
        index_in = OSupport(j_max);
        Support(i_max) = index_in;
        OSupport(j_max) = index_out;
        [x_s,f_val] = oracle(Support);
        temp = x_s==0;
        if sum(temp)>0
            OSupport = [OSupport,Support(temp)];    %#ok
            Support(temp) = [];
            x_s(temp) = [];            
        end
    end
    
    
    
     
 
            fun_all = [fun_all;f_val];
               tt = [tt;etime(clock(),t1)];
   
  
                
    if display
        fprintf('|%5.0f | %5.0f | %5.0f | %5.10f  |\n',iter,index_out, ...
                                                           index_in,f_val); 
    end
 
    
end
x_out = zeros(n,1);
x_out(Support) = x_s;
 
 
tend = cputime;
  time_spend = tend - tbegin;
  

fun_all = -fun_all;
end