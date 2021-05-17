
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Souce Code for smoothopt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, histroy,times]=smoothopt(x,HandleObj,maxIter,varargin)
randn('state',0);
[n,d] = size(x);
histroy=[];
iter=0;
times = [];
t1 = clock;

% gradmethod = 1; % Conjugate Gradient
gradmethod = 2; % L-bfgs

 linesearch = 1; % WolfeLineSearch
%  linesearch = 2; % ArmijoBacktrack

while 1
    iter=iter+1;
    if(iter>maxIter),
        fprintf('exit!(max iteration:%d)\n',iter);
        break;end
    [fobj_cur,grad_cur]=HandleObj(x);t2 = clock;
    second_spent = etime(t2,t1);
    times = [times;second_spent];
    
    
    histroy=[histroy;fobj_cur];
    fprintf('it:%d obj:%.2e norm:%.2e \n',iter,fobj_cur,norm(grad_cur,'fro'))
    
    if (gradmethod==1)
        if(iter==1)
            direction= - grad_cur ;
        else
            %      beta= mdot(grad_cur,grad_cur)/ mdot(grad_old,grad_old);                     % FR
            beta=  mdot(grad_cur,grad_cur-grad_old) / mdot(grad_old,grad_old);            % PRP
            beta(beta<0)=0;
            direction= - grad_cur + beta *   direction;
        end
    elseif(gradmethod==2)
        
        if(iter==1)
            direction= - grad_cur ;
            old_dirs = zeros(n,d,0);
            old_stps = zeros(n,d,0);
            Hdiag = 1;
        else
            
            [old_dirs,old_stps,Hdiag] = lbfgsUpdate(grad_cur-grad_old,x-x_old,old_dirs,old_stps,Hdiag);
            direction= lbfgs(-grad_cur,old_dirs,old_stps,Hdiag);

            
            
        end
        
        
    end
    

    gtd = mdot(direction,grad_cur);
    if(gtd>0),
        direction= - grad_cur;
    end
    if(abs(gtd)<1e-14),
        fprintf('exit!(gradient test too small:%f)\n',gtd);
        break;
    end
    
    
    
    grad_old = grad_cur;
    x_old=x;
    
    if(linesearch==1)
        [stepsize,f1,g1] = WolfeLineSearch(x,1,direction,fobj_cur,grad_cur,gtd,1e-4,0.9,2,0,25,1e-9,0,0,1,HandleObj,varargin{:});
    else
        [stepsize,x_new,f1,g1] = ArmijoBacktrack(x,1,direction,fobj_cur,fobj_cur,grad_cur,gtd,1e-4,2,1e-9,0,0,HandleObj,varargin{:});
    end
    
    
    if(stepsize<0 || stepsize<1e-14),
        fprintf('exit!(step size too small:%e)\n',stepsize);
        break;
    end
    
    
    optCond = max(max(abs(g1)));
    if optCond <= 1e-12
        fprintf('optCond too small:%e\n',optCond);
        break;
    end
    
    step = max(max(abs(stepsize*direction)));
    if step <= 1e-12
        fprintf('Step Size too small:%f\n',step);
        break;
    end
    
    
    x=x+stepsize*direction;
end;





function [d] = lbfgs(g,s,y,Hdiag)
% BFGS Search Direction
%
% This function returns the (L-BFGS) approximate inverse Hessian,
% multiplied by the gradient
% inv(H)*grad
% If you pass in all previous directions/sizes, it will be the same as full BFGS
% If you truncate to the k most recent directions/sizes, it will be L-BFGS
%
% s - previous search directions (n by k)
% y - previous step sizes (n by k)
% g - gradient (n by 1)
% Hdiag - value of initial Hessian diagonal elements (scalar)

[n,d,k] = size(s);

for i = 1:k
    ro(i,1) = 1 / mdot(y(:,:,i),s(:,:,i));
end

q = zeros(n,d,k+1);
r = zeros(n,d,k+1);
al =zeros(k,1);
be =zeros(k,1);

q(:,:,k+1) = g;

for i = k:-1:1
    al(i) = ro(i) * mdot(s(:,:,i),q(:,:,i+1));
    q(:,:,i) = q(:,:,i+1) - al(i)*y(:,:,i);
end

% Multiply by Initial Hessian
r(:,:,1) = Hdiag*q(:,:,1);

for i = 1:k
    be(i) = ro(i)*mdot(y(:,:,i),r(:,:,i));
    r(:,:,i+1) = r(:,:,i) + s(:,:,i)*(al(i)-be(i));
end
d=r(:,:,k+1);


function [old_dirs,old_stps,Hdiag] = lbfgsUpdate(y,s,old_dirs,old_stps,Hdiag)
corrections=10;
debug=0;

ys = mdot(y,s);
if ys > 1e-10,
    numCorrections = size(old_dirs,3);%看有多少列
    if numCorrections < corrections
        % Full Update
        old_dirs(:,:,numCorrections+1) = s;
        old_stps(:,:,numCorrections+1) = y;
    else
        % Limited-Memory Update
        
        for i=1:corrections-1,
            old_dirs(:,:,i)=old_dirs(:,:,i+1);
            old_stps(:,:,i)=old_stps(:,:,i+1);
            
        end;
        old_dirs(:,:,corrections)=s;
        old_stps(:,:,corrections)=y;
    end
    
    % Update scale of initial Hessian approximation
    Hdiag = ys/mdot(y,y);
else
    if(debug)
        fprintf('Skipping Update\n');
    end
end





function [t,f_new,g_new,funEvals] = WolfeLineSearch(x,t,d,f,g,gtd,c1,c2,LS_interp,LS_multi,maxLS,progTol,debug,doPlot,saveHessianComp,funObj,varargin)
%
% Bracketing Line Search to Satisfy Wolfe Conditions
%
% Inputs:
%   x: starting location
%   t: initial step size
%   d: descent direction
%   f: function value at starting location
%   g: gradient at starting location
%   gtd: directional derivative at starting location
%   c1: sufficient decrease parameter
%   c2: curvature parameter
%   debug: display debugging information
%   LS_interp: type of interpolation
%   maxLS: maximum number of iterations
%   progTol: minimum allowable step length
%   doPlot: do a graphical display of interpolation
%   funObj: objective function
%   varargin: parameters of objective function
%
% Outputs:
%   t: step length
%   f_new: function value at x+t*d
%   g_new: gradient value at x+t*d
%   funEvals: number function evaluations performed by line search


% Evaluate the Objective and Gradient at the Initial Step

[f_new,g_new] = funObj(x+t*d,varargin{:});

funEvals = 1;
gtd_new = mdot(g_new,d);

% Bracket an Interval containing a point satisfying the
% Wolfe criteria

LSiter = 0;
t_prev = 0;
f_prev = f;
g_prev = g;
gtd_prev =  (gtd);
nrmD = max(max(abs(d)));
done = 0;

while LSiter < maxLS
    
    %% Bracketing Phase
    if ~isLegal(f_new) || ~isLegal(g_new)
        if debug
            fprintf('Extrapolated into illegal region, switching to Armijo line-search\n');
        end
        t = (t + t_prev)/2;
        % Do Armijo
        
        [t,x_new,f_new,g_new,armijoFunEvals] = ArmijoBacktrack(x,t,d,f,f,g,gtd,c1,LS_interp,LS_multi,progTol,debug,doPlot,...
            funObj,varargin{:});
        
        funEvals = funEvals + armijoFunEvals;
        return;
    end
    
    
    if f_new > f + c1*t*gtd || (LSiter > 1 && f_new >= f_prev)
        bracket = [t_prev t];
        bracketFval = [f_prev f_new];
        %         bracketGval = [g_prev g_new];
        bracketGval(:,:,1) = g_prev;bracketGval(:,:,2) = g_new;
        break;
    elseif abs(gtd_new) <= -c2*gtd
        bracket = t;
        bracketFval = f_new;
        bracketGval(:,:,1) =  g_new;
        done = 1;
        break;
    elseif gtd_new >= 0
        bracket = [t_prev t];
        bracketFval = [f_prev f_new];
        bracketGval(:,:,1) =  g_prev  ;
        bracketGval(:,:,2) =g_new;
        break;
    end
    temp = t_prev;
    t_prev = t;
    minStep = t + 0.01*(t-temp);
    maxStep = t*10;
    if LS_interp <= 1
        if debug
            fprintf('Extending Braket\n');
        end
        t = maxStep;
    elseif LS_interp == 2
        if debug
            fprintf('Cubic Extrapolation\n');
        end
        t = polyinterp([temp f_prev gtd_prev; t f_new gtd_new],doPlot,minStep,maxStep);
    elseif LS_interp == 3
        t = mixedExtrap(temp,f_prev,gtd_prev,t,f_new,gtd_new,minStep,maxStep,debug,doPlot);
    end
    
    f_prev = f_new;
    g_prev = g_new;
    gtd_prev = gtd_new;
    
    
    [f_new,g_new] = funObj(x + t*d,varargin{:});
    
    funEvals = funEvals + 1;
    gtd_new = mdot(g_new,d);
    LSiter = LSiter+1;
end

if LSiter == maxLS
    bracket = [0 t];
    bracketFval = [f f_new];
    %     bracketGval = [g g_new];
    bracketGval(:,:,1) = g;
    bracketGval(:,:,2) = g_new;
    
end

%% Zoom Phase

% We now either have a point satisfying the criteria, or a bracket
% surrounding a point satisfying the criteria
% Refine the bracket until we find a point satisfying the criteria
insufProgress = 0;
Tpos = 2;
LOposRemoved = 0;
while ~done && LSiter < maxLS
    
    % Find High and Low Points in bracket
    [f_LO LOpos] = min(bracketFval);
    HIpos = -LOpos + 3;
    
    % Compute new trial value
    if LS_interp <= 1 || ~isLegal(bracketFval) || ~isLegal(bracketGval)
        if debug
            fprintf('Bisecting\n');
        end
        t = mean(bracket);
    elseif LS_interp == 2
        if debug
            fprintf('Grad-Cubic Interpolation\n');
        end
        t = polyinterp([bracket(1) bracketFval(1) mdot(bracketGval(:,:,1),d)
            bracket(2) bracketFval(2) mdot(bracketGval(:,:,2),d)],doPlot);
    else
        % Mixed Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nonTpos = -Tpos+3;
        if LOposRemoved == 0
            oldLOval = bracket(nonTpos);
            oldLOFval = bracketFval(nonTpos);
            oldLOGval = bracketGval(:,:,nonTpos);
        end
        t = mixedInterp(bracket,bracketFval,bracketGval,d,Tpos,oldLOval,oldLOFval,oldLOGval,debug,doPlot);
    end
    
    
    % Test that we are making sufficient progress
    if min(max(bracket)-t,t-min(bracket))/(max(bracket)-min(bracket)) < 0.1
        if debug
            fprintf('Interpolation close to boundary');
        end
        if insufProgress || t>=max(bracket) || t <= min(bracket)
            if debug
                fprintf(', Evaluating at 0.1 away from boundary\n');
            end
            if abs(t-max(bracket)) < abs(t-min(bracket))
                t = max(bracket)-0.1*(max(bracket)-min(bracket));
            else
                t = min(bracket)+0.1*(max(bracket)-min(bracket));
            end
            insufProgress = 0;
        else
            if debug
                fprintf('\n');
            end
            insufProgress = 1;
        end
    else
        insufProgress = 0;
    end
    
    % Evaluate new point

        [f_new,g_new] = funObj(x + t*d,varargin{:});

    funEvals = funEvals + 1;
    gtd_new = mdot(g_new,d);
    LSiter = LSiter+1;
    
    armijo = f_new < f + c1*t*gtd;
    if ~armijo || f_new >= f_LO
        % Armijo condition not satisfied or not lower than lowest
        % point
        bracket(HIpos) = t;
        bracketFval(HIpos) = f_new;
        bracketGval(:,:,HIpos) = g_new;
        Tpos = HIpos;
    else
        if abs(gtd_new) <= - c2*gtd
            % Wolfe conditions satisfied
            done = 1;
        elseif gtd_new*(bracket(HIpos)-bracket(LOpos)) >= 0
            % Old HI becomes new LO
            bracket(HIpos) = bracket(LOpos);
            bracketFval(HIpos) = bracketFval(LOpos);
            bracketGval(:,:,HIpos) = bracketGval(:,:,LOpos);
            if LS_interp == 3
                if debug
                    fprintf('LO Pos is being removed!\n');
                end
                LOposRemoved = 1;
                oldLOval = bracket(LOpos);
                oldLOFval = bracketFval(LOpos);
                oldLOGval = bracketGval(:,:,LOpos);
            end
        end
        % New point becomes new LO
        bracket(LOpos) = t;
        bracketFval(LOpos) = f_new;
        bracketGval(:,:,LOpos) = g_new;
        Tpos = LOpos;
    end
    
    if ~done && abs(bracket(1)-bracket(2))*nrmD < progTol
        if debug
            fprintf('Line-search bracket has been reduced below progTol\n');
        end
        break;
    end
    
end

%%
if LSiter == maxLS
    if debug
        fprintf('Line Search Exceeded Maximum Line Search Iterations\n');
    end
end

[f_LO LOpos] = min(bracketFval);
t = bracket(LOpos);
f_new = bracketFval(LOpos);
g_new = bracketGval(:,:,LOpos);




function [t] = mixedExtrap(x0,f0,g0,x1,f1,g1,minStep,maxStep,debug,doPlot);
alpha_c = polyinterp([x0 f0 g0; x1 f1 g1],doPlot,minStep,maxStep);
alpha_s = polyinterp([x0 f0 g0; x1 sqrt(-1) g1],doPlot,minStep,maxStep);
if alpha_c > minStep && abs(alpha_c - x1) < abs(alpha_s - x1)
    if debug
        fprintf('Cubic Extrapolation\n');
    end
    t = alpha_c;
else
    if debug
        fprintf('Secant Extrapolation\n');
    end
    t = alpha_s;
end


%%
function [t] = mixedInterp(bracket,bracketFval,bracketGval,d,Tpos,oldLOval,oldLOFval,oldLOGval,debug,doPlot);

% Mixed Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nonTpos = -Tpos+3;

gtdT = mdot(bracketGval(:,:,Tpos),d);
gtdNonT = mdot(bracketGval(:,:,nonTpos),d);
oldLOgtd = mdot(oldLOGval,d);
if bracketFval(Tpos) > oldLOFval
    alpha_c = polyinterp([oldLOval oldLOFval oldLOgtd
        bracket(Tpos) bracketFval(Tpos) gtdT],doPlot);
    alpha_q = polyinterp([oldLOval oldLOFval oldLOgtd
        bracket(Tpos) bracketFval(Tpos) sqrt(-1)],doPlot);
    if abs(alpha_c - oldLOval) < abs(alpha_q - oldLOval)
        if debug
            fprintf('Cubic Interpolation\n');
        end
        t = alpha_c;
    else
        if debug
            fprintf('Mixed Quad/Cubic Interpolation\n');
        end
        t = (alpha_q + alpha_c)/2;
    end
elseif gtdT'*oldLOgtd < 0
    alpha_c = polyinterp([oldLOval oldLOFval oldLOgtd
        bracket(Tpos) bracketFval(Tpos) gtdT],doPlot);
    alpha_s = polyinterp([oldLOval oldLOFval oldLOgtd
        bracket(Tpos) sqrt(-1) gtdT],doPlot);
    if abs(alpha_c - bracket(Tpos)) >= abs(alpha_s - bracket(Tpos))
        if debug
            fprintf('Cubic Interpolation\n');
        end
        t = alpha_c;
    else
        if debug
            fprintf('Quad Interpolation\n');
        end
        t = alpha_s;
    end
elseif abs(gtdT) <= abs(oldLOgtd)
    alpha_c = polyinterp([oldLOval oldLOFval oldLOgtd
        bracket(Tpos) bracketFval(Tpos) gtdT],...
        doPlot,min(bracket),max(bracket));
    alpha_s = polyinterp([oldLOval sqrt(-1) oldLOgtd
        bracket(Tpos) bracketFval(Tpos) gtdT],...
        doPlot,min(bracket),max(bracket));
    if alpha_c > min(bracket) && alpha_c < max(bracket)
        if abs(alpha_c - bracket(Tpos)) < abs(alpha_s - bracket(Tpos))
            if debug
                fprintf('Bounded Cubic Extrapolation\n');
            end
            t = alpha_c;
        else
            if debug
                fprintf('Bounded Secant Extrapolation\n');
            end
            t = alpha_s;
        end
    else
        if debug
            fprintf('Bounded Secant Extrapolation\n');
        end
        t = alpha_s;
    end
    
    if bracket(Tpos) > oldLOval
        t = min(bracket(Tpos) + 0.66*(bracket(nonTpos) - bracket(Tpos)),t);
    else
        t = max(bracket(Tpos) + 0.66*(bracket(nonTpos) - bracket(Tpos)),t);
    end
else
    t = polyinterp([bracket(nonTpos) bracketFval(nonTpos) gtdNonT
        bracket(Tpos) bracketFval(Tpos) gtdT],doPlot);
end



function [legal] = isLegal(v)
legal = sum(any(imag(v(:))))==0 & sum(isnan(v(:)))==0 & sum(isinf(v(:)))==0;



function [minPos,fmin] = polyinterp(points,doPlot,xminBound,xmaxBound)
% function [minPos] = polyinterp(points,doPlot,xminBound,xmaxBound)
%
%   Minimum of interpolating polynomial based on function and derivative
%   values
%
%   It can also be used for extrapolation if {xmin,xmax} are outside
%   the domain of the points.
%
%   Input:
%       points(pointNum,[x f g])
%       doPlot: set to 1 to plot, default: 0
%       xmin: min value that brackets minimum (default: min of points)
%       xmax: max value that brackets maximum (default: max of points)
%
%   set f or g to sqrt(-1) if they are not known
%   the order of the polynomial is the number of known f and g values minus 1

if nargin < 2
    doPlot = 0;
end

nPoints = size(points,1);
order = sum(sum((imag(points(:,2:3))==0)))-1;

xmin = min(points(:,1));
xmax = max(points(:,1));

% Compute Bounds of Interpolation Area
if nargin < 3
    xminBound = xmin;
end
if nargin < 4
    xmaxBound = xmax;
end

% Code for most common case:
%   - cubic interpolation of 2 points
%       w/ function and derivative values for both

if nPoints == 2 && order ==3 && doPlot == 0
    % Solution in this case (where x2 is the farthest point):
    %    d1 = g1 + g2 - 3*(f1-f2)/(x1-x2);
    %    d2 = sqrt(d1^2 - g1*g2);
    %    minPos = x2 - (x2 - x1)*((g2 + d2 - d1)/(g2 - g1 + 2*d2));
    %    t_new = min(max(minPos,x1),x2);
    [minVal minPos] = min(points(:,1));
    notMinPos = -minPos+3;
    d1 = points(minPos,3) + points(notMinPos,3) - 3*(points(minPos,2)-points(notMinPos,2))/(points(minPos,1)-points(notMinPos,1));
    d2 = sqrt(d1^2 - points(minPos,3)*points(notMinPos,3));
    if isreal(d2)
        t = points(notMinPos,1) - (points(notMinPos,1) - points(minPos,1))*((points(notMinPos,3) + d2 - d1)/(points(notMinPos,3) - points(minPos,3) + 2*d2));
        minPos = min(max(t,xminBound),xmaxBound);
    else
        minPos = (xmaxBound+xminBound)/2;
    end
    return;
end

% Constraints Based on available Function Values
A = zeros(0,order+1);
b = zeros(0,1);
for i = 1:nPoints
    if imag(points(i,2))==0
        constraint = zeros(1,order+1);
        for j = order:-1:0
            constraint(order-j+1) = points(i,1)^j;
        end
        A = [A;constraint];
        b = [b;points(i,2)];
    end
end

% Constraints based on available Derivatives
for i = 1:nPoints
    if isreal(points(i,3))
        constraint = zeros(1,order+1);
        for j = 1:order
            constraint(j) = (order-j+1)*points(i,1)^(order-j);
        end
        A = [A;constraint];
        b = [b;points(i,3)];
    end
end

% Find interpolating polynomial
[params,ignore] = linsolve(A,b);

% Compute Critical Points
dParams = zeros(order,1);
for i = 1:length(params)-1
    dParams(i) = params(i)*(order-i+1);
end

if any(isinf(dParams))
    cp = [xminBound;xmaxBound;points(:,1)].';
else
    cp = [xminBound;xmaxBound;points(:,1);roots(dParams)].';
end

% Test Critical Points
fmin = inf;
minPos = (xminBound+xmaxBound)/2; % Default to Bisection if no critical points valid
for xCP = cp
    if imag(xCP)==0 && xCP >= xminBound && xCP <= xmaxBound
        fCP = polyval(params,xCP);
        if imag(fCP)==0 && fCP < fmin
            minPos = real(xCP);
            fmin = real(fCP);
        end
    end
end

% Plot Situation
if doPlot
    clf; hold on;
    
    % Plot Points
    plot(points(:,1),points(:,2),'b*');
    
    % Plot Derivatives
    for i = 1:nPoints
        if isreal(points(i,3))
            m = points(i,3);
            b = points(i,2) - m*points(i,1);
            plot([points(i,1)-.05 points(i,1)+.05],...
                [(points(i,1)-.05)*m+b (points(i,1)+.05)*m+b],'c.-');
        end
    end
    
    % Plot Function
    x = min(xmin,xminBound)-.1:(max(xmax,xmaxBound)+.1-min(xmin,xminBound)+.1)/100:max(xmax,xmaxBound)+.1;
    for i = 1:length(x)
        f(i) = polyval(params,x(i));
    end
    plot(x,f,'y');
    axis([x(1)-.1 x(end)+.1 min(f)-.1 max(f)+.1]);
    
    % Plot Minimum
    plot(minPos,fmin,'g+');
    if doPlot == 1
        pause(1);
    end
end

function [t,x_new,f_new,g_new,funEvals] = ArmijoBacktrack(x,t,d,f,fr,g,gtd,c1,LS,tolX,debug,doPlot,funObj,varargin)
%
% Backtracking linesearch to satisfy Armijo condition
%
% Inputs:
%   x: starting location
%   t: initial step size
%   d: descent direction
%   f: function value at starting location
%   fr: reference function value (usually funObj(x))
%   gtd: directional derivative at starting location
%   c1: sufficient decrease parameter
%   debug: display debugging information
%   LS: type of interpolation
%   tolX: minimum allowable step length
%   doPlot: do a graphical display of interpolation
%   funObj: objective function
%   varargin: parameters of objective function
%
% Outputs:
%   t: step length
%   f_new: function value at x+t*d
%   g_new: gradient value at x+t*d
%   funEvals: number function evaluations performed by line search


% Evaluate the Objective and Gradient at the Initial Step
[f_new,g_new] = feval(funObj, x + t*d, varargin{:});
funEvals = 1;

while f_new > fr + c1*t*gtd || ~isLegal(f_new)
    
    temp = t;
    if LS == 0 || ~isLegal(f_new)
        % Backtrack w/ fixed backtracking rate
        if debug
            fprintf('Fixed BT\n');
        end
        t = 0.5*t;
    elseif LS == 2 && isLegal(g_new)
        % Backtracking w/ cubic interpolation w/ derivative
        if debug
            fprintf('Grad-Cubic BT\n');
        end
        t = polyinterp([0 f gtd; t f_new mdot(g_new,d)],doPlot);
    elseif funEvals < 2 || ~isLegal(f_prev)
        % Backtracking w/ quadratic interpolation (no derivative at new point)
        if debug
            fprintf('Quad BT\n');
        end
        t = polyinterp([0 f gtd; t f_new sqrt(-1)],doPlot);
    else%if LS == 1
        % Backtracking w/ cubic interpolation (no derivatives at new points)
        if debug
            fprintf('Cubic BT\n');
        end
        t = polyinterp([0 f gtd; t f_new sqrt(-1); t_prev f_prev sqrt(-1)],doPlot);
    end
    
    % Adjust if change in t is too small/large
    
    if t < temp*1e-3
        if debug
            fprintf('Interpolated Value Too Small, Adjusting\n');
        end
        t = temp*1e-3;
    elseif t > temp*0.6
        if debug
            fprintf('Interpolated Value Too Large, Adjusting\n');
        end
        t = temp*0.6;
    end
    
    f_prev = f_new;
    t_prev = temp;
    
    [f_new,g_new] = feval(funObj, x + t*d, varargin{:});
    
    funEvals = funEvals+1;
    
    % Check whether step size has become too small
    if sum(abs(t*d)) <= tolX
        if debug
            fprintf('Backtracking Line Search Failed\n');
        end
        t = 0;
        f_new = f;
        g_new = g;
        break;
    end
end


x_new = x + t*d;


function [r] = mdot(A,B)
r = sum(sum(A.*B));
