function cooper_day10()
    [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01,[0,10],1,.1);

    figure
    hold on
    fplot(@solution01, [0,10]);
    plot(t_list, X_list)

    [t_list2,X_list2,h_avg2, num_evals2] = explicit_midpoint_fixed_step_integration(@rate_func01,[0,10],1,.1);

    figure
    hold on
    fplot(@solution01, [0,10]);
    plot(t_list2, X_list2)




    [t_list3, X_list3, h_avg3, num_evals3] = fixed_step_integration(@rate_func01,@backward_euler_step,[0,10],1,.1);
    [t_list4, X_list4, h_avg4, num_evals4] = fixed_step_integration(@rate_func01,@implicit_midpoint_step,[0,10],1,.1);

    % 
    % 
    num = 100;
    hs = logspace(-5, -1, num);
    tref = 0.492;
    X_real = solution01(tref);

    %local
    analytical_diff = zeros(size(hs));
    loc_for = zeros(size(hs));
    loc_mid = zeros(size(hs));
    for i = 1:num
        X_next = solution01(tref+hs(i));
        analytical_diff(i) = norm(X_next-X_real);
        [XB,~] = forward_euler_step(@rate_func01,tref,X_real,hs(i));
        loc_for(i) = norm(XB-X_next);
        [XB,~] = explicit_midpoint_step(@rate_func01,tref,X_real,hs(i));
        loc_mid(i) = norm(XB-X_next);
    end

    [p_a,k_a] = loglog_fit(hs,analytical_diff);
    [p_lfor,k_lfor] = loglog_fit(hs,loc_for);
    [p_lmid,k_lmid] = loglog_fit(hs,loc_mid);
    p_a
    p_lfor
    p_lmid
    figure
    loglog(hs, analytical_diff,'mo','markerfacecolor','m','markersize',2);
    hold on
    loglog(hs, loc_for,'ro','markerfacecolor','r','markersize',2);
    loglog(hs, loc_mid,'bo','markerfacecolor','b','markersize',2);

    loglog(hs,k_a*hs.^p_a,'k','LineWidth',1);
    loglog(hs,k_lfor*hs.^p_lfor,'k','LineWidth',1);
    loglog(hs,k_lmid*hs.^p_lmid,'k','LineWidth',1);

    %global

    glob_for =  zeros(size(hs));
    glob_mid =  zeros(size(hs));
    h_real_for =  zeros(size(hs));
    h_real_mid =  zeros(size(hs));

    t0 = 0; tf = 2;
    X0 = solution01(t0);
    Xf = solution01(tf);
    for i = 1:num
        [~, X_list, hf, ~] = forward_euler_fixed_step_integration(@rate_func01,[t0,tf],X0,hs(i));
        glob_for(i) = norm(X_list(:,end)-Xf);
        h_real_for(i) = hf;
        [~, X_list, hm, ~] = explicit_midpoint_fixed_step_integration(@rate_func01,[t0,tf],X0,hs(i));
        glob_mid(i) = norm(X_list(:,end)-Xf);
        h_real_mid(i) = hm;
    end

    [p_gfor,k_gfor] = loglog_fit(h_real_for,glob_for);
    [p_gmid,k_gmid] = loglog_fit(h_real_mid,glob_mid);
    figure
    loglog(h_real_for, glob_for,'ro','markerfacecolor','r','markersize',2)
    hold on;
    % grid on;
    loglog(h_real_mid, glob_mid,'bo','markerfacecolor','b','markersize',2)
    % 
    y_fit_for = k_gfor * h_real_for.^p_gfor;
    y_fit_mid = k_gmid * h_real_mid.^p_gmid;
    % 
    % % Plot fitted lines
    loglog(h_real_for, y_fit_for, 'k-', 'LineWidth', 1.5)
    loglog(h_real_mid, y_fit_mid, 'k-', 'LineWidth', 1.5)

    p_gfor
    p_gmid
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dXdt = rate_func02(t,X)
    dXdt = [0,-1;1,0]*X;
end
function X = solution02(t)
    X = [cos(t);sin(t)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [XB,num_evals] = forward_euler_step(rate_func_in,t,XA,h)
    XB = XA + h*rate_func_in(t, XA);
    num_evals = 1;
end

function [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    ti = tspan(1);
    tf = tspan(2);
    num_steps = ceil((tf-ti)/h_ref);
    h_avg = (tf-ti)/num_steps;
    t_list = linspace(ti, tf, num_steps+1);
    XA = X0;
    X_list = zeros(length(X0),num_steps+1);
    num_evals = 0;
    X_list(:, 1) = XA;
    for i = 1:num_steps
        [XB,add_evals] = forward_euler_step(rate_func_in,t_list(i),XA,h_avg);
        XA = XB;
        num_evals = num_evals+add_evals;
        X_list(:, i+1) = XA;
    end
end

function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    Xhalf = XA + (h/2)*rate_func_in(t, XA);
    XB = XA + h*rate_func_in(t+h/2, Xhalf);
    num_evals = 2;
end

%Runs numerical integration using explicit midpoint approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = explicit_midpoint_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    ti = tspan(1);
    tf = tspan(2);
    num_steps = ceil((tf-ti)/h_ref);
    h_avg = (tf-ti)/num_steps;
    t_list = linspace(ti, tf, num_steps+1);
    XA = X0;
    X_list = zeros(length(X0),num_steps+1);
    num_evals = 0;
    X_list(:, 1) = XA;
    for i = 1:num_steps
        [XB,add_evals] = explicit_midpoint_step(rate_func_in,t_list(i),XA,h_avg);
        XA = XB;
        num_evals = num_evals+add_evals;
        X_list(:, i+1) = XA;
    end
end

function [XB,num_evals] = backward_euler_step(rate_func_in,t,XA,h)
    back_euler = @(xb) XA + h*rate_func_in(t+h, xb)-xb;
    [XB, ~, ~, ~, num_evals] = multi_newton_solver(back_euler, XA, 1e-14, 1e-14, 200, 1);
end

function [XB,num_evals] = implicit_midpoint_step(rate_func_in,t,XA,h)
    imp_mid = @(xb) XA + h*rate_func_in(t+h/2, .5*(XA+xb))-xb;
    [XB, ~, ~, ~, num_evals] = multi_newton_solver(imp_mid, XA, 1e-14, 1e-14, 200, 1);
end

function [t_list, X_list, h_avg, num_evals] = fixed_step_integration(rate_func_in,step_func,tspan,X0,h_ref)
    ti = tspan(1);
    tf = tspan(2);
    num_steps = ceil((tf-ti)/h_ref);
    h_avg = (tf-ti)/num_steps;
    t_list = linspace(ti, tf, num_steps);
    XA = X0;
    X_list = [];
    num_evals = 0;
    for i = 1:num_steps
        X_list(:, end+1) = XA;
        [XB,add_evals] = step_func(rate_func_in,t_list(i),XA,h_avg);
        XA = XB;
        num_evals = num_evals+add_evals;
    end
end




function [root, it, flag, glist, num_evals] = multi_newton_solver(FUN_both, X0, Athresh, Bthresh, maxit, numdiff)
    % 1 = numerical, 2 = analytical
    % fast convergence near     root
    % using slope to jump closer
    num_evals = 0;
    glist = [];                             % step 1: start xn list
    root = X0; it = 0; flag = 0;
    glist(:, end+1) = root;                    % step 1: save first guess

    if numdiff == 2
        try
            [FX, J] = FUN_both(root);  % Try to get both outputs (analytical)
            num_evals = num_evals+1;
        catch
            FX = FUN_both(root);       % Fallback if FUN_both only returns one output
            [J, n] = approximate_jacobian(FUN_both, root);
            num_evals = num_evals+n+1;
        end
    else
        FX = FUN_both(root);           % Numerical mode: always get 1 output
        [J, n] = approximate_jacobian(FUN_both, root);
        num_evals = num_evals+n+1;
    end

    if numdiff == 1
        [J, n] = approximate_jacobian(FUN_both, X0);
        num_evals = num_evals+n;
    end
    while norm(FX) > Bthresh && it < maxit
        if det(J*J') == 0
            flag = -2; return               % -2 = zero derivative
        end
        X_new = root - J\FX;              % updating step
        glist(:, end+1) = X_new;                % step 1: save iterate (kept as is)
        if abs(X_new - root) < Athresh
            root = X_new; flag = 1; return
        end
        root = X_new;
        if numdiff == 2
            try
                [FX, J] = FUN_both(root);  % Try to get both outputs (analytical)
                num_evals = num_evals+1;
            catch
                FX = FUN_both(root);       % Fallback if FUN_both only returns one output
                [J, n] = approximate_jacobian(FUN_both, root);
                num_evals = num_evals+n+1;
            end
        else
            FX = FUN_both(root);           % Numerical mode: always get 1 output
            [J, n] = approximate_jacobian(FUN_both, root);
            num_evals = num_evals+n+1;
        end
        it = it + 1;
    end
    flag = 1;
end

function [J, num_evals] = approximate_jacobian(FUN,X)
    num_evals = 0;
    h = 1e-6;
    J = [];
    varnum = length(X);
    for i = 1:varnum
        ei = zeros(varnum, 1);
        ei(i) = h/2;
        J(:,end+1) = (FUN(X+ei)-FUN(X-ei))./h;
        num_evals = num_evals+2;
    end
end

function [p,k] = loglog_fit(x_regression,y_regression,varargin)
%convert x_regression to a column vector if it's a row vector
if size(x_regression,1)==1
x_regression = abs(x_regression)';
end
%convert y_regression to a column vector if it's a row vector
if size(y_regression,1)==1
y_regression = abs(y_regression)';
end
%if filter_params has been provided, then filter the data points
if nargin==3
filter_params = varargin{1};
num_points = length(x_regression);
indices = 1:num_points;
filter_bool = ones(num_points,1)>=1;
if isfield(filter_params,'min_index')
filter_bool = filter_bool & indices>=filter_params.min_index;
end
if isfield(filter_params,'max_index')
filter_bool = filter_bool & indices<=filter_params.max_index;
end
if isfield(filter_params,'min_xval')
filter_bool = filter_bool & x_regression>=filter_params.min_xval;
end
if isfield(filter_params,'max_xval')
filter_bool = filter_bool & x_regression<=filter_params.max_xval;
end
if isfield(filter_params,'min_yval')
filter_bool = filter_bool & y_regression>=filter_params.min_yval;
end
if isfield(filter_params,'max_yval')
filter_bool = filter_bool & y_regression<=filter_params.max_yval;
end
x_regression = x_regression(filter_bool);
y_regression = y_regression(filter_bool);
end
%compute the logs of x_regression and y_regression
Y = log(y_regression);
X1 = log(x_regression);
%set up the regression
X2 = ones(length(X1),1);
%run the regression
coeff_vec = regress(Y,[X1,X2]);
%pull out the coefficients from the fit
p = coeff_vec(1);
k = exp(coeff_vec(2));
end