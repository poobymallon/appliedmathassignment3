function cooper_day10()
    [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01,[0,10],1,.35);

    figure
    hold on
    fplot(@solution01, [0,10]);
    plot(t_list, X_list)

    [t_list2,X_list2,h_avg2, num_evals2] = explicit_midpoint_fixed_step_integration(@rate_func01,[0,10],1,.1);

    figure
    hold on
    fplot(@solution01, [0,10]);
    plot(t_list2, X_list2)
    num = 100;
    hs = logspace(-5, -1, num);
    loc = 0.492;
    X_real = solution01(loc);

    %local
    loc_for = [];
    loc_mid = [];
    for i = 1:num
        [XB,~] = forward_euler_step(@rate_func01,loc,1,hs(i));
        loc_for(end+1) = XB-X_real;
        [XB,~] = explicit_midpoint_step(@rate_func01,loc,1,hs(i));
        loc_mid(end+1) = XB-X_real;
    end

    [p_lfor,k_lfor] = loglog_fit(hs,loc_for);
    [p_lmid,k_lmid] = loglog_fit(hs,loc_mid);
    figure
    loglog(hs, loc_for)
    hold on
    loglog(hs, loc_mid)

    %global
    glob_for = [];
    glob_mid = [];

    for i = 1:num
        [~, X_list, ~, ~] = forward_euler_fixed_step_integration(@rate_func01,[0,loc],1,hs(i));
        glob_for(end+1) = X_list(end)-X_real;
        [~, X_list, ~, ~] = explicit_midpoint_fixed_step_integration(@rate_func01,[0,loc],1,hs(i));
        glob_mid(end+1) = X_list(end)-X_real;
    end

    [p_gfor,k_gfor] = loglog_fit(hs,glob_for);
    [p_gmid,k_gmid] = loglog_fit(hs,glob_mid);
    figure
    loglog(hs, glob_for)
    hold on
    loglog(hs, glob_mid)
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
    t_list = linspace(ti, tf, num_steps);
    XA = X0;
    X_list = [];
    num_evals = 0;
    for i = 1:num_steps
        X_list(:, end+1) = XA;
        [XB,add_evals] = forward_euler_step(rate_func_in,t_list(i),XA,h_avg);
        XA = XB;
        num_evals = num_evals+add_evals;
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
    t_list = linspace(ti, tf, num_steps);
    XA = X0;
    X_list = [];
    num_evals = 0;
    for i = 1:num_steps
        X_list(:, end+1) = XA;
        [XB,add_evals] = explicit_midpoint_step(rate_func_in,t_list(i),XA,h_avg);
        XA = XB;
        num_evals = num_evals+add_evals;
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