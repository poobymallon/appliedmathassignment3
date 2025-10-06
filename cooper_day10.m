function cooper_day10()
    [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01,[0,10],1,.35);

    figure
    hold on
    fplot(@solution01, [0,10]);
    plot(t_list, X_list)

    [t_list2,X_list2,h_avg2, num_evals2] = forward_euler_fixed_step_integration(@rate_func02,[0,10],[1;0],.1);
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

%Runs numerical integration using forward Euler approximation
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