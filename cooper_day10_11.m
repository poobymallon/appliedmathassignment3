function cooper_day10_11()
    % part: menu to generate plots in guideline order (1..10)

    clc; close all;
    disp("choose which plots to generate (numbers refer to the submission order):");
    disp(" ");
    opts = [
        "1. fig 1: forward euler vs exact (several h, overlay)"
        "2. fig 2: explicit midpoint vs exact (several h, overlay)"
        "3. fig 3: solutions vs exact at h_ref=0.38 (4 subplots)"
        "4. fig 4: solutions vs exact at h_ref=0.45 (4 subplots)"
        "5. fig 5: local truncation (explicit only, with |X(t+h)-X(t)| and fit lines)"
        "6. fig 6: local truncation (all methods, no fit lines)"
        "7. fig 7: global truncation vs h (explicit only, with fit lines)"
        "8. fig 8: global truncation vs h (all methods, no fit lines)"
        "9. fig 9: global truncation vs # of f calls (explicit only, with fit lines)"
        "10. fig 10: global truncation vs # of f calls (all methods, no fit lines)"
    ];
    disp(join(opts, newline)); disp(" ");
    choice = input("enter plot number:");

    for i = reshape(choice,1,[])
        switch i
            case 1,  fig1_forward_euler_vs_exact();
            case 2,  fig2_explicit_midpoint_vs_exact();
            case 3,  fig3_stability(0.38);
            case 4,  fig4_stability(0.45);
            case 5,  fig5_local_explicit_with_fits();
            case 6,  fig6_local_all_no_fits();
            case 7,  fig7_global_explicit_with_fits();
            case 8,  fig8_global_all_no_fits();
            case 9,  fig9_global_calls_explicit_with_fits();
            case 10, fig10_global_calls_all_no_fits();
            otherwise, warning("invalid choice: %d", i);
        end
    end
    if isempty(choice), disp("no plots selected."); else, disp("done"); end
end

% FIGURES IN ORDSER OF SUBMISSIONS GUIDELINES 

% fig 1: forward euler vs exact (overlay, several h)
function fig1_forward_euler_vs_exact()
    % uses stable set of step sizes
    t0 = 0; tf = 10; X0 = 1;
    hs_show = [0.4 0.2 0.1];

    figure(1); clf; hold on; grid on; box on
    fplot(@solution01,[t0,tf],'k','LineWidth',1.8)
    colors = lines(numel(hs_show)); markers = {'o','s','^'};
    for i = 1:numel(hs_show)
        h = hs_show(i);
        [tF, XF] = forward_euler_fixed_step_integration(@rate_func01,[t0,tf],X0,h);
        plot(tF, XF, ['-' markers{i}], 'Color', colors(i,:), 'LineWidth', 1.0, 'MarkerSize', 4)
    end
    xlabel('t'); ylabel('x')
    title('fig 1: forward euler vs exact (overlay)')
    legend(['exact', compose('h=%.2g',hs_show)],'location','best')
end

% fig 2: explicit midpoint vs exact (overlay, several h)
function fig2_explicit_midpoint_vs_exact()
    t0 = 0; tf = 10; X0 = 1;
    hs_show = [0.4 0.2 0.1];

    figure(2); clf; hold on; grid on; box on
    fplot(@solution01,[t0,tf],'k','LineWidth',1.8)
    colors = lines(numel(hs_show)); markers = {'o','s','^'};
    for i = 1:numel(hs_show)
        h = hs_show(i);
        [tM, XM] = explicit_midpoint_fixed_step_integration(@rate_func01,[t0,tf],X0,h);
        plot(tM, XM, ['-' markers{i}], 'Color', colors(i,:), 'LineWidth', 1.0, 'MarkerSize', 4)
    end
    xlabel('t'); ylabel('x')
    title('fig 2: explicit midpoint vs exact (overlay)')
    legend(['exact', compose('h=%.2g',hs_show)],'location','best')
end

% fig 3: stability-style comparison at href=0.38 (4 subplots)
function fig3_stability(h_ref)
    t0 = 0; tf = 20; X0 = solution01(t0);

    [tF, XF] = forward_euler_fixed_step_integration(@rate_func01,[t0,tf],X0,h_ref);
    [tB, XB] = backward_euler_fixed_step_integration(@rate_func01,[t0,tf],X0,h_ref);
    [tE, XE] = explicit_midpoint_fixed_step_integration(@rate_func01,[t0,tf],X0,h_ref);
    [tI, XI] = implicit_midpoint_fixed_step_integration(@rate_func01,[t0,tf],X0,h_ref);

    figure(3); clf; tiledlayout(4,1,'Padding','compact','TileSpacing','compact')
    nexttile; hold on; grid on; box on; fplot(@solution01,[t0,tf]); plot(tF, XF,'o-')
    title(sprintf('forward euler, h=%.3f', h_ref)); xlabel('t'); ylabel('x'); legend('exact','num','location','best')
    nexttile; hold on; grid on; box on; fplot(@solution01,[t0,tf]); plot(tB, XB,'o-')
    title(sprintf('backward euler, h=%.3f', h_ref)); xlabel('t'); ylabel('x'); legend('exact','num','location','best')
    nexttile; hold on; grid on; box on; fplot(@solution01,[t0,tf]); plot(tE, XE,'o-')
    title(sprintf('explicit midpoint, h=%.3f', h_ref)); xlabel('t'); ylabel('x'); legend('exact','num','location','best')
    nexttile; hold on; grid on; box on; fplot(@solution01,[t0,tf]); plot(tI, XI,'o-')
    title(sprintf('implicit midpoint, h=%.3f', h_ref)); xlabel('t'); ylabel('x'); legend('exact','num','location','best')
    sgtitle('fig 3: numerical solutions vs exact, h_{ref}=0.38')
end

% fig 4: stability-style comparison at href=0.45 (4 subplots)
function fig4_stability(h_ref)
    t0 = 0; tf = 20; X0 = solution01(t0);

    [tF, XF] = forward_euler_fixed_step_integration(@rate_func01,[t0,tf],X0,h_ref);
    [tB, XB] = backward_euler_fixed_step_integration(@rate_func01,[t0,tf],X0,h_ref);
    [tE, XE] = explicit_midpoint_fixed_step_integration(@rate_func01,[t0,tf],X0,h_ref);
    [tI, XI] = implicit_midpoint_fixed_step_integration(@rate_func01,[t0,tf],X0,h_ref);

    figure(4); clf; tiledlayout(4,1,'Padding','compact','TileSpacing','compact')
    nexttile; hold on; grid on; box on; fplot(@solution01,[t0,tf]); plot(tF, XF,'o-')
    title(sprintf('forward euler, h=%.3f', h_ref)); xlabel('t'); ylabel('x'); legend('exact','num','location','best')
    nexttile; hold on; grid on; box on; fplot(@solution01,[t0,tf]); plot(tB, XB,'o-')
    title(sprintf('backward euler, h=%.3f', h_ref)); xlabel('t'); ylabel('x'); legend('exact','num','location','best')
    nexttile; hold on; grid on; box on; fplot(@solution01,[t0,tf]); plot(tE, XE,'o-')
    title(sprintf('explicit midpoint, h=%.3f', h_ref)); xlabel('t'); ylabel('x'); legend('exact','num','location','best')
    nexttile; hold on; grid on; box on; fplot(@solution01,[t0,tf]); plot(tI, XI,'o-')
    title(sprintf('implicit midpoint, h=%.3f', h_ref)); xlabel('t'); ylabel('x'); legend('exact','num','location','best')
    sgtitle('fig 4: numerical solutions vs exact, h_{ref}=0.45')
end

% fig 5: local truncation scaling for explicit methods (with |X(t+h)-X(t)| and fit lines)
function fig5_local_explicit_with_fits()
    num = 100; hs = logspace(-5,-1,num); tref = 0.492; Xref = solution02(tref); %Xref = solution01(tref);
    analytical = zeros(size(hs)); loc_fe = analytical; loc_em = analytical;
    for i=1:num
        h = hs(i); Xtrue = solution02(tref+h); %Xtrue = solution01(tref+h);
        analytical(i) = norm(Xtrue - Xref);
        [XB,~] = forward_euler_step(@rate_func02, tref, Xref, h);      loc_fe(i)=norm(XB-Xtrue);
        [XB,~] = explicit_midpoint_step(@rate_func02, tref, Xref, h);  loc_em(i)=norm(XB-Xtrue);
    end
    [p_ref,k_ref] = loglog_fit(hs, analytical);
    [p_fe,k_fe]   = loglog_fit(hs, loc_fe);
    [p_em,k_em]   = loglog_fit(hs, loc_em);
    fprintf('[local explicit fits] ref p≈%.3f | fe p≈%.3f | em p≈%.3f\n', p_ref, p_fe, p_em);

    figure(5); clf;  grid on; box on
    loglog(hs, analytical,'mo','markerfacecolor','m','markersize',3)
    hold on;
    loglog(hs, loc_fe,'ro','markerfacecolor','r','markersize',3)
    loglog(hs, loc_em,'bo','markerfacecolor','b','markersize',3)
    loglog(hs, k_ref*hs.^p_ref,'k-','linewidth',1)
    loglog(hs, k_fe*hs.^p_fe,'k-','linewidth',1)
    loglog(hs, k_em*hs.^p_em,'k-','linewidth',1)
    xlabel('h'); ylabel('local error')
    title('fig 5: local truncation (explicit only) with |X(t+h)-X(t)| and fits')
    legend('|X(t+h)-X(t)|','forward euler','explicit midpoint','fits','location','southeast')
end

% fig 6: local truncation scaling for all four methods (no fit lines)
function fig6_local_all_no_fits()
    num = 60; hs = logspace(-4,-1,num); tref = 0.492; Xref = solution02(tref); %Xref = solution01(tref);
    loc_fe=zeros(size(hs)); loc_be=loc_fe; loc_em=loc_fe; loc_im=loc_fe;
    for i=1:num
        h = hs(i); Xtrue = solution02(tref+h); %Xtrue = solution01(tref+h);
        [xB,~] = forward_euler_step(@rate_func02, tref, Xref, h);      loc_fe(i)=norm(xB-Xtrue);
        [xB,~] = backward_euler_step(@rate_func02, tref, Xref, h);     loc_be(i)=norm(xB-Xtrue);
        [xB,~] = explicit_midpoint_step(@rate_func02, tref, Xref, h);  loc_em(i)=norm(xB-Xtrue);
        [xB,~] = implicit_midpoint_step(@rate_func02, tref, Xref, h);  loc_im(i)=norm(xB-Xtrue);
    end
    [p_be,k_be] = loglog_fit(hs, loc_be)
    [p_im,k_im] = loglog_fit(hs, loc_im)
    figure(6); clf; grid on; box on
    loglog(hs, loc_fe,'o','markerfacecolor','auto','markersize',3)
    hold on; 
    loglog(hs, loc_be,'o','markerfacecolor','auto','markersize',3)
    loglog(hs, loc_em,'o','markerfacecolor','auto','markersize',3)
    loglog(hs, loc_im,'o','markerfacecolor','auto','markersize',3)
    xlabel('h'); ylabel('local error')
    title('fig 6: local truncation (all methods, no fit lines)')
    legend('forward euler','backward euler','explicit midpoint','implicit midpoint','location','southeast')
end

% fig 7: global truncation vs average h for explicit methods (with fits)
function fig7_global_explicit_with_fits()
    t0=0; tf=2; X0 = solution01(t0); Xf = solution01(tf);
    num=80; hs=logspace(-4,-1,num);
    glob_fe=zeros(size(hs)); glob_em=glob_fe; h_fe=glob_fe; h_em=glob_fe;
    for i=1:num
        h=hs(i);
        [~, XF, hf] = forward_euler_fixed_step_integration(@rate_func01,[t0,tf],X0,h);
        [~, XM, hm] = explicit_midpoint_fixed_step_integration(@rate_func01,[t0,tf],X0,h);
        glob_fe(i)=norm(XF(:,end)-Xf); h_fe(i)=hf;
        glob_em(i)=norm(XM(:,end)-Xf); h_em(i)=hm;
    end
    [p_gfe,k_gfe] = loglog_fit(h_fe, glob_fe);
    [p_gem,k_gem] = loglog_fit(h_em, glob_em);
    fprintf('[global vs h explicit fits] fe p≈%.3f | em p≈%.3f\n', p_gfe, p_gem);

    figure(7); clf;  grid on; box on
    loglog(h_fe, glob_fe,'ro','markerfacecolor','r','markersize',3)
    hold on;
    loglog(h_em, glob_em,'bo','markerfacecolor','b','markersize',3)
    loglog(h_fe, k_gfe*h_fe.^p_gfe,'k-','linewidth',1.25)
    loglog(h_em, k_gem*h_em.^p_gem,'k-','linewidth',1.25)
    xlabel('average h'); ylabel('global error at t_f')
    title('fig 7: global truncation vs h (explicit only, with fits)')
    legend('forward euler','explicit midpoint','fits','location','southeast')
end

% fig 8: global truncation vs average h for all methods (no fit lines)
function fig8_global_all_no_fits()
    t0=0; tf=2; X0 = solution02(t0); Xf = solution02(tf); %X0 = solution01(t0); Xf = solution01(tf);
    num=60; hs=logspace(-4,-1,num);
    glob_fe=zeros(size(hs)); glob_be=glob_fe; glob_em=glob_fe; glob_im=glob_fe;
    h_fe=glob_fe; h_be=glob_fe; h_em=glob_fe; h_im=glob_fe;
    for i=1:num
        h=hs(i);
        [~, XF, hf] = forward_euler_fixed_step_integration(@rate_func02,[t0,tf],X0,h);
        [~, XB, hb] = backward_euler_fixed_step_integration(@rate_func02,[t0,tf],X0,h);
        [~, XM, hm] = explicit_midpoint_fixed_step_integration(@rate_func02,[t0,tf],X0,h);
        [~, XI, hi] = implicit_midpoint_fixed_step_integration(@rate_func02,[t0,tf],X0,h);
        glob_fe(i)=norm(XF(:,end)-Xf); h_fe(i)=hf;
        glob_be(i)=norm(XB(:,end)-Xf); h_be(i)=hb;
        glob_em(i)=norm(XM(:,end)-Xf); h_em(i)=hm;
        glob_im(i)=norm(XI(:,end)-Xf); h_im(i)=hi;
    end
    [p_gfe,k_gfe] = loglog_fit(h_fe, glob_fe)
    [p_gbe,k_gbe] = loglog_fit(h_fe, glob_be)
    [p_gem,k_gem] = loglog_fit(h_fe, glob_em)
    [p_gim,k_gim] = loglog_fit(h_fe, glob_im)

    figure(8); clf;  grid on; box on
    loglog(h_fe, glob_fe,'o','markerfacecolor','auto','markersize',3)
    hold on;
    loglog(h_be, glob_be,'o','markerfacecolor','auto','markersize',3)
    loglog(h_em, glob_em,'o','markerfacecolor','auto','markersize',3)
    loglog(h_im, glob_im,'o','markerfacecolor','auto','markersize',3)
    xlabel('average h'); ylabel('global error at t_f')
    title('fig 8: global truncation vs h (all methods, no fit lines)')
    legend('forward euler','backward euler','explicit midpoint','implicit midpoint','location','southeast')
end

% fig 9: global truncation vs number of f calls for explicit methods (with fits)
function fig9_global_calls_explicit_with_fits()
    t0=0; tf=2; X0=solution01(t0); Xf=solution01(tf);
    num=60; hs=logspace(-4,-1,num);
    calls_fe=zeros(size(hs)); calls_em=calls_fe; glob_fe=calls_fe; glob_em=calls_fe;
    for i=1:num
        h=hs(i);
        [~, XF, ~, nfe] = forward_euler_fixed_step_integration(@rate_func01,[t0,tf],X0,h);
        [~, XM, ~, nem] = explicit_midpoint_fixed_step_integration(@rate_func01,[t0,tf],X0,h);
        calls_fe(i)=nfe; glob_fe(i)=norm(XF(:,end)-Xf);
        calls_em(i)=nem; glob_em(i)=norm(XM(:,end)-Xf);
    end
    [p_fe,k_fe] = loglog_fit(calls_fe, glob_fe);
    [p_em,k_em] = loglog_fit(calls_em, glob_em);
    fprintf('[global vs calls explicit fits] fe p≈%.3f | em p≈%.3f\n', p_fe, p_em);

    figure(9); clf;  grid on; box on
    loglog(calls_fe, glob_fe,'ro','markerfacecolor','r','markersize',3)
    hold on;
    loglog(calls_em, glob_em,'bo','markerfacecolor','b','markersize',3)
    loglog(calls_fe, k_fe*calls_fe.^p_fe,'k-','linewidth',1.25)
    loglog(calls_em, k_em*calls_em.^p_em,'k-','linewidth',1.25)
    xlabel('# of f calls'); ylabel('global error at t_f')
    title('fig 9: global truncation vs number of f calls (explicit only, with fits)')
    legend('forward euler','explicit midpoint','fits','location','southwest')
end

% fig 10: global truncation vs number of f calls for all methods (no fit lines)
function fig10_global_calls_all_no_fits()
    t0=0; tf=2; X0=solution02(t0); Xf=solution02(tf); %X0=solution01(t0); Xf=solution01(tf);
    num=60; hs=logspace(-4,-1,num);
    calls_fe=zeros(size(hs)); calls_be=calls_fe; calls_em=calls_fe; calls_im=calls_fe;
    glob_fe=calls_fe; glob_be=calls_fe; glob_em=calls_fe; glob_im=calls_fe;
    for i=1:num
        h=hs(i);
        [~, XF, ~, nfe] = forward_euler_fixed_step_integration(@rate_func02,[t0,tf],X0,h);
        [~, XB, ~, nbe] = backward_euler_fixed_step_integration(@rate_func02,[t0,tf],X0,h);
        [~, XM, ~, nem] = explicit_midpoint_fixed_step_integration(@rate_func02,[t0,tf],X0,h);
        [~, XI, ~, nim] = implicit_midpoint_fixed_step_integration(@rate_func02,[t0,tf],X0,h);
        calls_fe(i)=nfe; glob_fe(i)=norm(XF(:,end)-Xf);
        calls_be(i)=nbe; glob_be(i)=norm(XB(:,end)-Xf);
        calls_em(i)=nem; glob_em(i)=norm(XM(:,end)-Xf);
        calls_im(i)=nim; glob_im(i)=norm(XI(:,end)-Xf);
    end
    [p_gfe,k_gfe] = loglog_fit(calls_fe, glob_fe)
    [p_gbe,k_gbe] = loglog_fit(calls_be, glob_be)
    [p_gem,k_gem] = loglog_fit(calls_em, glob_em)
    [p_gim,k_gim] = loglog_fit(calls_im, glob_im)

    figure(10); clf;  grid on; box on
    loglog(calls_fe, glob_fe,'o','markerfacecolor','auto','markersize',3)
    hold on;
    loglog(calls_be, glob_be,'o','markerfacecolor','auto','markersize',3)
    loglog(calls_em, glob_em,'o','markerfacecolor','auto','markersize',3)
    loglog(calls_im, glob_im,'o','markerfacecolor','auto','markersize',3)
    xlabel('# of f calls'); ylabel('global error at t_f')
    title('fig 10: global truncation vs number of f calls (all methods, no fit lines)')
    legend('forward euler','backward euler','explicit midpoint','implicit midpoint','location','southwest')
end

% HELPERS AND METHODS

% rate and exact solution for eqn (5)
function dXdt = rate_func01(t,X), dXdt = -5*X + 5*cos(t) - sin(t); end
function X = solution01(t), X = cos(t); end

% rate and exact solution for second helper
function dXdt = rate_func02(t,X)
    dXdt = [0,-1;1,0]*X;
end
function X = solution02(t)
    X = [cos(t);sin(t)];
end

% forward euler step
function [XB,num_evals] = forward_euler_step(rate_func_in,t,XA,h)
    XB = XA + h*rate_func_in(t, XA); num_evals = 1;
end

% explicit midpoint step
function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    k1 = rate_func_in(t, XA);
    xh = XA + (h/2)*k1;
    k2 = rate_func_in(t + h/2, xh);
    XB = XA + h*k2; num_evals = 2;
end

% backward euler step (implicit)
function [XB,num_evals] = backward_euler_step(rate_func_in,t,XA,h)
    G = @(xb) XA + h*rate_func_in(t + h, xb) - xb;
    [XB, ~, ~, ~, num_evals] = multi_newton_solver(G, XA, 1e-12, 1e-12, 200, 1);
end

% implicit midpoint step
function [XB,num_evals] = implicit_midpoint_step(rate_func_in,t,XA,h)
    G = @(xb) XA + h*rate_func_in(t + h/2, 0.5*(XA + xb)) - xb;
    [XB, ~, ~, ~, num_evals] = multi_newton_solver(G, XA, 1e-12, 1e-12, 200, 1);
end

% shared fixed-step integrator
function [t_list, X_list, h_avg, num_evals] = fixed_step_integration(rate_func_in, step_func, tspan, X0, h_ref)
    ti = tspan(1); tf = tspan(2);
    N  = ceil((tf - ti)/h_ref); h_avg = (tf - ti)/N;
    t_list = linspace(ti, tf, N+1);
    nx = numel(X0); X_list = zeros(nx, N+1); X_list(:,1) = X0;
    num_evals = 0; XA = X0;
    for k = 1:N
        [XB, add_evals] = step_func(rate_func_in, t_list(k), XA, h_avg);
        num_evals = num_evals + add_evals;
        X_list(:,k+1) = XB; XA = XB;
    end
end

% wrappers
function [t_list,X_list,h_avg,num_evals] = forward_euler_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    [t_list,X_list,h_avg,num_evals] = fixed_step_integration(rate_func_in,@forward_euler_step,tspan,X0,h_ref);
end
function [t_list,X_list,h_avg,num_evals] = explicit_midpoint_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    [t_list,X_list,h_avg,num_evals] = fixed_step_integration(rate_func_in,@explicit_midpoint_step,tspan,X0,h_ref);
end
function [t_list,X_list,h_avg,num_evals] = backward_euler_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    [t_list,X_list,h_avg,num_evals] = fixed_step_integration(rate_func_in,@backward_euler_step,tspan,X0,h_ref);
end
function [t_list,X_list,h_avg,num_evals] = implicit_midpoint_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    [t_list,X_list,h_avg,num_evals] = fixed_step_integration(rate_func_in,@implicit_midpoint_step,tspan,X0,h_ref);
end

% newton solver + jacobian
function [root, it, flag, glist, num_evals] = multi_newton_solver(FUN, X0, Athresh, Bthresh, maxit, numdiff)
    num_evals = 0; it = 0; flag = 0; root = X0; glist = root;
    FX = FUN(root); num_evals = num_evals + 1;
    [J, n] = approximate_jacobian(FUN, root); num_evals = num_evals + n;
    while norm(FX) > Bthresh && it < maxit
        if rcond(J) < 1e-14, flag = -2; return; end
        X_new = root - J\FX; glist(:,end+1) = X_new; %#ok<AGROW>
        if norm(X_new - root) < Athresh, root = X_new; flag = 1; return; end
        root = X_new; FX = FUN(root); num_evals = num_evals + 1;
        [J, n] = approximate_jacobian(FUN, root); num_evals = num_evals + n; it = it + 1;
    end
    flag = (norm(FX) <= Bthresh);
end
function [J, num_evals] = approximate_jacobian(FUN, X)
    h = 1e-6; f0 = FUN(X); n = numel(X); m = numel(f0);
    J = zeros(m, n); num_evals = 1;
    for i = 1:n
        ei = zeros(n,1); ei(i) = h/2;
        fp = FUN(X + ei); fm = FUN(X - ei);
        J(:, i) = (fp - fm) / h; num_evals = num_evals + 2;
    end
end

% log-log linear fit helper
function [p,k] = loglog_fit(x_regression,y_regression,varargin)
    if size(x_regression,1)==1, x_regression = abs(x_regression)'; end
    if size(y_regression,1)==1, y_regression = abs(y_regression)'; end
    if nargin==3
        fp = varargin{1}; N = length(x_regression); idx = (1:N).';
        mask = true(N,1);
        if isfield(fp,'min_index'), mask = mask & idx>=fp.min_index; end
        if isfield(fp,'max_index'), mask = mask & idx<=fp.max_index; end
        if isfield(fp,'min_xval'),  mask = mask & x_regression>=fp.min_xval; end
        if isfield(fp,'max_xval'),  mask = mask & x_regression<=fp.max_xval; end
        if isfield(fp,'min_yval'),  mask = mask & y_regression>=fp.min_yval; end
        if isfield(fp,'max_yval'),  mask = mask & y_regression<=fp.max_yval; end
        x_regression = x_regression(mask); y_regression = y_regression(mask);
    end
    Y = log(y_regression); X = [log(x_regression), ones(length(x_regression),1)];
    coeff = regress(Y, X); p = coeff(1); k = exp(coeff(2));
end
