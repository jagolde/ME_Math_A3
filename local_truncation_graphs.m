function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end
function X = solution01(t)
X = cos(t);
end
clear;
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
close all


new_rate= @(t, X) rate_func01(t, X);
G_func1 = @(t, X0, h_ref) forward_euler_step(new_rate,t,X0,h_ref);
G_func2 = @(t, X0, h_ref) explicit_midpoint_step(new_rate,t,X0,h_ref);
G_func3 = @(t, X0, h_ref) implicit_midpoint_step(new_rate,t,X0,h_ref);
G_func4 = @(t, X0, h_ref) backward_euler_step(new_rate,t,X0,h_ref);

sol3= @(t) solution01(t);

h_vals = logspace(-5, 1, 100);
t_ref = 0.492;
e_list1= [];
e_list2= [];
e_list3= [];
e_list4= [];
for h = h_vals
    epsilon1 = local_truncation_error(G_func1, sol3, t_ref, h);
    e_list1(end+1) = epsilon1;
    epsilon2 = local_truncation_error(G_func2, sol3, t_ref, h);
    e_list2(end+1) = epsilon2;
    epsilon3 = local_truncation_error(G_func3, sol3, t_ref, h);
    e_list3(end+1) = epsilon3;
    epsilon4 = local_truncation_error(G_func4, sol3, t_ref, h);
    e_list4(end+1) = epsilon4;
end
errors1 = abs(e_list1);
errors2 = abs(e_list2);
errors3 = abs(e_list3);
errors4 = abs(e_list4);
loglog(h_vals, errors1, 'r--',displayname='Forward Euler', LineWidth=1.5); hold on
loglog(h_vals, errors2, 'k', displayname='Explicit Midpoint', LineWidth=1.5);
loglog(h_vals, errors3, 'g', displayname='Implicit Midpoint', LineWidth=1.5);
loglog(h_vals, errors4, 'b--',displayname='Backward Euler', LineWidth=1.5);


%fit lines
% errors = abs(e_list1);
% valid = (errors > 0) & isfinite(errors);
% 
% p = polyfit(log10(h_vals(valid)), log10(errors(valid)), 1);
% fit_errors = 10.^polyval(p, log10(h_vals));
% 
% 
% figure;
% loglog(h_vals, errors, 'r',displayname='Forward Euler'); hold on
% loglog(h_vals, fit_errors, 'b--', displayname=sprintf('Forward Euler Fit, slope = %.2f', p(1)));
% 
% 
% 
% errors = abs(e_list2);
% valid = (errors > 0) & isfinite(errors);
% 
% p = polyfit(log10(h_vals(valid)), log10(errors(valid)), 1);
% fit_errors = 10.^polyval(p, log10(h_vals));
% 
% loglog(h_vals, errors, 'k', displayname='Explicit Midpoint');
% loglog(h_vals, fit_errors, 'g--', displayname = sprintf('Explicit Midpoint Fit, slope = %.2f', p(1)));










legend('location', 'best')
xlabel('Step Size (~)')
ylabel('Measured Local Truncation Error $\epsilon$ (~)')
title('Comparison of Four Different Step Methods for (-5x + 5cos(t) - sin(t))')