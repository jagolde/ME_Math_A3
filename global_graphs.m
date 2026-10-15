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


G_func1 = @(rate_func, t, X, h) forward_euler_step(rate_func, t, X, h);
G_func2 =@(rate_func, t, X, h) explicit_midpoint_step(rate_func, t, X, h);
G_func3=@(rate_func, t, X, h) implicit_midpoint_step(rate_func, t, X, h);
G_func4=@(rate_func, t, X, h) backward_euler_step(rate_func, t, X, h);


new_rate  = @(t, X) rate_func01(t, X);
sol3      = @(t) solution01(t);

h_vals = logspace(-5, -1, 100);
t_ref  = [0, 10];

e_list1 = zeros(size(h_vals));
e_list2 = zeros(size(h_vals));
e_list3 = zeros(size(h_vals));
e_list4 = zeros(size(h_vals));

for i = 1:numel(h_vals)
    e_list1(i) = global_truncation_error(new_rate, G_func1, sol3, t_ref, h_vals(i));
    e_list2(i) = global_truncation_error(new_rate, G_func2, sol3, t_ref, h_vals(i));
    e_list3(i) = global_truncation_error(new_rate, G_func3, sol3, t_ref, h_vals(i));
    e_list4(i) = global_truncation_error(new_rate, G_func4, sol3, t_ref, h_vals(i));
end

errors1 = abs(e_list1);
errors2 = abs(e_list2);
errors3 = abs(e_list3);
errors4 = abs(e_list4);
loglog(h_vals, errors1, 'r-','MarkerSize',2, displayname='Forward Euler'); hold on
loglog(h_vals, errors2, 'mo','MarkerSize',2, displayname='Explicit Midpoint');
loglog(h_vals, errors3, 'b-','MarkerSize',2, displayname='Implicit Midpoint');
loglog(h_vals, errors4, 'go','MarkerSize',2,displayname='Backward Euler');

%loglog(h_vals, fit_errors, 'k','LineWidth',1.5);
xlabel('H Values (-)')
ylabel('Epsilon Values (-)')
% title('Log Plot of Rate Func 1')
title('Global Truncation Error vs Step Size Experiment: All Four Methods')
legend()


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