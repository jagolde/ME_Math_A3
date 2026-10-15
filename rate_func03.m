function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end
function X = solution01(t)
X = cos(t);
end
close all

step_func = @(rate_func, t, X, h) forward_euler_step(rate_func, t, X, h);
step_func2 = @(rate_func, t, X, h) explicit_midpoint_step(rate_func, t, X, h);
new_rate  = @(t, X) rate_func01(t, X);
sol3      = @(t) solution01(t);

h_vals = logspace(-5, -1, 100);
t_ref  = [0, 10];
e_list = zeros(size(h_vals));
e_list2 = zeros(size(h_vals));

for i = 1:numel(h_vals)
    e_list(i) = global_truncation_error(new_rate, step_func, sol3, t_ref, h_vals(i));
    e_list2(i) = global_truncation_error(new_rate, step_func2, sol3, t_ref, h_vals(i));
end

errors = abs(e_list);
errors2 = abs(e_list2);
valid  = (errors > 0) & isfinite(errors);
valid2  = (errors2 > 0) & isfinite(errors2);

p = polyfit(log10(h_vals(valid)), log10(errors(valid)), 1);
fit_errors = 10.^polyval(p, log10(h_vals));

p2 = polyfit(log10(h_vals(valid2)), log10(errors(valid2)), 1);
fit_errors2 = 10.^polyval(p2, log10(h_vals));

figure;

loglog(h_vals, errors,'ro','MarkerFaceColor','r','MarkerSize',3, displayname='Forward Euler'); hold on
loglog(h_vals, fit_errors, 'k','LineWidth',1.5, displayname=sprintf('Forward Euler Fit, slope = %.2f', p(1)));

loglog(h_vals, errors2,'go','MarkerFaceColor','g','MarkerSize',3, displayname='Explicit Midpoint'); hold on
loglog(h_vals, fit_errors2, 'blue','LineWidth',1.5, displayname=sprintf('Forward Euler Fit, slope = %.2f', p2(1)));

xlabel('H Values (-)')
ylabel('Epsilon Values (-)')
% title('Log Plot of Rate Func 1')
title('Global Truncation Error Experiment with Fit: Forward Euler and Explicit Midpoint')
legend()