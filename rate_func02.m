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
G_func = @(t, X0, h_ref) forward_euler_step(new_rate,t,X0,h_ref);

sol3= @(t) solution01(t);

h_vals = logspace(-5, 1, 100);
t_ref = 0.492;
e_list= [];
for h = h_vals
    epsilon = local_truncation_error(G_func, sol3, t_ref, h);
    e_list(end+1) = epsilon;
end

errors = abs(e_list);
valid = (errors > 0) & isfinite(errors);

p = polyfit(log10(h_vals(valid)), log10(errors(valid)), 1);
fit_errors = 10.^polyval(p, log10(h_vals));



G_func = @(t, X0, h_ref) explicit_midpoint_step(new_rate,t,X0,h_ref);

sol3= @(t) solution01(t);

h_vals = logspace(-5, 1, 100);
t_ref = 0.492;
e_list= [];
for h = h_vals
    epsilon = local_truncation_error(G_func, sol3, t_ref, h);
    e_list(end+1) = epsilon;
end

errors = abs(e_list);
valid = (errors > 0) & isfinite(errors);

p = polyfit(log10(h_vals(valid)), log10(errors(valid)), 1);
fit_errors = 10.^polyval(p, log10(h_vals));









figure;
loglog(h_vals, errors); hold on
loglog(h_vals, fit_errors, 'r-');

xlabel('H Values')
ylabel('Epsilon Values')
title('Log Plot of Rate Func 1')
legend('Computed error', sprintf('Fit: slope = %.2f', p(1)), ...
       'Location', 'best')



xlabel('Step Size (~)')
ylabel('Measured Local Truncation Error $\epsilon$ (~)')
title('Log Plot of Rate Func 1')