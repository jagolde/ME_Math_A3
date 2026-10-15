function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end
function X = solution01(t)
X = cos(t);
end
close all

G_func = @(t, X0, h_ref) forward_euler_fixed_step_integration(new_rate,t,X0,h_ref);

new_rate= @(t, X) rate_func01(t, X);
sol3= @(t) solution01(t);

h_vals = logspace(-5, 1, 100);
t_ref = 0.492;
e_list= [];
figure;
for h = h_vals
    epsilon = local_truncation_error(G_func, sol3, t_ref, h);
    e_list(end+1) = epsilon;
end

p = polyfit(h_vals, e_list, 1);



loglog(h_vals, e_list); hold on
loglog(h_vals, exp(polyval(p, h_vals)), 'b--', 'LineWidth', 2);


xlabel('H Values')
ylabel('Epsilon Values')
title('Log Plot of Rate Func 1')