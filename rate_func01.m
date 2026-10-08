function dXdt = rate_func0(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end
function X = solution01(t)
X = cos(t);
end

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

close all
X0=1;

fun= @(time, X) rate_func0(time, X);

tspan=[0,10];

h_ref = 0.1;

colors = jet(10); 
i=1;

figure;
for h_ref=linspace(0.01,0.1,10)
    [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(fun,tspan,X0,h_ref);
    plot(t_list,X_list, 'Color', colors(i, :), DisplayName="Numerical, h=" +h_ref); hold on;
    i=i+1;
end

t=linspace(0,10);
X = solution01(t);

plot(t,X, 'r--', LineWidth=2, DisplayName="Analytical")

legend()
title("Forward Euler Integration Plot -5x + 5cos(t) - sin(t);")
xlabel("Time (~)")
ylabel("X Value (~)")