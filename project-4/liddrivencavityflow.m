% Lid-Driven Cavity Flow Simulation using SIMPLE Algorithm
clc; clear; close all;

%% Parameters
L = 1;                  % Length of cavity (m)
H = 1;                  % Height of cavity (m)
Re = 100;               % Reynolds number
nu = 1/Re;              % Kinematic viscosity
n = 21;                 % Grid size (41x41)
dx = L / (n-1);         % Grid spacing in x
dy = H / (n-1);         % Grid spacing in y
u_top = 1;              % Lid velocity

% Relaxation parameters
alpha_u = 0.7;          % Under-relaxation for u
alpha_v = 0.7;          % Under-relaxation for v
alpha_p = 0.3;          % Under-relaxation for pressure

max_iter = 1000;        % Maximum iterations
tolerance = 1e-7;       % Convergence criterion

%% Initialization
u = zeros(n+1, n);      % u-velocity on staggered grid
v = zeros(n, n+1);      % v-velocity on staggered grid
p = zeros(n, n);        % Pressure on main grid
u_new = u; v_new = v;

%% Boundary Conditions
u(:, end) = u_top;      % Top wall (lid)
u(:, 1) = 0;            % Bottom wall
u(1, :) = 0;            % Left wall
u(end, :) = 0;          % Right wall
v(:, 1) = 0;            % Bottom wall
v(:, end) = 0;          % Top wall
v(1, :) = 0;            % Left wall
v(end, :) = 0;          % Right wall

%% SIMPLE Algorithm
for iter = 1:max_iter
    % Predictor step: Solve for u* and v*
    [u_star, v_star] = solve_momentum(u, v, p, nu, dx, dy);

    % Pressure correction equation
    [p_corr, u_corr, v_corr] = solve_pressure_correction(u_star, v_star, dx, dy, n);

    % Correct velocities and pressure
    u_new(2:end-1, :) = u_star(2:end-1, :) + u_corr(2:end-1, :);
    v_new(:, 2:end-1) = v_star(:, 2:end-1) + v_corr(:, 2:end-1);
    p = p + alpha_p * p_corr;

    % Check for convergence
    % res = max(abs(p_corr(:)));
    % if res < tolerance
    %     fprintf('Converged after %d iterations\n', iter);
    %     break;
    % end

    % Update variables
    u = alpha_u * u_new + (1 - alpha_u) * u;
    v = alpha_v * v_new + (1 - alpha_v) * v;
end

%% Post-Processing and Visualization
[X, Y] = meshgrid(linspace(0, L, n), linspace(0, H, n));

% Interpolate to cell centers
u_center = 0.5 * (u(1:end-1, :) + u(2:end, :));
v_center = 0.5 * (v(:, 1:end-1) + v(:, 2:end));
velocity_magnitude = sqrt(u_center.^2 + v_center.^2);

% Contour Plots
figure;
contourf(X, Y, u_center', 20, 'LineStyle', 'none');
title('X-Velocity (u)'); colorbar; axis equal tight;

figure;
contourf(X, Y, v_center', 20, 'LineStyle', 'none');
title('Y-Velocity (v)'); colorbar; axis equal tight;

figure;
contourf(X, Y, p', 20, 'LineStyle', 'none');
title('Pressure Field (p)'); colorbar; axis equal tight;

%% Plotting the variation along centerlines

% Create grid points
x = linspace(0, L, n);
y = linspace(0, H, n);

% Interpolated velocities and pressure on the main grid
u_center = 0.5 * (u(1:end-1, :) + u(2:end, :));
v_center = 0.5 * (v(:, 1:end-1) + v(:, 2:end));

% Extract data along centerlines
x_mid_index = round(n/2);
y_mid_index = round(n/2);

u_centerline_x = u_center(x_mid_index, :);
v_centerline_x = v_center(x_mid_index, :);
p_centerline_x = p(x_mid_index, :);

u_centerline_y = u_center(:, y_mid_index);
v_centerline_y = v_center(:, y_mid_index);
p_centerline_y = p(:, y_mid_index);

% Plotting x-velocity along x = 0.5
figure;
plot(y, u_centerline_x, '-o ');
xlabel('y'); ylabel('u'); title('X-Velocity along x = 0.5');
grid on;

% Plotting y-velocity along x = 0.5
figure;
plot(y, v_centerline_x, '-o');
xlabel('y'); ylabel('v'); title('Y-Velocity along x = 0.5');
grid on;

% Plotting pressure along x = 0.5
figure;
plot(y, p_centerline_x, '-o');
xlabel('y'); ylabel('Pressure'); title('Pressure along x = 0.5');
grid on;

% Plotting x-velocity along y = 0.5
figure;
plot(x, u_centerline_y, '-o');
xlabel('x'); ylabel('u'); title('X-Velocity along y = 0.5');
grid on;

% Plotting y-velocity along y = 0.5
figure;
plot(x, v_centerline_y, '-o');
xlabel('x'); ylabel('v'); title('Y-Velocity along y = 0.5');
grid on;

% Plotting pressure along y = 0.5
figure;
plot(x, p_centerline_y, '-o');
xlabel('x'); ylabel('Pressure'); title('Pressure along y = 0.5');
grid on;

%% Plotting Velocity Streamlines
figure;
[X_stream, Y_stream] = meshgrid(linspace(0, L, n), linspace(0, H, n));
u_center = 0.5 * (u(1:end-1, :) + u(2:end, :));
v_center = 0.5 * (v(:, 1:end-1) + v(:, 2:end));

streamslice(X_stream, Y_stream, u_center', v_center');
title('Velocity Streamlines');
xlabel('X');
ylabel('Y');
axis equal tight;
grid on;


%% Function Definitions

% Solve momentum equations
function [u_star, v_star] = solve_momentum(u, v, p, nu, dx, dy)
    [n_u, m_u] = size(u);
    [n_v, m_v] = size(v);
    u_star = u;
    v_star = v;

    % X-Momentum equation
    for i = 2:n_u-1
        for j = 2:m_u-1
            aE = nu / dx^2; aW = nu / dx^2; aN = nu / dy^2; aS = nu / dy^2;
            aP = aE + aW + aN + aS;
            u_star(i, j) = (1 / aP) * ...
                           (aE * u(i+1, j) + aW * u(i-1, j) + ...
                            aN * u(i, j+1) + aS * u(i, j-1) - ...
                            (p(i, j) - p(i-1, j)) / dx);
        end
    end

    % Y-Momentum equation
    for i = 2:n_v-1
        for j = 2:m_v-1
            aE = nu / dx^2; aW = nu / dx^2; aN = nu / dy^2; aS = nu / dy^2;
            aP = aE + aW + aN + aS;
            v_star(i, j) = (1 / aP) * ...
                           (aE * v(i, j+1) + aW * v(i, j-1) + ...
                            aN * v(i+1, j) + aS * v(i-1, j) - ...
                            (p(i, j) - p(i, j-1)) / dy);
        end
    end
end

% Solve pressure correction equation
function [p_corr, u_corr, v_corr] = solve_pressure_correction(u_star, v_star, dx, dy, n)
    p_corr = zeros(n, n);
    u_corr = zeros(size(u_star));
    v_corr = zeros(size(v_star));

    % Define areas
    A_E = dy; A_W = dy; A_N = dx; A_S = dx;

    % Continuity correction
    for i = 2:n-1
        for j = 2:n-1
            div = (u_star(i+1, j) * A_E - u_star(i, j) * A_W) + (v_star(i, j+1) * A_N - v_star(i, j) * A_S);
            p_corr(i, j) = -div / (2 * (1 / dx^2 + 1 / dy^2));
        end
    end

    % Correct velocities
    for i = 2:size(u_star, 1)-1
        for j = 2:size(u_star, 2)-1
            u_corr(i, j) = (p_corr(i, j) - p_corr(i-1, j)) / dx;
        end
    end
    for i = 2:size(v_star, 1)-1
        for j = 2:size(v_star, 2)-1
            v_corr(i, j) = (p_corr(i, j) - p_corr(i, j-1)) / dy;
        end
    end
end