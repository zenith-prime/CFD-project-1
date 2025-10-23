clc;

% Parameters
Nx = 41;
Ny = 81;
x = linspace(0, 1, Nx);
y = linspace(0, 1, Ny);
dx = 1 / (Nx - 1);
dy = 1 / (Ny - 1);

% Fill the source term S(x,y)
S = @(x, y) 50000 * exp(-50 * ((1 - x).^2 + y.^2)) .* (100 * ((1 - x).^2 + y.^2) - 2);

% Boundary conditions
phi_left = @(y) 500 * exp(-50 * (1 + y.^2)); % Left boundary (x = 0)
phi_right = @(y) 100 * (1 - y) + 500 * exp(-50 * y.^2); % Right boundary (x = 1)
phi_bottom = @(x) 100 * x + 500 * exp(-50 * (1 - x).^2); % Bottom boundary (y = 0)
phi_top = @(x) 500 * exp(-50 * ((1 - x).^2 + 1)); % Top boundary (y = 1)

% Initial guess
phi = zeros(Nx, Ny);

% Apply the boundary conditions
phi(1, :) = phi_left(y);         % Bottom boundary (y = 0)
phi(Nx, :) = phi_right(y);       % Top boundary (y = 1)
phi(:, 1) = phi_bottom(x);       % Left boundary (x = 0)
phi(:, Ny) = phi_top(x);         % Right boundary (x = 1)
phi2 = phi;
phi3 = phi;

% Precompute source term
S_values = zeros(Nx, Ny);
for i = 1:Nx
    for j = 1:Ny
        S_values(i, j) = S(x(i), y(j));
    end
end

% Set the tolerance
tolerance = 1e-6;
max_iterations = 10000;

% Iterative solution using ADI method
main_diag_x = -2 * (1 + (dx/dy)^2) * ones(Nx-2, 1);
sub_diag_x = ones(Nx-3, 1);
super_diag_x = ones(Nx-3, 1);

main_diag_y = -2 * (1 + (dy/dx)^2) * ones(Ny-2, 1);
sub_diag_y = ones(Ny-3, 1);
super_diag_y = ones(Ny-3, 1);

residuals = zeros(2 * max_iterations, 1); % Store residuals
residuals2 = zeros(max_iterations, 1); % Store residuals
residuals3 = zeros(max_iterations, 1); % Store residuals

tic;
for iter = 1:max_iterations
    % Store old values of phi for convergence checking
    phi_old = phi;

    % Column sweep
    phi = column_sweep(phi, S_values, x, y, dx, dy, Nx, Ny, main_diag_x, sub_diag_x, super_diag_x);

    % Calculate residual for column sweep
    residual_column = norm(phi - phi_old, inf);
    residuals(2*iter - 1) = residual_column;

    phi_old = phi;

    % Row sweep
    phi = row_sweep(phi, S_values, x, y, dx, dy, Nx, Ny, main_diag_y, sub_diag_y, super_diag_y);

    % Calculate residual for row sweep
    residual_row = norm(phi - phi_old, inf);
    residuals(2*iter) = residual_row;

    % Check for convergence
    if residual_row < tolerance
        fprintf('ADI method converged in %d iterations\n', 2*iter);
        residuals = residuals(1:2*iter); % Trim residuals array
        break;
    end
    if iter == max_iterations
        warning('ADI method did not converge in the maximum number of iterations');
    end
end
cpu_runtime = toc;
disp(cpu_runtime);

[X, Y] = meshgrid(x, y);
figure
contourf(X, Y, phi', 20); % Transpose phi for correct contour plot
colorbar;
xlabel('x');
ylabel('y');
title('Contour Plot of \phi for the 2D Steady-State Diffusion Equation');

% For row sweep
for iter = 1:max_iterations
    phi_old = phi2; % Store old values of phi for convergence checking
    phi2 = row_sweep(phi2, S_values, x, y, dx, dy, Nx, Ny, main_diag_y, sub_diag_y, super_diag_y);
    residual = norm(phi2 - phi_old, inf);
    residuals2(iter) = residual;
    if residual < tolerance
        fprintf('Row sweep method converged in %d iterations\n', iter);
        residuals2 = residuals2(1:iter); % Trim residuals array
        break;
    end
    if iter == max_iterations
        warning('Row sweep method did not converge in the maximum number of iterations');
    end
end

% For column sweep
for iter = 1:max_iterations
    phi_old = phi3; % Store old values of phi for convergence checking
    phi3 = column_sweep(phi3, S_values, x, y, dx, dy, Nx, Ny, main_diag_x, sub_diag_x, super_diag_x);
    residual = norm(phi3 - phi_old, inf);
    residuals3(iter) = residual;
    if residual < tolerance
        fprintf('Column sweep method converged in %d iterations\n', iter);
        residuals3 = residuals3(1:iter); % Trim residuals array
        break;
    end
    if iter == max_iterations
        warning('Column sweep method did not converge in the maximum number of iterations');
    end
end

% Plot residuals vs. number of iterations
figure;
hold on; % Hold on to plot both data sets on the same graph
plot(1:2:length(residuals), residuals(1:2:end), 'b-', 'DisplayName', 'ADI Residuals - Column Sweep');
plot(2:2:length(residuals), residuals(2:2:end), 'r-', 'DisplayName', 'ADI Residuals - Row Sweep');
plot(1:length(residuals2), residuals2, 'g-', 'DisplayName', 'Row-Sweep Residuals');
plot(1:length(residuals3), residuals3, 'k-', 'DisplayName', 'Column-Sweep Residuals');
title('Residuals vs. Number of Iterations for Different Methods');
xlabel('Iterations');
ylabel('Residuals');
ylim([0, 1])
legend('show');
grid on;
hold off;

% Functions for row and column sweeps
function phi = column_sweep(phi, S_values, x, y, dx, dy, Nx, Ny, main_diag, sub_diag, super_diag)
    for j = 2:Ny-1
        rhs = zeros(Nx-2, 1);
        for i = 2:Nx-1
            rhs(i-1) = S_values(i, j) * dx^2 - (phi(i, j+1) + phi(i, j-1)) * (dx/dy)^2;
            if i == 2
                rhs(i-1) = rhs(i-1) - phi(i-1, j);
            elseif i == Nx-1
                rhs(i-1) = rhs(i-1) - phi(i+1, j);
            end
        end
        phi(2:Nx-1, j) = tdma_solver(sub_diag, main_diag, super_diag, rhs);
    end
end

function phi = row_sweep(phi, S_values, x, y, dx, dy, Nx, Ny, main_diag, sub_diag, super_diag)
    for i = 2:Nx-1
        rhs = zeros(Ny-2, 1);
        for j = 2:Ny-1
            rhs(j-1) = S_values(i, j) * dy^2 - (phi(i+1, j) + phi(i-1, j)) * (dy/dx)^2;
            if j == 2
                rhs(j-1) = rhs(j-1) - phi(i, j-1);
            elseif j == Ny-1
                rhs(j-1) = rhs(j-1) - phi(i, j+1);
            end
        end
        phi(i, 2:Ny-1) = tdma_solver(sub_diag, main_diag, super_diag, rhs);
    end
end

% TDMA solver function
function x = tdma_solver(a, b, c, d)
    n = length(b);  % number of equations
    
    % Forward elimination
    for i = 2:n
        w = a(i-1) / b(i-1);
        b(i) = b(i) - w * c(i-1);
        d(i) = d(i) - w * d(i-1);
    end
    
    % Backward substitution
    x = zeros(n, 1);
    x(n) = d(n) / b(n);
    for i = n-1:-1:1
        x(i) = (d(i) - c(i) * x(i+1)) / b(i);
    end
end
