clc;
gridsize = [41, 81, 161];
cpu_runtime = zeros(3, 1);

% Fill the source term S(x, y)
S = @(x, y) 50000 * exp(-50 * ((1 - x)^2 + y^2)) * (100 * ((1 - x)^2 + y^2) - 2);

% Boundary conditions
phi_left = @(y) 500 * exp(-50 * (1 + y.^2)); % Left boundary (x = 0)
phi_right = @(y) 100 * (1 - y) + 500 * exp(-50 * y.^2); % Right boundary (x = 1)
phi_bottom = @(x) 100 * x + 500 * exp(-50 * (1 - x).^2); % Bottom boundary (y = 0)
phi_top = @(x) 500 * exp(-50 * ((1 - x).^2 + 1)); % Top boundary (y = 1)


figure;
hold on; % Hold on to plot multiple lines on the same plot
legends = {}; % To store legend entries

for p = 1:length(gridsize)
    N = gridsize(p);
    % Parameters
    h = 1 / (N - 1); % Grid spacing
    x = linspace(0, 1, N);
    y = linspace(0, 1, N);

    % Initialize phi
    [X, Y] = meshgrid(x, y);
    phi = zeros(N, N);
    phi(1, :) = phi_left(y);         % Bottom boundary (y = 0)
    phi(N, :) = phi_right(y);       % Top boundary (y = 1)
    phi(:, 1) = phi_bottom(x);       % Left boundary (x = 0)
    phi(:, N) = phi_top(x);         % Right boundary (x = 1)

    % Solve using Gauss-Seidel method
    tic;
    [phi, iter, residuals_gs] = gauss_seidel_solver(phi, S, h, N);
    cpu_runtime(p) = toc;

    % Plot residuals vs. number of iterations
    plot(1:length(residuals_gs), residuals_gs, '-', 'DisplayName', ['Grid size ' num2str(N)])
    legends{end+1} = ['Grid size ' num2str(N)]; % Add legend entry

    % Plot the solution as a contour plot
    if p == length(gridsize)
        figure
        contourf(X, Y, phi.', 20);
        colorbar;
        xlabel('x');
        ylabel('y');
        title('Contour Plot of \phi for the 2D Steady-State Diffusion Equation');
        save('gs_161.mat','residuals_gs')
    end
end

figure
plot(gridsize, cpu_runtime, '-o')
title('CPU Runtime vs. Gridsize')
xlabel('Gridsize')
ylabel('CPU Runtime (in s)')
grid on
disp(cpu_runtime)

% Finalize the residuals plot
figure(1);
title('Residuals vs. Number of Iterations for Different Grid Sizes')
xlabel('Iterations')
ylabel('Residuals')
ylim([0, 1])
legend(legends)
grid on
hold off

function [phi, iter, residuals_gs] = gauss_seidel_solver(phi, S, h, N, tol, max_iter)
    if nargin < 5
        tol = 1e-6;
    end
    if nargin < 6
        max_iter = 40000;
    end

    residuals_gs = zeros(max_iter, 1);
    for iter = 1:max_iter
        phi_old = phi; % Store the old solution
        
        % Gauss-Seidel iteration for interior points
        for i = 2:N-1
            for j = 2:N-1
                phi(i, j) = 0.25 * (phi(i+1, j) + phi(i-1, j) + phi(i, j+1) + phi(i, j-1) - h^2 * S((i-1)*h, (j-1)*h));
            end
        end

        % Compute the residual (infinity norm)
        residuals_gs(iter) = norm(phi(:) - phi_old(:), inf);
        
        % Check for convergence
        if residuals_gs(iter) < tol
            residuals_gs = residuals_gs(1:iter);
            return;
        end
    end
    residuals_gs = residuals_gs(1:iter);
    warning('Gauss-Seidel iterative method did not converge within the maximum number of iterations')
end
