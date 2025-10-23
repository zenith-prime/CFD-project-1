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
    % Parameters
    N = gridsize(p); 
    h = 1 / (N - 1); % Grid spacing
    x = linspace(0, 1, N);
    y = linspace(0, 1, N);
    
    % Initial guess
    phi = zeros(N, N);
    % Apply the boundary conditions
    phi(1, :) = phi_bottom(x);      % Bottom boundary (y = 0) 
    phi(N, :) = phi_top(x);         % Top boundary (y = 1) 
    phi(:, 1) = phi_left(y);        % Left boundary (x = 0)  
    phi(:, N) = phi_right(y);       % Right boundary (x = 1)   
    
    % Precompute source term
    S_values = zeros(N, N);
    for i = 1:N
        for j = 1:N
            S_values(i, j) = S(y(j), x(i));
        end
    end
    
    % Set the tolerance
    tolerance = 1e-6;
    max_iterations = 20000;
    
    % Iterative solution
    main_diag = -4 * ones(N-2, 1);
    sub_diag = ones(N-3, 1);
    super_diag = ones(N-3, 1);
    
    residuals = zeros(max_iterations, 1); % Store residuals
    
    tic;
    for iter = 1:max_iterations
        phi_old = phi; % Store old values of phi for convergence checking
        for i = 2:N-1
            rhs = zeros(N-2, 1);
            for j = 2:N-1
                rhs(j-1) = S_values(i, j) * h^2 - phi(i+1, j) - phi(i-1, j);
                if j == 2
                    rhs(j-1) = rhs(j-1) - phi(i, j-1);
                elseif j == N-1
                    rhs(j-1) = rhs(j-1) - phi(i, j+1);
                end
            end
            phi(i, 2:N-1) = tdma_solver(sub_diag, main_diag, super_diag, rhs);
        end
        
        % Calculate residual
        residual = norm(phi(:) - phi_old(:), inf);
        residuals(iter) = residual;
        
        % Check for convergence
        if residual < tolerance
            fprintf('Converged in %d iterations\n', iter);
            residuals = residuals(1:iter); % Trim residuals array
            break;
        end
        if iter == max_iterations
            warning('Method did not converge in the maximum number of iterations');
        end
    end
    cpu_runtime(p) = toc;
    disp(cpu_runtime)
    [X, Y] = meshgrid(x, y);
    
    % Plot residuals vs. number of iterations
    figure(1);
    plot(1:length(residuals), residuals, '-', 'DisplayName', ['Grid size ' num2str(N)])
    legends{end+1} = ['Grid size ' num2str(N)]; % Add legend entry
    if p == length(gridsize)
        save('tdma_161.mat','residuals')
    end
end

figure(2)
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

figure
contourf(X, Y, phi, 20);
colorbar;
xlabel('x');
ylabel('y');
title('Contour Plot of \phi for the 2D Steady-State Diffusion Equation');

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
