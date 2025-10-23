clc;
% Parameters
gridsize = [21,41,81]; 
cpu_runtime=zeros(3,1); 
% Fill the source term S(x,y)
S = @(x, y) 50000 * exp(-50 * ((1 - x)^2 + y^2)) * (100 * ((1 - x)^2 + y^2) - 2);

% Boundary conditions
phi_left = @(y) 500 * exp(-50 * (1 + y.^2)); % Left boundary (x = 0)
phi_right = @(y) 100 * (1 - y) + 500 * exp(-50 * y.^2); % Right boundary (x = 1)
phi_bottom = @(x) 100 * x + 500 * exp(-50 * (1 - x).^2); % Bottom boundary (y = 0)
phi_top = @(x) 500 * exp(-50 * ((1 - x).^2 + 1)); % Top boundary (y = 1)

for p=1:1:length(gridsize)
    N=gridsize(p);

    % Initialize A and b
    A = zeros(N*N, N*N);
    b = zeros(N*N, 1);
    h = 1 / (N - 1); % Grid spacing
    x = linspace(0, 1, N);
    y = linspace(0, 1, N);
    [X,Y] = meshgrid(x,y);
    % Mapping function to convert 2D index (i, j) to 1D index
    index = @(i, j) (i - 1) * N + j;

    % Fill matrices A and b
    for i = 1:N
        for j = 1:N
            idx = index(i, j);
            xi = x(i);
            yj = y(j);
            
            if i == 1 % Left boundary (x = 0)
                A(idx, idx) = 1;
                b(idx) = phi_left(yj);
            elseif i == N % Right boundary (x = 1)
                A(idx, idx) = 1;
                b(idx) = phi_right(yj);
            elseif j == 1 % Bottom boundary (y = 0)
                A(idx, idx) = 1;
                b(idx) = phi_bottom(xi);
            elseif j == N % Top boundary (y = 1)
                A(idx, idx) = 1;
                b(idx) = phi_top(xi);
            else % Interior points
                A(idx, idx) = -4; % Center
                A(idx, index(i-1, j)) = 1; % Left
                A(idx, index(i+1, j)) = 1; % Right
                A(idx, index(i, j-1)) = 1; % Bottom
                A(idx, index(i, j+1)) = 1; % Top
                b(idx) = h^2 * S(xi, yj); % Source term
            end
        end
    end
    tic;
    % The system A*phi_flattened = b is now ready to be solved using Gaussian elimination
    phi_flattened = gaussElimination(A,b);
    cpu_runtime(p)=toc;
    % Reshape the solution into a 2D grid
    phi_2D = reshape(phi_flattened, [N, N]);
    % Plot the solution as a contour plot
    if p == length(gridsize)
        figure
        contourf(X, Y, phi_2D, 20);
        colorbar;
        xlabel('x');
        ylabel('y');
        title('Contour Plot of \phi for the 2D Steady-State Diffusion Equation');
    end
    disp(cpu_runtime)
end

figure
plot(gridsize,cpu_runtime, '-o')
title('CPU Runtime vs. Gridsize')
xlabel('Gridsize')
ylabel('CPU Runtime (in s)')
grid on
disp(cpu_runtime)

function x = gaussElimination(A, b)
    % Augment matrix A with vector b
    [n, ~] = size(A);
    AugmentedMatrix = [A, b];
    
    % Forward elimination
    for k = 1:n
        % Pivot
        [~, i_max] = max(abs(AugmentedMatrix(k:n, k)));
        i_max = i_max + k - 1;
        
        % Swap rows
        if k ~= i_max
            AugmentedMatrix([k, i_max], :) = AugmentedMatrix([i_max, k], :);
        end
        
        % Eliminate
        for i = k+1:n
            factor = AugmentedMatrix(i, k) / AugmentedMatrix(k, k);
            AugmentedMatrix(i, k:end) = AugmentedMatrix(i, k:end) - factor * AugmentedMatrix(k, k:end);
        end
    end
    
    % Back substitution
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (AugmentedMatrix(i, end) - AugmentedMatrix(i, i+1:end-1) * x(i+1:end)) / AugmentedMatrix(i, i);
    end
end


