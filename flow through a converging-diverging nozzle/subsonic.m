clc;
close all;

% Functions for the vector conserved variables and the flux vector
conserved_vars = @(q1, q2, q3) [q1; q2; q3];
flux_vector = @(f1, f2, f3) [f1; f2; f3];

% Parameters
nx = 100; % Axial gridpoint

% Set check to false to then change max_iter to observe instabilities
max_iter = 1.7e4; % Maximum Iterations
check=true; 

gamma = 1.4; % Adiabatic Ratio
x = linspace(0, 3, nx); 
dx = x(2) - x(1); 
A = ones(1, nx) + 2.2 * (x - 1.5).^2; % Area function
throat=min(find(A==min(A)));

% Initial conditions
rho = ones(1, nx);
T = ones(1,nx);
V = 0.3 ./ (rho .* A);

% Conservative variables
Q = conserved_vars(rho .* A, rho .* A .* V, rho .* A .* (T / (gamma - 1) + gamma * V.^2 / 2));

tic;
for iter = 1:max_iter
    if check && abs(V(throat)/T(throat)-1)<4.5e-2
        break;
    end
    Q_old = Q;

    % Calculate speed of sound and determine time step
    a = sqrt(gamma * T); % Sound speed
    dt = 0.3 * dx / max(V + a); % Courant condition with C = 0.3

    % Predictor step
    F = flux_vector(Q(2, :), ...
        (Q(2, :).^2 ./ Q(1, :)) + (1 - 1 / gamma) * (Q(3, :) - (gamma / 2) * (Q(2, :).^2 ./ Q(1, :))), ...
        (gamma * Q(2, :) .* Q(3, :) ./ Q(1, :)) - (gamma * (gamma - 1) / 2) * (Q(2, :).^3 ./ Q(1, :).^2));

    % Initialize predictors
    dQpdt = zeros(3, nx);

    for i = 2:nx-1
        % Calculate source term for momentum equation
        s2 = (1 / gamma) * rho(i) * T(i) * (A(i+1) - A(i)) / dx;

        % Compute predictor time derivatives
        dQpdt(:, i) = - (F(:, i+1) - F(:, i)) / dx;
        dQpdt(2, i) = dQpdt(2, i) + s2;

        % Update Q values for the predictor step
        Q(:, i) = Q_old(:, i) + dQpdt(:, i) * dt;
    end

    % Update primitive variables after predictor step
    rho = Q(1, :) ./ A;
    V = Q(2, :) ./ Q(1, :);
    T = (Q(3, :) ./ Q(1, :) - (gamma / 2) * V.^2) * (gamma - 1);

    % Update fluxes for corrector step
    F = flux_vector(Q(2, :), ...
        (Q(2, :).^2 ./ Q(1, :)) + (1 - 1 / gamma) * (Q(3, :) - (gamma / 2) * (Q(2, :).^2 ./ Q(1, :))), ...
        (gamma * Q(2, :) .* Q(3, :) ./ Q(1, :)) - (gamma * (gamma - 1) / 2) * (Q(2, :).^3 ./ Q(1, :).^2));

    % Corrector step
    for i = 2:nx-1
        % Calculate source term for the momentum equation
        s2 = (1 / gamma) * rho(i) * T(i) * (A(i) - A(i-1)) / dx;

        % Compute corrector time derivatives
        dQcdt = -(F(:, i) - F(:, i-1)) / dx;
        dQcdt(2) = dQcdt(2) + s2;

        % Average time derivatives and update Q for the corrector step
        Q(:, i) = Q_old(:, i) + 0.5 * (dQpdt(:, i) + dQcdt) * dt;
    end

    % Boundary conditions
    % Inlet
    rho(1) = 1;
    T(1) = 1;
    Q(:, 1) = conserved_vars(rho(1) * A(1), 2 * Q(2, 2) - Q(2, 3), rho(1)...
        * A(1) * (T(1) / (gamma - 1) + gamma * (Q(2, 1) / Q(1, 1))^2 / 2));

    % Outlet
    Q(1, nx) = 2 * Q(1, nx-1) - Q(1, nx-2);
    Q(2, nx) = 2 * Q(2, nx-1) - Q(2, nx-2);
    rho_nx=Q(1,nx)/A(1,nx);
    P_nx=0.99;
    T_nx=P_nx/rho_nx;
    Q(3, nx) = rho_nx * A(1,nx) * (T_nx / (gamma - 1) + gamma * (Q(2, nx) / Q(1, nx))^2 / 2);
    
    % Update primitive variables
    rho = Q(1, :) ./ A;
    V = Q(2, :) ./ Q(1, :);
    T = (Q(3, :) ./ Q(1, :) - (gamma / 2) * V.^2) * (gamma - 1);
    P = rho .* T;
end
cpu_runtime=toc;

[rho_rho0,T_T0,P_P0,M_a]=analytical_subsonic();
% Results
figure
title('Variation of \rho, T, P and M with the Axial Coordinate')
subplot(221)
hold on
plot(1.5*ones(1,100),linspace(0.6,1,100),'r','LineStyle',':','LineWidth',1.2)
plot(x, rho,'r','LineWidth',1.2)
plot(x,rho_rho0,'LineStyle','--','LineWidth',1.2)
ylabel('\rho / \rho_{o}')
xlabel('Axial Coordinate')
title('Variation of Non-Dimensionalised Density')
legend('Throat','Numerical Solution','Analytical Solution');
grid on
hold off

subplot(222)
hold on
plot(1.5*ones(1,100),linspace(0.8,1,100),'g','LineStyle',':','LineWidth',1.2)
plot(x, T,'g','LineWidth',1.2)
plot(x,T_T0,'LineStyle','--','LineWidth',1.2)
ylabel('T / T_{o}')
xlabel('Axial Coordinate')
title('Variation of Non-Dimensionalised Temperature')
legend('Throat','Numerical Solution','Analytical Solution');
grid on
hold off

subplot(223)
hold on
plot(1.5*ones(1,100),linspace(0.5,1,100),'b','LineStyle',':','LineWidth',1.2)
plot(x, P,'b','LineWidth',1.2)
plot(x,P_P0,'LineStyle','--','LineWidth',1.2)
ylabel('P / P_{o}')
xlabel('Axial Coordinate')
title('Variation of Non-Dimensionalised Pressure')
legend('Throat','Numerical Solution','Analytical Solution');
grid on
hold off

subplot(224)
M_n=V ./ sqrt(T);
hold on
plot(1.5*ones(1,100),linspace(0,1,100),'k','LineStyle',':','LineWidth',1.2)
plot(x, M_n,'k','LineWidth',1.2)
plot(x,M_a,'LineStyle','--','LineWidth',1.2)
ylabel('M')
xlabel('Axial Coordinate')
legend('Throat','Numerical Solution','Analytical Solution');
title('Variation of Mach Number')
grid on
hold off
sgtitle('Flow Properties along the Nozzle (Numerical and Analytical Solutions)');

error=[rho_rho0-rho;T_T0-T;P_P0-P;M_a-M_n];
for i=1:4
    rmse(i)=sqrt(sum(error(i,:).^2)/nx);
end
disp('Root Mean Squared Error= ');
rmse
% Analytical Solution Function
function [rho_rho0,T_T0,P_P0,M]=analytical_subsonic()
    % Given parameters
    gamma = 1.4;  % Specific heat ratio for air
    x = linspace(0, 3, 100);  % 100 points on the axis from entrance to exit
    A = 1 + 2.2 * (x - 1.5).^2;  % Area function A/A*
    
    % Initialize arrays for the results
    M = zeros(size(x));
    P_P0 = zeros(size(x));
    T_T0 = zeros(size(x));
    rho_rho0 = zeros(size(x));
    
    % Define an initial guess for Mach number to solve iteratively
    M_guess = 0.1;  % Initial guess for subsonic region
    options = optimset('Display', 'off');  % Suppress fsolve output
    
    % Solve for Mach number at each point along the duct
    for i = 1:length(x)
        % Define the equation to solve based on the area-Mach relation
        area_ratio = A(i);
        eqn = @(M) (1/M) * ((2/(gamma + 1)) * (1 + (gamma - 1)/2 * M^2))^((gamma + 1)/(2 * (gamma - 1))) - area_ratio;
        
        % Use fsolve to find M at each x location
        M(i) = fsolve(eqn, M_guess, options);
        
        % Update the guess for the next point based on continuity (if needed)
        M_guess = M(i);
        
        % Calculate P/P0, T/T0, and rho/rho0 using the analytical expressions
        P_P0(i) = (1 + (gamma - 1)/2 * M(i)^2)^(-gamma / (gamma - 1));
        T_T0(i) = (1 + (gamma - 1)/2 * M(i)^2)^-1;
        rho_rho0(i) = (1 + (gamma - 1)/2 * M(i)^2)^(-1 / (gamma - 1));
    end
end