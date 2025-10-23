clc;
close all

% Grid Parameters
nx = 10000; ny = 200;
L = 1; H = 0.2;
dx = L / (nx - 1); dy = H / (ny - 1);

% Parameters
Re = 1e4; U_inf = 1;
nu = U_inf * L / Re;

tolerance = 1e-6;
max_iter = 10000;

% Initialize u and v velocity fields
u = U_inf * ones(ny, nx); v = zeros(ny, nx);

% Boundary Conditions
u(:, 1) = U_inf; u(1,:) = 0; u(end, :) = U_inf;

tic;
for j = 2:nx
    u(:,j) = u(:,j-1);
    p = 0;
    while (true)
        u_old = u(:,j);
        for i = 2:ny-1
            u(i,j) = ((u(i,j-1)^2 / dx) - u(i+1,j) * (v(i,j-1) / (2 * dy) ...
              - nu / dy^2) + u(i-1,j) * (v(i,j-1) / (2 * dy) ...
              + nu / dy^2)) / (u(i,j-1) / dx + 2 * nu / dy^2);
        end
        residual = sqrt(sum((u(:,j) - u_old).^2));
        p = p + 1;
        if(residual < tolerance || p >= max_iter)
            break
        end
    end

    v(:,j) = v(:,j-1);
    p = 0;
    while (true)
        v_old = v(:,j);
        for i = 2:ny-1
            v(i,j) = v(i-1,j) + (u(i,j) - u(i,j-1)) * dy / dx;
        end
        residual = sqrt(sum((v(:,j) - v_old).^2));
        p = p + 1;
        if(residual < tolerance || p >= max_iter)
            break;
        end
    end
end

% Arrays for plotting
l1 = [[],0];
l2 = [[],1];

for j = 2:nx
    for i = 2:ny-1
        if u(i,j) >= 0.9999 * U_inf
            l1 = [l1,i];
            l2 = [l2,j];
            break;
        end
    end
end
cpu_runtime = toc;

% Plotting u-velocity contour
figure;
contourf(linspace(0, L, nx), linspace(0, H, ny), u, 20, 'EdgeColor','none', 'EdgeColor','none');
colorbar;
hold on;
plot((L / nx) * l2, (H / ny) * l1, 'Color','#000000');
title('Contour of u-velocity');
xlabel('x');
ylabel('y');
grid on;

% Plotting v-velocity contour
figure;
contourf(linspace(0, L, nx), linspace(0, H, ny), v, 'EdgeColor','none');
colorbar;
title('Contour of v-velocity');
xlabel('x');
ylabel('y');
grid on;

% Compute and plot normalized x-velocity (F') as a function of similarity variable (eta)
eta = linspace(0, H * sqrt(U_inf / nu / L), ny);
F_prime = u(:, end) / U_inf;

figure;
plot(eta, F_prime, '-','Color', "g");
hold on;
eta_max = 10; % Analytical Solution for F'
eta = linspace(0, eta_max, nx);
blasius_ode = @(eta, F) [F(2); F(3); -0.5 * F(1) * F(3)];
F0 = [0, 0, 0.332];
[eta, F] = ode45(blasius_ode, eta, F0);
plot(eta, F(:, 2));
title('Normalized x-velocity F''(\eta) vs. Similarity Variable \eta');
xlabel('Similarity Variable \eta');
ylabel('Normalized x-velocity F''(\eta)');
legend('Simulated F''', 'Theoretical F''')
grid on;

% Compute and plot boundary layer thickness (delta)
delta_x = 4.91 * sqrt(nu * linspace(0, L, nx) / U_inf);
figure;
plot(linspace(0, L, nx), delta_x, 'Color', "g");
hold on;
plot((L/nx) * l2, (H/ny) * l1, 'Color',"r");
title('Boundary Layer Thickness \delta vs. x');
xlabel('x');
ylabel('Boundary Layer Thickness \delta');
legend('Analytical \delta', 'Simulated \delta')
grid on;
hold off;

disp(cpu_runtime)