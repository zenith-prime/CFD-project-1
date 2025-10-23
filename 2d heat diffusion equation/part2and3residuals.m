% Load the data
gs_161 = load('gs_161.mat');
tdma_161 = load('tdma_161.mat');

% Extract residuals
rgs = gs_161.residuals_gs; % Assuming the variable is named 'data'
rtd = tdma_161.residuals; % Assuming the variable is named 'data'

% Create a figure for the plot
figure;
hold on; % Hold on to plot both data sets on the same graph

% Plot both residuals
plot(rgs, 'b-', 'DisplayName', 'Gauss-Seidel Residuals'); % 'b-' is for blue line
plot(rtd, 'r-', 'DisplayName', 'TDMA Residuals'); % 'r-' is for red line

% Add titles and labels
title('Residuals vs. Number of Iterations for Different Methods');
xlabel('Iterations');
ylabel('Residuals');
ylim([0,0.5]);

% Add legend
legend('show');

% Enable grid
grid on;

% Hold off to stop adding to the current plot
hold off;
