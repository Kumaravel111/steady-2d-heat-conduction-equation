clc;            % Clear the command window
clear;          % Clear all variables from the workspace
close all;      % Close all figure windows
% grid poins for L = 1m
m = 21; %grid points along x direction
n = 21; %grid points along y direction
dx = 1.0 / (m - 1); %grid size along x direction
dy = 1.0 / (n - 1); %grid size along y direction

% matricxs for store temperature values
T_old = zeros(m, n);
T_new = zeros(m, n);

% Coefficients of temperature
S = 1.0 / dy^2;
W = 1.0 / dx^2;
P = 2.0 * ((1.0 / dx^2) + (1.0 / dy^2));
E = 1.0 / dx^2;
N = 1.0 / dy^2;

% Initialize error and iteration count
e = 1; % Start with a large error to enter the loop
itr = 0;% to know how much itration was took

% Apply boundary conditions
for i = 1:m
    for j = 1:n
        if i == 1
            T_new(i, j) = 500.0; % Right boundary
        elseif i == m
            T_new(i, j) = 1000.0; % Bottom boundary
        else
           T_new(i, j) = 0.0; 
        end
    end
end

% solve using jacobi iteration method
while e > 1e-8
  
    T_old = T_new; % replace with old temperature with new temperature

    
    for i = 2:(m - 1) 
        for j = 2:(n - 1) 
            T_new(i, j) = (1.0 / P) * ( ...
                + S * T_old(i, j - 1) ...
                + W * T_old(i - 1, j) ...
                + E * T_old(i + 1, j) ...
                + N * T_old(i, j + 1));
        end
    end

    % Apply Neumann boundary condition at the top and bottom

  T_new(:,n) = T_new(:,n-1);
  T_new(:,1) = T_new(:,2);


    % Calculate the error
    e = 0; 
    for i = 1:m
        for j = 1:n
            e = e + (T_new(i, j) - T_old(i, j))^2;
        end
    end
    e = sqrt(e / (m * n));

    % Increment iteration count
    itr = itr + 1;
end

% Display results

disp(['Converged in ', num2str(itr), ' iterations.']);
disp(['Final Error: ', num2str(e)]);

% Generate grid for plotting
x = linspace(0, 1, m);
y = linspace(0, 1, n);
[X, Y] = meshgrid(x, y);

% Plot temperature contours as filled color regions
figure;
contourf(X, Y, T_new', 50, 'LineStyle', 'none'); % Transpose T_new for correct orientation
colorbar;
caxis([min(T_new(:)), max(T_new(:))]); % Set color scale to data range
title('Temperature Distribution');
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
grid on;
