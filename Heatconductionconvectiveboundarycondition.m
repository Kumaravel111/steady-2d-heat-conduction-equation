clc;            % Clear the command window
clear;          % Clear all variables from the workspace
close all;      % Close all figure windows
% grid poins for L = 1
m = 21; %grid pints along x direction
n = 21; %grid points along y direction
dx = 1.0 / (m - 1); % grid size along x direction
dy = 1.0 / (n - 1); % grid size along y direction
h = 10;  % convective heat transfer coefficient  
k = 1;   % conductive heat transfer coefficient  
tf = 300; % free stream temperature 

%  matrixcs for store the temperature values
T_old = zeros(m, n);
T_new = zeros(m, n);

% Coefficients temperature S,W,P,E,N
S = 1.0 / dy^2;
W = 1.0 / dx^2;
P = 2.0 * ((1.0 / dx^2) + (1.0 / dy^2));
E = 1.0 / dx^2;
N = 1.0 / dy^2;

% Initialize error and iteration count
e = 1; % Start with a large error to enter the loop
itr = 0;% To Know how much iteration was took

% Apply boundary conditions
for i = 1:m
    for j = 1:n
        if j == n
            T_new(i, j) = 500.0; % Right boundary
        elseif i == 1
            T_new(i, j) = 500.0; % Bottom boundary
        elseif i == m
            T_new(i, j) = 500.0; % Lower part boundary
        else
            T_new(i, j) = 0.0; 
        end
    end
end

% solve using jocobi iteration method

while e > 1e-8
  
    T_old = T_new; % replace old temperature with new temperature 

    
    for i = 2:(m - 1) 
        for j = 2:(n - 1) 
            T_new(i, j) = (1.0 / P) * ( ...
                + S * T_old(i, j - 1) ...
                + W * T_old(i - 1, j) ...
                + E * T_old(i + 1, j) ...
                + N * T_old(i, j + 1));
        end
    end

    
    % Mixed or convective boundary condition at the bottom wall (j = 1)
for i = 1:m
    T_new(i, 1) = (T_new(i, 2)+(h*dy/k)*tf) / (1 +  (h * dy / k));
end

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
