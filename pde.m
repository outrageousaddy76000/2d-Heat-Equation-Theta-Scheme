% Ask the user for input values
X = input('Enter the value of X: ');
Y = input('Enter the value of Y: ');
T = input('Enter the value of T: ');
alpha = input('Enter the value of alpha: ');
theta = input('Enter the value of theta: ');
nx = input('Enter the value of nx: ');
ny = input('Enter the value of ny: ');
nt = input('Enter the value of nt: ');
nx=nx+2;
ny=ny+2;
% Calculate grid spacing and time step
dx = X / nx;
dy = Y / ny;
dt = T / nt;

% Calculate stability condition
stability_condition = (1 - 2 * theta) * alpha * dt * ((1 / (dx^2)) + (1 / (dy^2)));

% Check if the scheme is stable
if stability_condition <= 0.5
    disp('The scheme is stable.');
else
    disp('The scheme is unstable. Please provide different configurations.');
    return; % End the program if unstable
end

% Initialize grid
x = linspace(0, X, nx);
y = linspace(0, Y, ny);
t = linspace(0, T, nt);
U = zeros(nx, ny, nt);

% Set initial condition (U(:,:,1)) and boundary conditions
[x, y] = meshgrid(linspace(0, X, nx), linspace(0, Y, ny));
%U(:,:,1) = exp(-((x - X/2).^2 / (X/5)^2 + (y - Y/2).^2 / (Y/5)^2));
U(:,:,1) = sin(pi * x / X) .* sin(pi * y / Y);
%U(1:floor(nx/2), :) = 1;
%U(:,:,1) = rand(nx, ny);
%U(:,:,1) = square(2 * pi * x / X) .* square(2 * pi * y / Y);

% Neumann boundary conditions
U(1,:,:) = U(3,:,:); % Left boundary
U(nx,:,:) = U(nx-2,:,:); % Right boundary
U(:,1,:) = U(:,3,:); % Bottom boundary
U(:,ny,:) = U(:,ny-2,:); % Top boundary

% Create the coefficient matrix A
A = sparse(nx * ny, nx * ny);
lambda = alpha * dt * ((1 / (dx^2)) + (1 / (dy^2)));
diagonals = [ones(nx*ny, 1), -lambda*theta*ones(nx*ny, 1), -lambda*theta*ones(nx*ny, 1), ...
             -lambda*theta*ones(nx*ny, 1), -lambda*theta*ones(nx*ny, 1)];
offsets = [0, -1, 1, -ny, ny];

% Populate the matrix A
A = spdiags(diagonals, offsets, nx * ny, nx * ny);
A = A + spdiags(ones(nx * ny, 1) * (1 + 4 * lambda * theta), 0, nx * ny, nx * ny);
%disp(full(A));
% Time-stepping loop to solve the equation AU^{(n+1)} = F(U^{(n)})
for k = 2:nt
    % Initialize F as a zeros matrix with the same size as U(:,:,k)
    F = zeros(nx, ny);
    
    for i = 2:nx-1
        for j = 2:ny-1
            % Calculate F for each cell using the provided formula
            F(i, j) = U(i, j, k-1) + (1 - theta) * lambda * (U(i+1, j, k-1) - 4 * U(i, j, k-1) + U(i-1, j, k-1) + U(i, j+1, k-1) + U(i, j-1, k-1));
        end
    end
    
    % Flatten the 2D array for matrix multiplication
    U_temp = reshape(U(:,:,k-1), [], 1);
    
    % Solve the linear system A * U_temp = F
    U_temp = A \ F(:);
    
    % Reshape the 1D array back to 2D
    U(:,:,k) = reshape(U_temp, nx, ny);
end


% Plot temperature distribution at T/2
figure;
subplot(1, 2, 1);
time_index_T_over_2 = round(nt / 2); % Index for T/2
surf(x, y, U(:,:,time_index_T_over_2));
xlabel('X');
ylabel('Y');
zlabel('Temperature');
title('Temperature Distribution at T/2');

% Plot temperature distribution at T
subplot(1, 2, 2);
surf(x, y, U(:,:,nt));
xlabel('X');
ylabel('Y');
zlabel('Temperature');
title('Temperature Distribution at T');

% Adjust the subplot layout
sgtitle('Temperature Distributions at Different Time Steps');
% Calculate grid spacing and time step for refined grid
dx_refined = X / (20 * nx);
dy_refined = Y / (20 * ny);

% Initialize the refined grid
x_refined = linspace(0, X, 20 * nx);
y_refined = linspace(0, Y, 20 * ny);
U_refined = zeros(20 * nx, 20 * ny, nt);

% Set initial condition for refined grid
[x_refined, y_refined] = meshgrid(linspace(0, X, 20 * nx), linspace(0, Y, 20 * ny));
U_refined(:,:,1) = exp(-((x_refined - X/2).^2 / (X/5)^2 + (y_refined - Y/2).^2 / (Y/5)^2));

% Dirichlet boundary conditions for refined grid
U_refined(1,:,:) = 0; % Left boundary
U_refined(20 * nx,:,:) = 0; % Right boundary
U_refined(:,1,:) = 0; % Bottom boundary
U_refined(:,20 * ny,:) = 0; % Top boundary

% Create the coefficient matrix A for refined grid
A_refined = sparse(20 * nx * 20 * ny, 20 * nx * 20 * ny);
lambda_refined = alpha * dt * ((1 / (dx_refined^2)) + (1 / (dy_refined^2)));

% Populate the matrix A for refined grid
diagonals_refined = [ones(20 * nx * 20 * ny, 1), -lambda_refined * theta * ones(20 * nx * 20 * ny, 1), -lambda_refined * theta * ones(20 * nx * 20 * ny, 1), ...
             -lambda_refined * theta * ones(20 * nx * 20 * ny, 1), -lambda_refined * theta * ones(20 * nx * 20 * ny, 1)];
offsets_refined = [0, -1, 1, -20 * ny, 20 * ny];

A_refined = spdiags(diagonals_refined, offsets_refined, 20 * nx * 20 * ny, 20 * nx * 20 * ny);
A_refined = A_refined + spdiags(ones(20 * nx * 20 * ny, 1) * (1 + 4 * lambda_refined * theta), 0, 20 * nx * 20 * ny, 20 * nx * 20 * ny);

% Time-stepping loop for refined grid
for k = 2:nt
    % Initialize F for refined grid
    F_refined = zeros(20 * nx, 20 * ny);
    
    for i = 2:20 * nx - 1
        for j = 2:20 * ny - 1
            % Calculate F for each cell using the provided formula
            F_refined(i, j) = U_refined(i, j, k-1) + (1 - theta) * lambda_refined * (U_refined(i+1, j, k-1) - 4 * U_refined(i, j, k-1) + U_refined(i-1, j, k-1) + U_refined(i, j+1, k-1) + U_refined(i, j-1, k-1));
        end
    end
    
    % Flatten the 2D array for matrix multiplication
    U_temp_refined = reshape(U_refined(:,:,k-1), [], 1);
    
    % Solve the linear system A_refined * U_temp_refined = F_refined
    U_temp_refined = A_refined \ F_refined(:);
    
    % Reshape the 1D array back to 2D for refined grid
    U_refined(:,:,k) = reshape(U_temp_refined, 20 * nx, 20 * ny);
end

% Calculate the error matrix
error_matrix = zeros(nx, ny, nt);
for k = 1:nt
    for i = 1:nx
        for j = 1:ny
            % Calculate the error at each spatial point and time step
            error_matrix(i, j, k) = abs(U_refined((20*(i-1))+1, (20*(j-1))+1, k) - U(i, j, k));
        end
    end
end
% Calculate the error at time T
error_at_T = error_matrix(:, :, nt);
error_at_T(:,nx)=0;
error_at_T(ny,:)=0;
error_at_T(1,:)=0;
error_at_T(:,1)=0;
% Plot the error_matrix at time T
figure;
surf(x, y, error_at_T);
xlabel('X');
ylabel('Y');
zlabel('Error');
title('Error at Time T');

% Initialize a matrix to store the differences for the specified row and column
differences = zeros(nt, 1);

% Specify the row and column of interest
row_of_interest = floor(nx/2);
column_of_interest = floor(ny/2);

for k = 3:nt
    % Calculate the logarithmic rate of convergence for the specified row and column
    error_time = error_matrix(row_of_interest, column_of_interest, k);
    error_time_minus_1 = error_matrix(row_of_interest, column_of_interest, k-1);
    error_time_minus_2 = error_matrix(row_of_interest, column_of_interest, k-2);
    
    % Check if the denominators are zero (to avoid division by zero)
    if error_time_minus_1 == 0 || error_time_minus_2 == 0
        differences(k) = 0;
    else
        % Calculate the logarithmic rate of convergence
        differences(k) = log(error_time / error_time_minus_1) / log(error_time_minus_1 / error_time_minus_2);
    end
end

% Specify the time and y-coordinate of interest
time_of_interest = floor(nt / 2);
y_coordinate_of_interest = floor(ny / 2);
x_coordinate_of_interest = floor(nx / 2);
% Initialize a matrix to store the differences for the specified row and column
dx_differences = zeros(nx, 1);

for i = 3:nx
    % Calculate the error at the specified time and y-coordinate for the current and previous grids
    error_current = error_matrix(i, y_coordinate_of_interest, time_of_interest);
    error_previous = error_matrix(i-1, y_coordinate_of_interest, time_of_interest);
    error_previous_2 = error_matrix(i-2, y_coordinate_of_interest, time_of_interest);
    % Check if the denominators are zero (to avoid division by zero)
    if error_previous == 0 || error_previous_2==0
        dx_differences(i) = 0;
    else
        % Calculate the order of convergence in dx
        dx_differences(i) = abs(log(error_current / error_previous) / log(error_previous / error_previous_2));
    end
end

% Initialize a matrix to store the differences for the specified row and column
dy_differences = zeros(ny, 1);

for j = 3:ny
    % Calculate the error at the specified time and x-coordinate for the current and previous grids
    error_current = error_matrix(x_coordinate_of_interest, j, time_of_interest);
    error_previous = error_matrix(x_coordinate_of_interest, j-1, time_of_interest);
    error_previous_2 = error_matrix(x_coordinate_of_interest, j-2, time_of_interest);
    % Check if the denominators are zero (to avoid division by zero)
    if error_previous == 0 || error_previous_2 == 0
        dy_differences(j) = 0;
    else
        % Calculate the order of convergence in dy
        dy_differences(j) = abs(log(error_current / error_previous) / log(error_previous / error_previous_2));
    end
end

% Calculate and plot the logarithmic rate of convergence for the specified row and column
figure;
plot(3:nt, differences(3:end));
xlabel('Time Step');
ylabel('Logarithmic Rate of Convergence');
title('Logarithmic Rate of Convergence for Specified Row and Column');
% Calculate and plot the mean value of differences
mean_differences = mean(differences(3:end));
hold on;
plot([3, nt], [mean_differences, mean_differences], '--r', 'LineWidth', 2);
legend('Logarithmic Rate of Convergence', 'Mean Value');
hold off;

% Calculate and plot the order of convergence in dx
figure;
plot(3:nx, dx_differences(3:end));
xlabel('Number of Grid Points (nx)');
ylabel('Order of Convergence in dx');
title('Order of Convergence in dx for Specified Time and Y-coordinate');
% Calculate and plot the mean value of dx_differences
mean_dx_differences = mean(dx_differences(3:end));
hold on;
plot([3, nx], [mean_dx_differences, mean_dx_differences], '--r', 'LineWidth', 2);
legend('Order of Convergence in dx', 'Mean Value');
hold off;

% Calculate and plot the order of convergence in dy
figure;
plot(3:ny, dy_differences(3:end));
xlabel('Number of Grid Points (ny)');
ylabel('Order of Convergence in dy');
title('Order of Convergence in dy for Specified Time and X-coordinate');
% Calculate and plot the mean value of dy_differences
mean_dy_differences = mean(dy_differences(3:end));
hold on;
plot([3, ny], [mean_dy_differences, mean_dy_differences], '--r', 'LineWidth', 2);
legend('Order of Convergence in dy', 'Mean Value');
hold off;