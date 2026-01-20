points = importdata('simulated_points.csv');
x_points = points(:, 1);
y_points = points(:, 2);
z_points = points(:, 3);
% 绘制点集
figure;
scatter3(x_points, y_points, z_points, 'b', 'filled');
hold on;

% 绘制斜截面
slope_angle = -pi/6; 
slope_height = 30;

slope_equation = @(x, y) tan(slope_angle) * (x - 0.5) + slope_height;

[X, Y] = meshgrid(linspace(min(x_points), max(x_points), 100), linspace(min(y_points), max(y_points), 100));
Z = slope_equation(X, Y);
surf(X, Y, Z, 'FaceAlpha', 0.5);

% 取出斜截面上的点
d = 5; % 取出点的距离
slope_points = [X(:), Y(:), slope_equation(X(:), Y(:))];
A = tan(slope_angle);
B = 0;
C = -1;
D = -0.5*tan(slope_angle)+slope_height;

distances = abs(A * x_points + B * y_points + C * z_points + D) / sqrt(A^2 + B^2 + C^2);

% 找到距离小于给定距离的点
selected_points = points(distances < d, :);
writematrix(selected_points, 'm_slope.csv');
writematrix(selected_points, '../data/m_slope.csv');

% 绘制取出的点
scatter3(selected_points(:,1), selected_points(:,2), selected_points(:,3), 'r', 'filled');

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Points, Sloped Plane, and Selected Points');
legend('Point Set', 'Sloped Plane', 'Selected Points');
grid on;
