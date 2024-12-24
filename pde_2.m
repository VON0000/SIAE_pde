% 参数设置
a = 1; % 热扩散系数
T = 5; % 总时间
Nr = 20; % 径向网格数
Nt = 5000; % 时间步数
dr = 1 / (Nr - 1); % 径向步长
dt = T / Nt; % 时间步长
lambda = a^2 * dt / dr^2;

% 稳定性检查
if lambda > 0.5
    error('稳定性条件不满足，请减小 dt 或增加 dr！');
end

% 径向网格
rho = linspace(0, 1, Nr); % 径向距离 [0, 1]

% 初值
u = sin(rho.^2); % 初始条件 u(r, 0)
u_next = zeros(size(u)); % 下一时间步的解

% 时间步循环
for n = 1:Nt
    t = n * dt; % 当前时间

    % 内部网格更新
    for i = 2:Nr-1
        u_next(i) = u(i) + lambda * ( ...
            (rho(i+1)^2 * u(i+1) - 2 * rho(i)^2 * u(i) + rho(i-1)^2 * u(i-1)) / (rho(i)^2 * dr^2) ...
        );
    end

    % 边界条件
    u_next(1) = 0; % 在原点处的边界条件
    u_next(end) = t * sin(1); % 球面处的边界条件

    % 更新解
    u = u_next;

    % 可视化（每隔若干步绘制一次）
    if mod(n, 500) == 0
        % 绘制球体内部的径向分布
        [sphereX, sphereY, sphereZ] = sphere(50);
        surf(sphereX, sphereY, sphereZ, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % 透明球面
        hold on;
        
        % 绘制径向分布
        plot3(rho .* sin(linspace(0, pi, Nr)), ...
              rho .* cos(linspace(0, pi, Nr)), ...
              u, 'r', 'LineWidth', 2); % 径向分布
        hold off;
        title(['Time: ', num2str(t)]);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        colorbar;
        drawnow;
    end
end
