% 参数设置
a = 1; % 热扩散系数
T = 5; % 总时间
L = 1; % 空间范围 [-L, L]
Nh = 21; % 空间网格数
Nt = 5000; % 时间步数
dh = 2 * L / (Nh - 1); % 空间步长
dt = T / Nt; % 时间步长
lambda = a^2 * dt / dh^2;

% 稳定性检查
if lambda > 1/6
    error('稳定性条件不满足，请减小 dt 或增加 dh！');
end

% 网格
x = linspace(-L, L, Nh);
y = linspace(-L, L, Nh);
z = linspace(-L, L, Nh);
[X, Y, Z] = meshgrid(x, y, z);

% 初值
u = sin(X.^2 + Y.^2 + Z.^2); % 初值 u(x, y, z, 0)
u_next = zeros(size(u)); % 下一时间步的解

% 可视化参数
num_plots = 4; % 子图数量，表示展示几个时间点的结果
time_intervals = round(linspace(1, Nt, num_plots)); % 在时间步中均匀选取展示时间点
figure; % 创建新图窗

% 时间步循环
for n = 1:Nt
    t = n * dt; % 当前时间
    
    % 内部网格的更新（只对球体内部点进行更新）
    for i = 2:Nh-1
        for j = 2:Nh-1
            for k = 2:Nh-1
                % 判断当前点是否在球体内部
                if X(i, j, k)^2 + Y(i, j, k)^2 + Z(i, j, k)^2 <= 1
                    % 检查 x, y, z 方向上的相邻点是否超出边界
                    if X(i+1, j, k)^2 + Y(i+1, j, k)^2 + Z(i+1, j, k)^2 > 1
                        u_ip1 = t * sin(X(i+1, j, k)^2 + Y(i+1, j, k)^2);
                    else
                        u_ip1 = u(i+1, j, k);
                    end
                    
                    if X(i-1, j, k)^2 + Y(i-1, j, k)^2 + Z(i-1, j, k)^2 > 1
                        u_im1 = t * sin(X(i-1, j, k)^2 + Y(i-1, j, k)^2);
                    else
                        u_im1 = u(i-1, j, k);
                    end
                    
                    if X(i, j+1, k)^2 + Y(i, j+1, k)^2 + Z(i, j+1, k)^2 > 1
                        u_jp1 = t * sin(X(i, j+1, k)^2 + Y(i, j+1, k)^2);
                    else
                        u_jp1 = u(i, j+1, k);
                    end
                    
                    if X(i, j-1, k)^2 + Y(i, j-1, k)^2 + Z(i, j-1, k)^2 > 1
                        u_jm1 = t * sin(X(i, j-1, k)^2 + Y(i, j-1, k)^2);
                    else
                        u_jm1 = u(i, j-1, k);
                    end
                    
                    if X(i, j, k+1)^2 + Y(i, j, k+1)^2 + Z(i, j, k+1)^2 > 1
                        u_kp1 = t * sin(X(i, j, k+1)^2 + Y(i, j, k+1)^2);
                    else
                        u_kp1 = u(i, j, k+1);
                    end
                    
                    if X(i, j, k-1)^2 + Y(i, j, k-1)^2 + Z(i, j, k-1)^2 > 1
                        u_km1 = t * sin(X(i, j, k-1)^2 + Y(i, j, k-1)^2);
                    else
                        u_km1 = u(i, j, k-1);
                    end

                    % 更新公式
                    u_next(i, j, k) = u(i, j, k) + lambda * ( ...
                        u_ip1 + u_im1 + ...
                        u_jp1 + u_jm1 + ...
                        u_kp1 + u_km1 - ...
                        6 * u(i, j, k));
                end
            end
        end
    end
    
    % 处理边界条件（直接赋值球面点）
    for i = 1:Nh
        for j = 1:Nh
            for k = 1:Nh
                if abs(X(i, j, k)^2 + Y(i, j, k)^2 + Z(i, j, k)^2 - 1) < 1e-6
                    u_next(i, j, k) = t * sin(X(i, j, k)^2 + Y(i, j, k)^2);
                end
            end
        end
    end
    
    % 更新解
    u = u_next;

    % 可视化（每隔若干步绘制一次）
    if ismember(n, time_intervals)
        subplot(2, 2, find(time_intervals == n)); % 创建 2x2 的子图布局
        % 绘制切片
        slice(X, Y, Z, u, 0, 0, 0); % 绘制三个平面上的切片
        title(['Time: ', num2str(t)]);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        colorbar;
        hold on;

        % 绘制透明球 x^2 + y^2 + z^2 = 1
        [sphereX, sphereY, sphereZ] = sphere(50); % 生成球面坐标，细分为 50x50
        sphereX = sphereX * L; % 缩放到半径为 L 的球
        sphereY = sphereY * L;
        sphereZ = sphereZ * L;
         
        % 绘制球
        h = surf(sphereX, sphereY, sphereZ, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % 设置透明度为 0.3
        set(h, 'FaceColor', 'flat');

        hold off;
    end
end
