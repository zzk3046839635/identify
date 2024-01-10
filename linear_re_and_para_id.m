clc;
clear;

%% 加载数据
load 3V1Hz正弦波物理模型输入(1000).mat
load 3V1Hz正弦波物理模型输出(1000).mat
u=in.signals.values;
y=out1.signals.values(:,1);

%% 梯度下降法
% 参数初始化
l = length(u);
yo = zeros(l, 1);
keseio = zeros(l, 1);
thetacell = cell(l, 1);
alpha = 40^2;
F = alpha * eye(6);
thetaini = zeros(6, 1);
thetacell{1} = thetaini;

% 开始迭代
for i = 1:l-3
    phi = [y(i+2); y(i+1); y(i); u(i+2); u(i+1); u(i)];
    yo(i) = thetacell{i}' * phi;
    keseio(i) = y(i+3) - yo(i);
    thetacell{i+1} = thetacell{i} + F * phi * keseio(i) / (1 + phi' * F * phi);
end
theta_p1 = thetacell{i};

% 使用估计参数进行预测
yo_estimated1 = zeros(size(y));
for i = 1:length(y)-3
    phi_estimated = [y(i+2); y(i+1); y(i); u(i+2); u(i+1); u(i)];
    yo_estimated1(i) = theta_p1' * phi_estimated;
end

% 评估性能
y_true = y(4:end);
y_look = yo_estimated1(1:end-3);
MAE = mean(abs(y_true - y_look));
RMSE = sqrt(mean((y_true - y_look).^2));

disp(['梯度下降法:']);
disp(['MAE: ', num2str(MAE)]);
disp(['RMSE: ', num2str(RMSE)]);

%% 最小二乘法
N = 998;
phiN = zeros(N, 6);
yN = zeros(N, 1);

% 构建最小二乘法矩阵
for i = 1:N
    phi = [y(i+2); y(i+1); y(i); u(i+2); u(i+1); u(i)]';
    phiN(i, :) = phi;
    yN(i) = y(i+3);
end

% 使用最小二乘法求解参数
thetals = inv(phiN' * phiN) * phiN' * yN;
theta_p2 = thetals;

% 使用估计参数进行预测
yo_estimated2 = zeros(size(y));
for i = 1:length(y)-3
    phi_estimated = [y(i+2); y(i+1); y(i); u(i+2); u(i+1); u(i)];
    yo_estimated2(i) = theta_p2' * phi_estimated;
end

% 评估性能
y_look = yo_estimated2(1:end-3);
MAE = mean(abs(y_true - y_look));
RMSE = sqrt(mean((y_true - y_look).^2));

disp(['最小二乘法:']);
disp(['MAE: ', num2str(MAE)]);
disp(['RMSE: ', num2str(RMSE)]);

%% 迭代最小二乘法
% 初始化
thetap = zeros(l, 6);
kesaio = zeros(l, 1);
alpha = 40^2;
F = alpha * eye(6);

% 迭代
for i = 3:l-3
    phi = [y(i); y(i-1); y(i-2); u(i); u(i-1); u(i-2)];
    kesaio(i+1) = y(i+1) - thetap(i, :) * phi;
    F = F - F * phi * phi' * F / (1 + phi' * F * phi);
    thetap(i+1, :) = thetap(i, :) + (F * phi * kesaio(i+1))';
end
theta_p3 = thetap(i, :);

% 使用估计参数进行预测
yo_estimated3 = zeros(size(y));
for i = 1:length(y)-3
    phi_estimated = [y(i+2); y(i+1); y(i); u(i+2); u(i+1); u(i)];
    yo_estimated3(i) = theta_p3 * phi_estimated;
end

% 评估性能
y_look = yo_estimated3(1:end-3);
MAE = mean(abs(y_true - y_look));
RMSE = sqrt(mean((y_true - y_look).^2));

disp(['迭代最小二乘法:']);
disp(['MAE: ', num2str(MAE)]);
disp(['RMSE: ', num2str(RMSE)]);

t=in.time;
ylook=yo_estimated2;
% 绘制真实值和估计值
figure;
plot(t, y, 'b-', 'LineWidth', 2);  % 绘制真实值，使用蓝色实线
hold on;  % 保持当前图形，使得下一个图形可以绘制在同一坐标系下
plot(t, ylook, 'r--', 'LineWidth', 2);  % 绘制估计值，使用红色虚线

% 添加标签和图例
xlabel('u');  % x轴标签
ylabel('y');  % y轴标签
title('真实值和估计值的绘制');  % 图标题
legend('真实值', '估计值');  % 图例

% 显示网格
grid on;

% 结束绘图
hold off;
