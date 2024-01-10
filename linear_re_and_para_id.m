clc;
clear;

%% ��������
load 3V1Hz���Ҳ�����ģ������(1000).mat
load 3V1Hz���Ҳ�����ģ�����(1000).mat
u=in.signals.values;
y=out1.signals.values(:,1);

%% �ݶ��½���
% ������ʼ��
l = length(u);
yo = zeros(l, 1);
keseio = zeros(l, 1);
thetacell = cell(l, 1);
alpha = 40^2;
F = alpha * eye(6);
thetaini = zeros(6, 1);
thetacell{1} = thetaini;

% ��ʼ����
for i = 1:l-3
    phi = [y(i+2); y(i+1); y(i); u(i+2); u(i+1); u(i)];
    yo(i) = thetacell{i}' * phi;
    keseio(i) = y(i+3) - yo(i);
    thetacell{i+1} = thetacell{i} + F * phi * keseio(i) / (1 + phi' * F * phi);
end
theta_p1 = thetacell{i};

% ʹ�ù��Ʋ�������Ԥ��
yo_estimated1 = zeros(size(y));
for i = 1:length(y)-3
    phi_estimated = [y(i+2); y(i+1); y(i); u(i+2); u(i+1); u(i)];
    yo_estimated1(i) = theta_p1' * phi_estimated;
end

% ��������
y_true = y(4:end);
y_look = yo_estimated1(1:end-3);
MAE = mean(abs(y_true - y_look));
RMSE = sqrt(mean((y_true - y_look).^2));

disp(['�ݶ��½���:']);
disp(['MAE: ', num2str(MAE)]);
disp(['RMSE: ', num2str(RMSE)]);

%% ��С���˷�
N = 998;
phiN = zeros(N, 6);
yN = zeros(N, 1);

% ������С���˷�����
for i = 1:N
    phi = [y(i+2); y(i+1); y(i); u(i+2); u(i+1); u(i)]';
    phiN(i, :) = phi;
    yN(i) = y(i+3);
end

% ʹ����С���˷�������
thetals = inv(phiN' * phiN) * phiN' * yN;
theta_p2 = thetals;

% ʹ�ù��Ʋ�������Ԥ��
yo_estimated2 = zeros(size(y));
for i = 1:length(y)-3
    phi_estimated = [y(i+2); y(i+1); y(i); u(i+2); u(i+1); u(i)];
    yo_estimated2(i) = theta_p2' * phi_estimated;
end

% ��������
y_look = yo_estimated2(1:end-3);
MAE = mean(abs(y_true - y_look));
RMSE = sqrt(mean((y_true - y_look).^2));

disp(['��С���˷�:']);
disp(['MAE: ', num2str(MAE)]);
disp(['RMSE: ', num2str(RMSE)]);

%% ������С���˷�
% ��ʼ��
thetap = zeros(l, 6);
kesaio = zeros(l, 1);
alpha = 40^2;
F = alpha * eye(6);

% ����
for i = 3:l-3
    phi = [y(i); y(i-1); y(i-2); u(i); u(i-1); u(i-2)];
    kesaio(i+1) = y(i+1) - thetap(i, :) * phi;
    F = F - F * phi * phi' * F / (1 + phi' * F * phi);
    thetap(i+1, :) = thetap(i, :) + (F * phi * kesaio(i+1))';
end
theta_p3 = thetap(i, :);

% ʹ�ù��Ʋ�������Ԥ��
yo_estimated3 = zeros(size(y));
for i = 1:length(y)-3
    phi_estimated = [y(i+2); y(i+1); y(i); u(i+2); u(i+1); u(i)];
    yo_estimated3(i) = theta_p3 * phi_estimated;
end

% ��������
y_look = yo_estimated3(1:end-3);
MAE = mean(abs(y_true - y_look));
RMSE = sqrt(mean((y_true - y_look).^2));

disp(['������С���˷�:']);
disp(['MAE: ', num2str(MAE)]);
disp(['RMSE: ', num2str(RMSE)]);

t=in.time;
ylook=yo_estimated2;
% ������ʵֵ�͹���ֵ
figure;
plot(t, y, 'b-', 'LineWidth', 2);  % ������ʵֵ��ʹ����ɫʵ��
hold on;  % ���ֵ�ǰͼ�Σ�ʹ����һ��ͼ�ο��Ի�����ͬһ����ϵ��
plot(t, ylook, 'r--', 'LineWidth', 2);  % ���ƹ���ֵ��ʹ�ú�ɫ����

% ��ӱ�ǩ��ͼ��
xlabel('u');  % x���ǩ
ylabel('y');  % y���ǩ
title('��ʵֵ�͹���ֵ�Ļ���');  % ͼ����
legend('��ʵֵ', '����ֵ');  % ͼ��

% ��ʾ����
grid on;

% ������ͼ
hold off;
