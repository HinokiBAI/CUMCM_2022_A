% 问题1-1：求解垂荡浮子运动方程
clc;clear all;
% 1. 求解浮子的运动数值解
tspan = [0,30];
y0 = [1 0 -1 0];
opts = odeset('RelTol', 1e-8);
[t, x] = ode45(@ode1, tspan, y0, opts);
x1 = x(:, 1:2);
x2 = x(:, 3:4);
x_f = x1(:, 1);
v_f = x1(:, 2);
x_z = x2(:, 1);
v_z = x2(:, 2);

% 3. 图像绘制
figure(1);
subplot(2, 1, 1);
yyaxis left
plot(t, x1(:,1),'linewidth',2), xlabel('t(s)'), ylabel('x(m)');
yyaxis right
plot(t, x1(:,2),'linewidth',2), xlabel('t(s)'), ylabel('v(m/s)'),legend('位移x','速度v')
title('浮子运动位移速度曲线')

subplot(2, 1, 2);
yyaxis left
plot(t, x2(:,1),'linewidth',2), xlabel('t(s)'), ylabel('x(m)');
yyaxis right
plot(t, x2(:,2),'linewidth',2), xlabel('t(s)'), ylabel('v(m/s)'),legend('位移x','速度v')

title('振子运动位移速度曲线')

% 4. 表格数据计算
T = 1 / 1.4005;
t_ = 0;
data = [];
for i=1:40
    % 获取浮子的位移和速度
    [~, index] = min(abs(t-t_));
    temp_x = x_f(index);
    temp_v = v_f(index);
    % 获取振子的位移和速度
    temp_x_2 = x_z(index);
    temp_v_2 = v_z(index);
    data = [data;t_ temp_x temp_v temp_x_2 temp_v_2];
    t_ = t_ + 0.2;
end
% 
% 5. 关键数据计算
key_data = [];
for t_ = [10 20 40 60 100]
    % 获取浮子的位移和速度
    [~, index] = min(abs(t-t_));
    temp_x = x_f(index);
    temp_v = v_f(index);
    % 获取振子的位移和速度
    temp_x_2 = x_z(index);
    temp_v_2 = v_z(index);
    key_data = [key_data;t_ temp_x temp_v temp_x_2 temp_v_2];
end

function dx=ode1(t, x)
m = 4866;
m_z = 2433;
A = 1335.535;
B = 656.3616;
C = 10000;
K = 80000;
e = 1025;
g = 9.8;
S = pi;
f = 6250;
w = 1.4005;
dx = zeros(4, 1);
dx(1) = x(2);
dx(2) = (-B*x(2)+C*(x(4)-x(2))+K*(x(3)-x(1))-(e*g*S)*x(1)+f*cos(w*t))/(m+A);
dx(3) = x(4);
dx(4) = (-C*(x(4)-x(2))- K*(x(3)-x(1)))/m_z;
end

