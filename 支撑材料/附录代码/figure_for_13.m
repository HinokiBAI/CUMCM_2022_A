% 求解垂荡与纵摇直线运动对比图
clc;clear all;
% 1. 求解浮子的运动数值解
tspan = [0,30];
y0 = [0 0 0 0 0 0 0 0];
% opts = odeset('RelTol', 1e-8);
[t1, x_] = ode45(@ode1, tspan, y0);
x3 = x_(:, 5:6);
x4 = x_(:, 7:8);

y0 = [0 0 0 0];
[t2, x_] = ode45(@ode2, tspan, y0);
x1 = x_(:, 1:2);


% 3. 图像绘制
figure(1);
subplot(2, 1, 1);
yyaxis left
plot(t2, x1(:,1),'linewidth',2), xlabel('t(s)'), ylabel('x(m)');
yyaxis right
plot(t2, x1(:,2),'linewidth',2), xlabel('t(s)'), ylabel('v(m/s)'),legend('位移x','速度v')
title('浮子垂荡运动位移速度曲线')

subplot(2, 1, 2);
yyaxis left
plot(t1, x3(:,1),'linewidth',2), xlabel('t(s)'), ylabel('x(m)');
yyaxis right
plot(t1, x3(:,2),'linewidth',2), xlabel('t(s)'), ylabel('v(m/s)'),legend('位移x','速度v')
title('浮子垂荡运动修正位移速度曲线')


% 计算旋转的角度
function dx=ode1(t, x)
m_z = 2433;
M = 4688;
L = 1/12*4866*3^2;
Ap = 7001.914/12; % 附加转动惯量
Bp = 654.3383; % 纵摇兴波阻尼系数
Cp = 1000;  % 旋转阻尼器
Kp = 250000;
A = 1028.876;
B = 654.3383;
C = 10000;
K = 80000;
e = 1025;
g = 9.8;
S = pi;
f = 3640;
fp = 1690;
w = 1.7152;
Cc = 8890.7;
dx = zeros(8, 1);
l1 = 0.5 + (x(3)-x(1)) + 0.5;
l2 = 0.5 + (x(3)-x(1));
L_z = 1/3*(l1-l2)*(l1^2+l1*l2+l2^2);
dx(1) = x(2);
dx(2) = (-Bp*x(2)+(Cp*(x(4)-x(2))+Kp*(x(3)-x(1)))*x(1)-(Cc)*x(1)+fp*cos(w*t))/(L+Ap);
dx(3) = x(4);
dx(4) = (-Cp*(x(4)-x(2))- Kp*(x(3)-x(1)))/L_z;
dx(5) = x(6);
dx(6) = (-B*x(6)+(C*(x(8)-x(6))+K*(x(7)-x(5)))*x(1)-(e*g*S)*x(5)+f*cos(w*t))/(M+A);
dx(7) = x(8);
dx(8) = (-C*(x(8)-x(6))- K*(x(7)-x(5)))/m_z;
end

function dx=ode2(t, x)
m = 4866;
m_z = 2433;
A = 1028.876;
B = 654.3383;
C = 10000;
K = 80000;
e = 1025;
g = 9.8;
S = pi;
f = 3640;
w = 1.7152;
dx = zeros(4, 1);
dx(1) = x(2);
dx(2) = (-B*x(2)+C*(x(4)-x(2))+K*(x(3)-x(1))-(e*g*S)*x(1)+f*cos(w*t))/(m+A);
dx(3) = x(4);
dx(4) = (-C*(x(4)-x(2))- K*(x(3)-x(1)))/m_z;
end

