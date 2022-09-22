% ����3�������ҡ�����˶�����
clc;clear all;
% 1. ��⸡�ӵ��˶���ֵ��
tspan = [0,30];
y0 = [0 0 0 0 0 0 0 0];
% opts = odeset('RelTol', 1e-8);
[t1, x_] = ode45(@ode1, tspan, y0);
x1 = x_(:, 1:2);
x2 = x_(:, 3:4);
x3 = x_(:, 5:6);
x4 = x_(:, 7:8);


% 3. ͼ�����
figure(1);
subplot(2, 1, 1);
yyaxis left
plot(t1, x1(:,1)*360/2/pi,'linewidth',2), xlabel('t(s)'), ylabel('\theta(��)');
yyaxis right
plot(t1, x1(:,2)*360/2/pi,'linewidth',2), xlabel('t(s)'), ylabel('\Omega(��/s)'),legend('�Ƕ�\theta','���ٶ�\Omega')
title('������ҡ�˶�λ���ٶ�����')

subplot(2, 1, 2);
yyaxis left
plot(t1, x2(:,1)*360/2/pi,'linewidth',2), xlabel('t(s)'), ylabel('\theta(��)');
yyaxis right
plot(t1, x2(:,2)*360/2/pi,'linewidth',2), xlabel('t(s)'), ylabel('\Omega(��/s)'),legend('�Ƕ�\theta','���ٶ�\Omega')
title('������ҡ�˶�λ���ٶ�����')

figure(2);
subplot(2, 1, 1);
yyaxis left
plot(t1, x3(:,1),'linewidth',2), xlabel('t(s)'), ylabel('x(m)');
yyaxis right
plot(t1, x3(:,2),'linewidth',2), xlabel('t(s)'), ylabel('v(m/s)'),legend('λ��x','�ٶ�v')
title('���Ӵ����˶�λ���ٶ�����')

subplot(2, 1, 2);
yyaxis left
plot(t1, x4(:,1),'linewidth',2), xlabel('t(s)'), ylabel('x(m)');
yyaxis right
plot(t1, x4(:,2),'linewidth',2), xlabel('t(s)'), ylabel('v(m/s)'),legend('λ��x','�ٶ�v')
title('���Ӵ����˶�λ���ٶ�����')

% 4. ������ݼ���
T = 1 / 1.4005;
t_ = 0;
data = [];
for i=1:40
    % ��ȡ���ӵ�λ�ƺ��ٶ�
    [~, index] = min(abs(t1-t_));
    temp_x = x_(index, 5);
    temp_v = x_(index, 6);
    temp_x_2 = x_(index, 1);
    temp_v_2 = x_(index, 2);
    temp_x_3 = x_(index, 7);
    temp_v_3 = x_(index, 8);
    temp_x_4 = x_(index, 3);
    temp_v_4 = x_(index, 4);
    data = [data;t_ temp_x temp_v temp_x_2 temp_v_2 temp_x_3 temp_v_3 temp_x_4 temp_v_4];
    t_ = t_ + 0.2;
end
% 
% 5. �ؼ����ݼ���
key_data = [];
for t_ = [10 20 40 60 100]
    [~, index] = min(abs(t1-t_));
    temp_x = x_(index, 5);
    temp_v = x_(index, 6);
    temp_x_2 = x_(index, 1);
    temp_v_2 = x_(index, 2);
    temp_x_3 = x_(index, 7);
    temp_v_3 = x_(index, 8);
    temp_x_4 = x_(index, 3);
    temp_v_4 = x_(index, 4);
    key_data = [key_data;t_ temp_x temp_v temp_x_2 temp_v_2 temp_x_3 temp_v_3 temp_x_4 temp_v_4];
    t_ = t_ + 0.2;
end

% ������ת�ĽǶ�
function dx=ode1(t, x)
m_z = 2433;
M = 4688;
L = 1/12*4866*3^2;
Ap = 7001.914/12; % ����ת������
Bp = 654.3383; % ��ҡ�˲�����ϵ��
Cp = 1000;  % ��ת������
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

