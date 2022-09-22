% 问题2-1：求解线性阻尼下最优功率
clear all;clc;
Ps = [];
Cs = [];
Pmax = 0;
Cmax = 0;
for C=linspace(1, 100000, 1000)
    [v_f, v_z, t] = get_v(C);
    k = 0;
    P_sum = 0;
    for t1 = 70:200
        k = k + 1;
        [~, index] = min(abs(t-t1));
        v_f_1 = v_f(index);
        v_z_1 = v_z(index);
        delta_v = abs(v_f_1-v_z_1);
        P_sum = P_sum + C*delta_v^2;
    end
    P_eval = P_sum/k;
    if P_eval > Pmax
       Pmax = P_eval;
       Cmax = C;
    end
    Ps = [Ps; P_eval];
    Cs = [Cs C];
end
plot(Cs, Ps, 'linewidth',2), xlabel('阻尼系数（N·s/m）'), ylabel('功率（w）')
hold on;scatter(Cmax, Pmax, 200, 'filled', 'p');
hold on;plot([Cmax, Cmax, 0], [0, Pmax, Pmax])

function [v_f, v_z, t]=get_v(C)
    tspan = [70,200];
    y0 = [0 0 0 0];
    [t, x] = ode45(@(t, x) ode1(t, x, C), tspan, y0);
    v_f = x(:, 2);
    v_z = x(:, 4);
end

function dx=ode1(t, x, C)
    m = 4866;
    m_z = 2433;
    A = 1165.992;
    B = 167.8395;
    K = 80000;
    e = 1025;
    g = 9.8;
    S = pi;
    f = 4890;
    w = 2.2143;
    dx = zeros(4, 1);
    dx(1) = x(2);
    dx(2) = (-B*x(2)+C*(x(4)-x(2))+K*(x(3)-x(1))-(e*g*S)*x(1)+f*cos(w*t))/(m+A);
    dx(3) = x(4);
    dx(4) = (-C*(x(4)-x(2))- K*(x(3)-x(1)))/m_z;
end
