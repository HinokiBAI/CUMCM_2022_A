% 问题2-2：求解非线性阻尼下最优功率
clear all;clc;
Pmax = 0;
Cmax = 0;
for alpha = [0.4]
    Ps = [];
    Cs = [];
    for C=linspace(1, 100000, 1000)
        [v_f, v_z, t] = get_v(C, alpha);
        P_sum = 0;
        for t1 = 70:200
            [~, index] = min(abs(t-t1));
            v_z_1 = v_z(index);
            v_f_1 = v_f(index);
            delta_v = abs(v_z_1-v_f_1);
            P_sum = P_sum + C*(delta_v^alpha)*delta_v^2;
        end
        P_eval = P_sum/131;
        Ps = [Ps; P_eval];
        Cs = [Cs C];
       if P_eval > Pmax
       Pmax = P_eval;
       Cmax = C;
        end
    end
    if alpha ~= 0
        hold on;
    end
    plot(Cs, Ps, 'linewidth',2), xlabel('阻尼系数（N·s/m）'), ylabel('功率（w）')
    disp('完成一个')
end


function [v_f, v_z, t]=get_v(C, alpha)
    tspan = [70,200];
    y0 = [0 0 0 0];
    [t, x] = ode45(@(t, x) ode1(t, x, C, alpha), tspan, y0);
    v_f = x(:, 2);
    v_z = x(:, 4);
end

function dx=ode1(t, x, C, alpha)
    m = 4866;
    m_z = 2433;
    A = 1165.992;
    B = 167.8395;
    e = 1025;
    g = 9.8;
    S = pi;
    f = 4890;
    w = 2.2143;
    K = 80000;
    dx = zeros(4, 1);
    dx(1) = x(2);
    dx(2) = (-B*x(2)+C*(abs(x(2)-x(4))^alpha)*(x(4)-x(2))+K*(x(3)-x(1))-(e*g*S)*x(1)+f*cos(w*t))/(m+A);
    dx(3) = x(4);
    dx(4) = (-C*(abs(x(2)-x(4))^alpha)*(x(4)-x(2))- K*(x(3)-x(1)))/m_z;
end
