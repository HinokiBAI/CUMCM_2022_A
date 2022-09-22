% 问题4-2：蒙特卡洛求解最大功率
% 遍历 C
Cs = [];
Ps = [];
Pmax = 0;
Cmax = 0;
Cpmax = 0;
x = 0;
cs = 100000.*rand(2000,1);
cps = 100000.*rand(2000,1);
for i = 1:2000
    fprintf('%d\n', i);
    C = cs(i);
    Cp = cps(i);
    P_sum_v = 0;
    [t2, v_f, v_z] = get_v(C, Cp);
    [Pw] = get_Pw(Cp);
    for t = 70:80
    % 2. 求解数值位移
        [~, index] = min(abs(t-t2));
        v = v_z(index) - v_f(index);
        Pv = C*v^2;
        P_sum_v = P_sum_v + Pv;
    end
    % 3. 求解平均
    P_eval_v = P_sum_v/11;
    P_eval = P_eval_v + Pw;
    Ps = [Ps; P_eval];
    if P_eval > Pmax
        Pmax = P_eval;
        Cmax = C;
        Cpmax = Cp;
    end
end



function [tempP] = get_Pw(Cp)
    w =  1.9806;
    B = 1655.905;
    l1 = 0.5;
    l2 = 1;
    m = 2433*1/3*(l1^2+l2^2+l1*l2);
    A = 7142.493;
    fl = 2140;
    tempP = ((Cp/2)*(fl*fl)) / ((B+Cp)^2 + (w*m+w*A-Cp/w)^2);
end

function [t2, v_f, v_z] = get_v(C, Cp)
    tspan = [70,80];
    y0 = [0 0 0 0 0 0 0 0];
    [t2, x1] = ode45(@(t, x) ode3(t, x, C, Cp), tspan, y0);
    v_f = x1(:, 6);
    v_z = x1(:, 8);
end

% 计算旋转的角度
function dx=ode3(t, x, C, Cp)
    m_z = 2433;
    M = 4688;
    L = 1/12*4866*3^2;
    Ap = 7142.293; % 附加转动惯量
    Bp = 1655.905; % 纵摇兴波阻尼系数
    Kp = 250000;
    A = 1091.099;
    B = 528.5018;
    K = 80000;
    e = 1025;
    g = 9.8;
    S = pi;
    f = 1760;
    fp = 2140;
    w = 1.9806;
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


