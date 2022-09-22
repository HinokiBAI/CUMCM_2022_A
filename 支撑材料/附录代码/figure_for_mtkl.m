clc;clear all;
load('data_mtkl.mat')

figure(1);
scatter(cps, cs, 8, 'red', 'filled'),
xlabel('旋转阻尼系数'),ylabel('直线阻尼系数'), title('蒙特卡洛随机投点图');

figure(2);
s = scatter3(cps, cs, Ps, 8, 'red', 'filled'), xlabel('旋转阻尼系数（N·m·s）'),
ylabel('直线阻尼系数（N·s/m）'), zlabel('功率（w）');
