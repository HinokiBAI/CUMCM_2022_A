% 求解热力图和三维图
clc;clear all;
load('data_50.mat')
figure(1);
surfc(X,Y,Ps)
xlabel('旋转阻尼器阻尼系数（N·m·s）')
ylabel('直线阻尼系数（N·s/m）'),
zlabel('功率（w）')

figure(2);
x2 = linspace(0, 100000, 50);
h = heatmap(Ps);
colormap('jet')
h.XLabel = '旋转阻尼系数';
h.YLabel = '直线阻尼系数';

YourYticklabel=cell(size(h.YDisplayLabels))
[YourYticklabel{:}]=deal('');
h.YDisplayLabels=YourYticklabel

YourYticklabel=cell(size(h.XDisplayLabels))
[YourYticklabel{:}]=deal('');
h.XDisplayLabels=YourYticklabel
