% �������ͼ����άͼ
clc;clear all;
load('data_50.mat')
figure(1);
surfc(X,Y,Ps)
xlabel('��ת����������ϵ����N��m��s��')
ylabel('ֱ������ϵ����N��s/m��'),
zlabel('���ʣ�w��')

figure(2);
x2 = linspace(0, 100000, 50);
h = heatmap(Ps);
colormap('jet')
h.XLabel = '��ת����ϵ��';
h.YLabel = 'ֱ������ϵ��';

YourYticklabel=cell(size(h.YDisplayLabels))
[YourYticklabel{:}]=deal('');
h.YDisplayLabels=YourYticklabel

YourYticklabel=cell(size(h.XDisplayLabels))
[YourYticklabel{:}]=deal('');
h.XDisplayLabels=YourYticklabel
