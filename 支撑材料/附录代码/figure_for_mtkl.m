clc;clear all;
load('data_mtkl.mat')

figure(1);
scatter(cps, cs, 8, 'red', 'filled'),
xlabel('��ת����ϵ��'),ylabel('ֱ������ϵ��'), title('���ؿ������Ͷ��ͼ');

figure(2);
s = scatter3(cps, cs, Ps, 8, 'red', 'filled'), xlabel('��ת����ϵ����N��m��s��'),
ylabel('ֱ������ϵ����N��s/m��'), zlabel('���ʣ�w��');
