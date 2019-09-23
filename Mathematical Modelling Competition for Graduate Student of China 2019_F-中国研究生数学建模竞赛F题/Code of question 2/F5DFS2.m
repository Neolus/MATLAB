%% 2019年第十六届中国研究生数学建模竞赛F题
%  题目：多约束条件下智能飞行器航迹快速规划
%  问题二代码
%  Copyright by @Neolus, @Jia, @Team 19102480059 and @Mathwork Co.
%  Date：2019/09
%  

clear all
close all
clc
%% original data read , Pretreatment and Visualization

dataset_1 = xlsread('附件1：数据集1-终稿.xlsx');
% dataset_1 = xlsread('附件2：数据集2-终稿.xlsx');
horizontal_check = dataset_1(dataset_1(:,5) == 0,:);
vertical_check = dataset_1(dataset_1(:,5) == 1,:);
end_point = dataset_1([1,end],:);

% figure (1),
% scatter3(horizontal_check(:,2)',horizontal_check(:,3)',horizontal_check(:,4)','b','k')
% hold on
% scatter3(vertical_check(:,2)',vertical_check(:,3)',vertical_check(:,4)','g','k')
% hold on
% scatter3(end_point(:,2)',end_point(:,3)',end_point(:,4)','ro','MarkerFaceColor','r')
% title('original data scatter','FontName','Times New Roman','FontWeight','Bold','FontSize',16)
% xlabel('x','FontName','Times New Roman','FontSize',14)
% ylabel('y','FontName','Times New Roman','FontSize',14,'Rotation',0)
% ylabel('y','FontName','Times New Roman','FontSize',14,'Rotation',0)
% hold on 
% plot3(end_point(:,2),end_point(:,3),end_point(:,4),'-r.','LineWidth',2)
%legend('a','Location','best')
% hold off

%% relative distance matrix creation
distmat = squareform( pdist(dataset_1 (:,[2,3,4])) );
m=size(distmat,1)
R=200; %转弯半径
%%

global disMat2;
disMat2=distmat;
for k=1:m
    disMat2(k,k)=inf;
end


global Position;
Position=dataset_1(:,2:4);

%数据集1，计算速度很快
RV_VError=25; RV_HError=15; %垂直校正参数 
RH_VErrror=20; RH_HError=25;%水平校正参数
theta=30; delta=0.001;

% % 数据集二，计算非常慢
% RV_VError=20; RV_HError=10;
% RH_VErrror=15; RH_HError=20;
% theta=20; delta=0.001;

%最大飞行距离，邻接矩阵
% MaxDisFly=theta/delta;
% AdjMat=disMat2<MaxDisFly;
% m=size(AdjMat,1);
% for k=1:m
%     AdjMat(k,k)=0;
% end
% global VisFlag;
global ResFlag;
VisFlag=zeros(1,m); %访问标记
ResFlag=dataset_1 (:,5); %1表示垂直误差校正点，0表示水平误差校正点
global V;
global TR;
global fobj;
% 
V=14; TR=10;%设飞行速度为14m/s,一次校正10s
fobj=10*distmat(1,m); %目标函数初始值，10倍距离
% fobj=10*distmat(1,m)/V; %距离加校正次数模型


%%
RootAB=zeros(1,m);
PointsC=zeros(m,3);
global RootABG; %保存最优路径
RootABG=zeros(1,m);
global PointsCG; %保存所有切点
PointsCG=zeros(m,3);

VelocityIni=[inf,inf,inf]; %飞往下一点的速度方向
VerA= struct('No',{1},'PreNo',{1},'Dis',{0},'NR',{0},'Obj',{0},'HError',{0},'VError',{0},'Pos',{Position(1,:)},'Velo',{VelocityIni},'DisPre2Now',inf); %初始化A点

dfsDrive(VerA,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,R,PointsC); %开始深度搜索
%%
DisBIdx=find(RootABG==m);
SA=RootABG(1:DisBIdx); %提取路径
PointsRev=dataset_1(SA,2:4); %校正点坐标
PointsRad=PointsCG; %脱离点坐标
% PointsRad=PointsRad(1:size(PointsRev,1),:);
PointsRad(1,:)=PointsRev(1,:); 
PointsRad(end,:)=PointsRev(end,:);

RecNo=size(PointsRad,1);

PointsAll=zeros(2*RecNo,3);
for k=1:RecNo %根据运动路径重新排序
    PointsAll(2*k-1,:)=PointsRev(k,:);
    PointsAll(2*k,:)=PointsRad(k,:);
end
%%

figure (2),
scatter3(horizontal_check(:,2)',horizontal_check(:,3)',horizontal_check(:,4)','b','k')
hold on
scatter3(vertical_check(:,2)',vertical_check(:,3)',vertical_check(:,4)','g','k')
hold on
scatter3(end_point(:,2)',end_point(:,3)',end_point(:,4)','ro','MarkerFaceColor','c')
title('original data scatter','FontName','Times New Roman','FontWeight','Bold','FontSize',16)
xlabel('x','FontName','Times New Roman','FontSize',14)
ylabel('y','FontName','Times New Roman','FontSize',14,'Rotation',0)
ylabel('y','FontName','Times New Roman','FontSize',14,'Rotation',0)
hold on 
plot3(end_point(:,2),end_point(:,3),end_point(:,4),'-r.','LineWidth',2)

for k=1:RecNo-1
    plot3(PointsAll(2*k-1:2*k,1),PointsAll(2*k-1:2*k,2),PointsAll(2*k-1:2*k,3),'-r.','LineWidth',2) %绘制路径，颜色符号，请参照相关要求
    plot3(PointsAll(2*k:2*k+1,1),PointsAll(2*k:2*k+1,2),PointsAll(2*k:2*k+1,3),'-b.','LineWidth',2) %绘制路径，颜色符号，请参照相关要求
end
% plot3(PointsAll(:,1),PointsAll(:,2),PointsAll(:,3),'-b.','LineWidth',2) %绘制路径，颜色符号，请参照相关要求
hold off
% plot3(dataset_1(SA,2),dataset_1(SA,3),dataset_1(SA,4),'-b.','LineWidth',3) %绘制路径，颜色符号，请参照相关要求

figure (3),
scatter(horizontal_check(:,2)',horizontal_check(:,3)','b','k')
hold on
scatter(vertical_check(:,2)',vertical_check(:,3)','g','k')
hold on
scatter(end_point(:,2)',end_point(:,3)','ro','MarkerFaceColor','c')
title('original data scatter','FontName','Times New Roman','FontWeight','Bold','FontSize',16)
xlabel('x','FontName','Times New Roman','FontSize',14)
ylabel('y','FontName','Times New Roman','FontSize',14,'Rotation',0)
ylabel('y','FontName','Times New Roman','FontSize',14,'Rotation',0)
hold on 

plot(end_point(:,2),end_point(:,3),'-r.','LineWidth',2)

for k=1:RecNo-1
    plot(PointsAll(2*k-1:2*k,1),PointsAll(2*k-1:2*k,2),'-r.','LineWidth',2) %绘制路径，颜色符号，请参照相关要求
    plot(PointsAll(2*k:2*k+1,1),PointsAll(2*k:2*k+1,2),'-b.','LineWidth',2) %绘制路径，颜色符号，请参照相关要求
end
% plot3(PointsAll(:,1),PointsAll(:,2),PointsAll(:,3),'-b.','LineWidth',2) %绘制路径，颜色符号，请参照相关要求
hold off

figure(5),
plot(end_point(:,2),end_point(:,3),'-r.','LineWidth',2)
hold on
for k=1:RecNo-1
    plot(PointsAll(2*k-1:2*k,1),PointsAll(2*k-1:2*k,2),'-r.','LineWidth',2) %绘制路径，颜色符号，请参照相关要求
    plot(PointsAll(2*k:2*k+1,1),PointsAll(2*k:2*k+1,2),'-b.','LineWidth',2) %绘制路径，颜色符号，请参照相关要求
end
hold off
