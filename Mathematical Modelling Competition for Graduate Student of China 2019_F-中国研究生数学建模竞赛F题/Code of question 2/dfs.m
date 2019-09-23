function dfs(Ver,Pre,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC)
%%
% VER当前顶点
% PRE上一个顶点

%全局变量 
global disMat2; %距离矩阵
global ResFlag; %校正点类型
global V;       %飞行器速度
global TR;      %一次校正所需时间
global fobj;    %目标函数值
global RootABG; %最优路径
global PointsCG; %最优路径对应的切点
global Position;

VerNo=Ver.No; %当前点序号
% if VerNo==1
%     disP=0;
% else
%     disP=disMat2(VerNo,Pre.No);   
% end

VisFlag(Ver.No)=1; %标记访问
Ver.PreNo=Pre.No;  %上一个点的序号

% disP=disMat2(VerNo,Pre.No); %两点间的距离
disP=Ver.DisPre2Now;%上一点到当前点的距离


ErrorIncr=disP*delta; %增加的误差
Ver.Dis=Pre.Dis+disP; %距离增加

Ver.NR=Pre.NR+1; %路径长度加1
Ver.Pos=Position(VerNo,:); %坐标值

RootAB(Ver.NR+1)=VerNo; %保存路径



% Ver.Obj=Ver.Dis/V+Ver.NR*TR; %距离和次数的目标函数，建议也选用这个，但还没有测试，请你们自己实验。
Ver.Obj=Ver.Dis; %取距离为目标函数

Ver.HError=Pre.HError+ErrorIncr; %水平误差增加
Ver.VError=Pre.VError+ErrorIncr; %垂直误差增加

if (Ver.HError>=theta||Ver.VError>=theta) %大于误差直接返回
    return;
end


%误差校正
if ResFlag(VerNo)==1 %1表示垂直误差校正点
    if(Ver.HError<=RV_HError&&Ver.VError<=RV_VError)
        Ver.VError=0;
    end
end

if ResFlag(VerNo)==0 %0表示水平误差校正点
    if(Ver.HError<=RH_HError&&Ver.VError<=RH_VErrror)
        Ver.HError=0;
    end
end

if Ver.Obj>=fobj %到达该点的目标函数值已经大于当前最优目标函数值
    fprintf(sprintf('Revise points %d, Objection value of vertex: %f\n',Ver.NR,Ver.Obj));
    return;
end

if(VerNo==m) %到达终点
    if(Ver.Obj<fobj) %目标函数值变小才更新
        fobj=Ver.Obj;
        fprintf(sprintf('Revise points %d, Object function value: %f\n',Ver.NR,fobj));
        RootABG=RootAB;
        RootABG(1:Ver.NR+1) %界面显示，建议写成Txt记录，包括目标函数值及路径长度，用于写论文
        PointsCG=PointsC;
        PointsCG=PointsCG(1:Ver.NR+1,:);
    end
    return;
end

ErrorMax=theta-max(Ver.HError,Ver.VError); %
DisMax=ErrorMax/delta; %能再次飞行的最远距离

DisVer=disMat2(VerNo,:);
DisToB=disMat2(VerNo,m); %直接到B的目标函数增量，还有优化的余地，如最短路径等
% DisToB=disMat2(VerNo,m)/V; %距离加次数模型，选用这个，或者设计其它的

if((Ver.Obj+DisToB)>=fobj) %如果当前点的目标函数值加上直接到终点的已经大于，已有目标函数值，剪枝，停止向下搜索，非常关键！！！！
%     fprintf(sprintf('Revise points %d, Objection value to B: %f, Object function value%f\n',Ver.NR,Ver.Obj+DisToB,fobj));
    return;
end

AdjVers=find(DisVer<DisMax&VisFlag==0); %可到达的且没有被访问过的点，粗略计算直线能够飞行的最远距离的点
AdjVersN=length(AdjVers);

% AdjVersEff=zeros(AdjVersN,1); %有效的邻点
% PointCEff=zeros(AdjVersN,3);
% DisNEff=zeros(AdjVersN,1);
% cntAdjEff=0;
% for j=1:AdjVersN  %根据飞行距离再次筛选
%        PosNext=Position(AdjVers(j),:);
%        [DisFly,PointC]=ComputeSphDis(Ver.Pos,PosNext,Ver.Velo,r);
%        if DisFly<DisMax&&DisFly>0
%            cntAdjEff=cntAdjEff+1;
%            AdjVersEff(cntAdjEff)=AdjVers(j);
%            PointCEff(cntAdjEff,:)=PointC;
%            DisNEff(cntAdjEff)=DisFly;
%        end  
%     
% end
% AdjVersEff=AdjVersEff(1:cntAdjEff);
% AdjVersN=length(AdjVersEff);
% PointCEff=PointCEff(1:cntAdjEff,:);
% DisNEff=DisNEff(1:cntAdjEff);
% fprintf(sprintf('Revise points %d, Adjacent vertex number: %d, distance %f\n',Ver.NR,AdjVersN,Ver.Dis));

MaxAdjNo=max(AdjVers);
if(MaxAdjNo==m) %最大点在相邻点内，只计算这一种情况
    PosNext=Position(m,:);
    Velocity=PosNext-Ver.Pos; %飞往下一点的速度方向
    [DisFly,PointC]=ComputeSphDis(Ver.Pos,PosNext,Ver.Velo,r);    
    DisNext=DisFly; %到下一个点的距离
    PointsC(Ver.NR+1,:)=PointC; %保存切点
    
    VerNext= struct('No',{MaxAdjNo},'PreNo',{VerNo},'Dis',{inf},'NR',{Ver.NR+1},'Obj',{inf},'HError',{inf},'VError',{inf},'Pos',{PosNext},'Velo',{Velocity},'DisPre2Now',DisNext);
    dfs(VerNext,Ver,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC);  %递归
    %     VisFlag(AdjVers(k))=0; %删除访问标记
    return;
else
    for k=1:AdjVersN
        PosNext=Position(AdjVers(k),:);
        Velocity=PosNext-Ver.Pos; %飞往下一点的速度方向
        [DisFly,PointC]=ComputeSphDis(Ver.Pos,PosNext,Ver.Velo,r);
        DisNext=DisFly; %到下一个点的距离
        PointsC(Ver.NR+1,:)=PointC; %保存切点
        VerNext= struct('No',{AdjVers(k)},'PreNo',{VerNo},'Dis',{inf},'NR',{Ver.NR+1},'Obj',{inf},'HError',{inf},'VError',{inf},'Pos',{PosNext},'Velo',{Velocity},'DisPre2Now',DisNext);
        dfs(VerNext,Ver,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC); %递归
%         VisFlag(AdjVers(k))=0; %删除访问标记
%        return;
    end
    
end

% VisFlag(Ver.No)=0;  %删除访问标记
% return;

end