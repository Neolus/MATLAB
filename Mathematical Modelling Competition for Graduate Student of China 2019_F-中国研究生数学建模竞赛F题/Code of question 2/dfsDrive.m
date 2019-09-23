function dfsDrive(Ver,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC)
%%
%VER当前顶点
%M终点序号
%THETA，DELTA飞行参数
%RV_VError,RV_HError,RH_VErrror,RH_HError校正参数
%RootAB记录路径%
%r球半径
%PointsC记录切点坐标

global disMat2; %距离矩阵
global Position;

VerNo=Ver.No;  %初始化A

VisFlag(Ver.No)=1; 
Ver.PreNo=0;
% disP=0;
% ErrorIncr=0;
Ver.Dis=0;

Ver.NR=0; 
RootAB(Ver.NR+1)=VerNo; %保存路径
PointsC(Ver.NR+1,:)=Position(VerNo); %保存切点
Ver.Obj=0;

Ver.HError=0; %水平误差
Ver.VError=0; %垂直误差
Ver.Pos=Position(VerNo,:); %坐标值

ErrorMax=theta-max(Ver.HError,Ver.VError);
DisMax=ErrorMax/delta;

DisVer=disMat2(VerNo,:); %当前点的与其余点的距离

AdjVers=find(DisVer<DisMax&VisFlag==0); %可到达的且没有被访问过的点
AdjVersN=length(AdjVers);

MaxAdjNo=max(AdjVers);
if(MaxAdjNo==m) %最大点在相邻点内，只计算这一种情况
    PosNext=Position(MaxAdjNo,:);
    Velocity=PosNext-Ver.Pos; %飞往下一点的速度方向
    DisNext=sqrt(sum(Velocity.^2));%到下一个点的距离
    VerNext= struct('No',{MaxAdjNo},'PreNo',{VerNo},'Dis',{inf},'NR',{Ver.NR+1},'Obj',{inf},'HError',{inf},'VError',{inf},'Pos',{PosNext},'Velo',{Velocity},'DisPre2Now',DisNext);
    dfs(VerNext,Ver,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC);
%     VisFlag(AdjVers(k))=0; %删除访问标记
    return;
else
    for k=1:AdjVersN %搜索所有A的相邻点
        fprintf(sprintf('Search of the %d leaf in %d leafs\n',k,AdjVersN));
        PosNext=Position(AdjVers(k),:);
        Velocity=PosNext-Ver.Pos; %飞往下一点的速度方向
        DisNext=sqrt(sum(Velocity.^2)); %到下一个点的距离
        VerNext= struct('No',{AdjVers(k)},'PreNo',{VerNo},'Dis',{inf},'NR',{Ver.NR+1},'Obj',{inf},'HError',{inf},'VError',{inf},'Pos',{PosNext},'Velo',{Velocity},'DisPre2Now',DisNext);
        dfs(VerNext,Ver,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC);
%         VisFlag(AdjVers(k))=0; %删除访问标记
%         return;
    end
    
end

% VisFlag(Ver.No)=0;  %删除访问标记
return;

end