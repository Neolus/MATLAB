function dfsDrive(Ver,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC)
%%
%VER��ǰ����
%M�յ����
%THETA��DELTA���в���
%RV_VError,RV_HError,RH_VErrror,RH_HErrorУ������
%RootAB��¼·��%
%r��뾶
%PointsC��¼�е�����

global disMat2; %�������
global Position;

VerNo=Ver.No;  %��ʼ��A

VisFlag(Ver.No)=1; 
Ver.PreNo=0;
% disP=0;
% ErrorIncr=0;
Ver.Dis=0;

Ver.NR=0; 
RootAB(Ver.NR+1)=VerNo; %����·��
PointsC(Ver.NR+1,:)=Position(VerNo); %�����е�
Ver.Obj=0;

Ver.HError=0; %ˮƽ���
Ver.VError=0; %��ֱ���
Ver.Pos=Position(VerNo,:); %����ֵ

ErrorMax=theta-max(Ver.HError,Ver.VError);
DisMax=ErrorMax/delta;

DisVer=disMat2(VerNo,:); %��ǰ����������ľ���

AdjVers=find(DisVer<DisMax&VisFlag==0); %�ɵ������û�б����ʹ��ĵ�
AdjVersN=length(AdjVers);

MaxAdjNo=max(AdjVers);
if(MaxAdjNo==m) %���������ڵ��ڣ�ֻ������һ�����
    PosNext=Position(MaxAdjNo,:);
    Velocity=PosNext-Ver.Pos; %������һ����ٶȷ���
    DisNext=sqrt(sum(Velocity.^2));%����һ����ľ���
    VerNext= struct('No',{MaxAdjNo},'PreNo',{VerNo},'Dis',{inf},'NR',{Ver.NR+1},'Obj',{inf},'HError',{inf},'VError',{inf},'Pos',{PosNext},'Velo',{Velocity},'DisPre2Now',DisNext);
    dfs(VerNext,Ver,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC);
%     VisFlag(AdjVers(k))=0; %ɾ�����ʱ��
    return;
else
    for k=1:AdjVersN %��������A�����ڵ�
        fprintf(sprintf('Search of the %d leaf in %d leafs\n',k,AdjVersN));
        PosNext=Position(AdjVers(k),:);
        Velocity=PosNext-Ver.Pos; %������һ����ٶȷ���
        DisNext=sqrt(sum(Velocity.^2)); %����һ����ľ���
        VerNext= struct('No',{AdjVers(k)},'PreNo',{VerNo},'Dis',{inf},'NR',{Ver.NR+1},'Obj',{inf},'HError',{inf},'VError',{inf},'Pos',{PosNext},'Velo',{Velocity},'DisPre2Now',DisNext);
        dfs(VerNext,Ver,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC);
%         VisFlag(AdjVers(k))=0; %ɾ�����ʱ��
%         return;
    end
    
end

% VisFlag(Ver.No)=0;  %ɾ�����ʱ��
return;

end