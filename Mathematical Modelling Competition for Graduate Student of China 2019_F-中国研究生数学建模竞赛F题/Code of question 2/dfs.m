function dfs(Ver,Pre,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC)
%%
% VER��ǰ����
% PRE��һ������

%ȫ�ֱ��� 
global disMat2; %�������
global ResFlag; %У��������
global V;       %�������ٶ�
global TR;      %һ��У������ʱ��
global fobj;    %Ŀ�꺯��ֵ
global RootABG; %����·��
global PointsCG; %����·����Ӧ���е�
global Position;

VerNo=Ver.No; %��ǰ�����
% if VerNo==1
%     disP=0;
% else
%     disP=disMat2(VerNo,Pre.No);   
% end

VisFlag(Ver.No)=1; %��Ƿ���
Ver.PreNo=Pre.No;  %��һ��������

% disP=disMat2(VerNo,Pre.No); %�����ľ���
disP=Ver.DisPre2Now;%��һ�㵽��ǰ��ľ���


ErrorIncr=disP*delta; %���ӵ����
Ver.Dis=Pre.Dis+disP; %��������

Ver.NR=Pre.NR+1; %·�����ȼ�1
Ver.Pos=Position(VerNo,:); %����ֵ

RootAB(Ver.NR+1)=VerNo; %����·��



% Ver.Obj=Ver.Dis/V+Ver.NR*TR; %����ʹ�����Ŀ�꺯��������Ҳѡ�����������û�в��ԣ��������Լ�ʵ�顣
Ver.Obj=Ver.Dis; %ȡ����ΪĿ�꺯��

Ver.HError=Pre.HError+ErrorIncr; %ˮƽ�������
Ver.VError=Pre.VError+ErrorIncr; %��ֱ�������

if (Ver.HError>=theta||Ver.VError>=theta) %�������ֱ�ӷ���
    return;
end


%���У��
if ResFlag(VerNo)==1 %1��ʾ��ֱ���У����
    if(Ver.HError<=RV_HError&&Ver.VError<=RV_VError)
        Ver.VError=0;
    end
end

if ResFlag(VerNo)==0 %0��ʾˮƽ���У����
    if(Ver.HError<=RH_HError&&Ver.VError<=RH_VErrror)
        Ver.HError=0;
    end
end

if Ver.Obj>=fobj %����õ��Ŀ�꺯��ֵ�Ѿ����ڵ�ǰ����Ŀ�꺯��ֵ
    fprintf(sprintf('Revise points %d, Objection value of vertex: %f\n',Ver.NR,Ver.Obj));
    return;
end

if(VerNo==m) %�����յ�
    if(Ver.Obj<fobj) %Ŀ�꺯��ֵ��С�Ÿ���
        fobj=Ver.Obj;
        fprintf(sprintf('Revise points %d, Object function value: %f\n',Ver.NR,fobj));
        RootABG=RootAB;
        RootABG(1:Ver.NR+1) %������ʾ������д��Txt��¼������Ŀ�꺯��ֵ��·�����ȣ�����д����
        PointsCG=PointsC;
        PointsCG=PointsCG(1:Ver.NR+1,:);
    end
    return;
end

ErrorMax=theta-max(Ver.HError,Ver.VError); %
DisMax=ErrorMax/delta; %���ٴη��е���Զ����

DisVer=disMat2(VerNo,:);
DisToB=disMat2(VerNo,m); %ֱ�ӵ�B��Ŀ�꺯�������������Ż�����أ������·����
% DisToB=disMat2(VerNo,m)/V; %����Ӵ���ģ�ͣ�ѡ��������������������

if((Ver.Obj+DisToB)>=fobj) %�����ǰ���Ŀ�꺯��ֵ����ֱ�ӵ��յ���Ѿ����ڣ�����Ŀ�꺯��ֵ����֦��ֹͣ�����������ǳ��ؼ���������
%     fprintf(sprintf('Revise points %d, Objection value to B: %f, Object function value%f\n',Ver.NR,Ver.Obj+DisToB,fobj));
    return;
end

AdjVers=find(DisVer<DisMax&VisFlag==0); %�ɵ������û�б����ʹ��ĵ㣬���Լ���ֱ���ܹ����е���Զ����ĵ�
AdjVersN=length(AdjVers);

% AdjVersEff=zeros(AdjVersN,1); %��Ч���ڵ�
% PointCEff=zeros(AdjVersN,3);
% DisNEff=zeros(AdjVersN,1);
% cntAdjEff=0;
% for j=1:AdjVersN  %���ݷ��о����ٴ�ɸѡ
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
if(MaxAdjNo==m) %���������ڵ��ڣ�ֻ������һ�����
    PosNext=Position(m,:);
    Velocity=PosNext-Ver.Pos; %������һ����ٶȷ���
    [DisFly,PointC]=ComputeSphDis(Ver.Pos,PosNext,Ver.Velo,r);    
    DisNext=DisFly; %����һ����ľ���
    PointsC(Ver.NR+1,:)=PointC; %�����е�
    
    VerNext= struct('No',{MaxAdjNo},'PreNo',{VerNo},'Dis',{inf},'NR',{Ver.NR+1},'Obj',{inf},'HError',{inf},'VError',{inf},'Pos',{PosNext},'Velo',{Velocity},'DisPre2Now',DisNext);
    dfs(VerNext,Ver,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC);  %�ݹ�
    %     VisFlag(AdjVers(k))=0; %ɾ�����ʱ��
    return;
else
    for k=1:AdjVersN
        PosNext=Position(AdjVers(k),:);
        Velocity=PosNext-Ver.Pos; %������һ����ٶȷ���
        [DisFly,PointC]=ComputeSphDis(Ver.Pos,PosNext,Ver.Velo,r);
        DisNext=DisFly; %����һ����ľ���
        PointsC(Ver.NR+1,:)=PointC; %�����е�
        VerNext= struct('No',{AdjVers(k)},'PreNo',{VerNo},'Dis',{inf},'NR',{Ver.NR+1},'Obj',{inf},'HError',{inf},'VError',{inf},'Pos',{PosNext},'Velo',{Velocity},'DisPre2Now',DisNext);
        dfs(VerNext,Ver,m,theta,delta,VisFlag,RV_VError,RV_HError,RH_VErrror,RH_HError,RootAB,r,PointsC); %�ݹ�
%         VisFlag(AdjVers(k))=0; %ɾ�����ʱ��
%        return;
    end
    
end

% VisFlag(Ver.No)=0;  %ɾ�����ʱ��
% return;

end