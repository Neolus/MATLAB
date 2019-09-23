function [DisFly,PointC]=ComputeSphDis(Sta,Des,Vel,r)
%% 计算当前校正点到下一校正点的飞行路径长度
% input: STA 起点坐标
%        DES 终点坐标
%        VEL起点的运动方向
%        R球半径
%output:DISFLY;飞行路径长度=球面弧长+脱离点到终点的直线距离
%       POINTC;脱离点的坐标

% r球半径
PointC=inf*ones(1,3);

x1=Sta(1); %起点坐标A
y1=Sta(2);
z1=Sta(3);

x3=Des(1);%终点坐标B
y3=Des(2);
z3=Des(3);

nx=Vel(1);%起点法向量
ny=Vel(2);
nz=Vel(3);

%%计算终点在球心平面的投影  
%求解方程组：nx(x-x1)+ny(y-y1)+nz(z-z1)==0和(x-x3)/nx=(y-y3)/ny=(z-z3)/nz=lambda
lamb=-(nx*(x3-x1)+ny*(y3-y1)+nz*(z3-z1))/(nx*nx+ny*ny+nz*nz+eps);
x2=lamb*nx+x3;
y2=lamb*ny+y3;
z2=lamb*nz+z3;

%%求解球心坐标
lambr=r/sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
xo1=x1+lambr*(x2-x1); %两组
yo1=y1+lambr*(y2-y1);
zo1=z1+lambr*(z2-z1);

xo2=x1-lambr*(x2-x1);
yo2=y1-lambr*(y2-y1);
zo2=z1-lambr*(z2-z1);


OB1=(xo1-x3)^2+(yo1-y3)^2+(zo1-z3)^2;
OB2=(xo2-x3)^2+(yo2-y3)^2+(zo2-z3)^2;

if OB1<OB2 %取到终点距离小的点为球心O
    Xo=xo1; 
    Yo=yo1;
    Zo=zo1;
else
    Xo=xo2; 
    Yo=yo2;
    Zo=zo2;   
end

%记代求切点为C,记切线所在圆环的圆心为D
DisOB=sqrt((Xo-x3)^2+(Yo-y3)^2+(Zo-z3)^2);
DisBC=sqrt(DisOB^2-r^2);
DisBD=DisBC^2/DisOB;

lambD=DisBD/DisOB;

xd=x3+lambD*(Xo-x3);
yd=y3+lambD*(Yo-y3);
zd=z3+lambD*(Zo-z3);

%计算A到切点所在平面的投影坐标E
mx=Xo-x3;
my=Yo-y3;
mz=Zo-z3;

lame=-(mx*(x1-xd)+my*(y1-yd)+mz*(z1-zd))/(mx*mx+my*my+mz*mz+eps);
xe=lame*mx+x1;
ye=lame*my+y1;
ze=lame*mz+z1;


%计算脱离点C
a=(xe-xd)^2+(ye-yd)^2+(ze-zd)^2;
b=2*((xe-xd)*(xd-Xo)+(ye-yd)*(yd-Yo)+(ze-zd)*(zd-Zo));
c=(xd-Xo)^2+(yd-Yo)^2+(zd-Zo)^2-r^2;

if b^2-4*a*c<=0
    DisFly=inf;
    return
end

lam1=(-b+sqrt(b^2-4*a*c))/(2*a);
lam2=(-b-sqrt(b^2-4*a*c))/(2*a);


xc1=xd+lam1*(xe-xd);
yc1=yd+lam1*(ye-yd);
zc1=zd+lam1*(ze-zd);

disAC1=sqrt((xc1-x1)^2+(yc1-y1)^2+(zc1-z1)^2);

xc2=xd+lam2*(xe-xd);
yc2=yd+lam2*(ye-yd);
zc2=zd+lam2*(ze-zd);
disAC2=sqrt((xc2-x1)^2+(yc2-y1)^2+(zc2-z1)^2);

if(disAC1<disAC2)
    XC=xc1;
    YC=yc1;
    ZC=zc1;
    DisAC=disAC1;
else
    XC=xc2;
    YC=yc2;
    ZC=zc2;   
    DisAC=disAC2;
end
PointC(1)=XC;
PointC(2)=YC;
PointC(3)=ZC;

radAC=2*asin(DisAC/(2*r));

DisFly=radAC*r+DisBC; 





end