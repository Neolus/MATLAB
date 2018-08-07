function coef=Halir(samples,BW,nh)
%HALIR Fits ellipses with a least square method modified by Halir.
%  COEF=HALIR(SAMPLES,BW,NH) fits ellipses based SMAPLES which is
%  a result of Function:CMHTN or CMPN. 
%  BW the binary image. 
%  NH the neighbour-hood of samples. 
%  COEF are coefficients of fitted ellispes.
%  SIZE(COEF)=[NE,6]. NE is the number of ellispes.
%  Every ellipse have 6 coefficients. See line 35.
if nargin<3
    nh=27; % Set the default neighbour-hood
end
[M,N,ne]=size(samples); % NE is the number of ellipses
%Constrain matrix, C1=[0,0,2; 0,-1,0; 2,0,0]; Cinv=inv(C1)
Cinv=[0, 0, 0.5; 0, -1, 0; 0.5, 0, 0];
coef=zeros(ne,6);
hold on
for k=1:ne
    XY=[];
    for i=1:M
        san=neighbour(samples(i,1,k),samples(i,2,k),nh,BW);
        XY=[XY;san];
    end
    N=length(XY);
    D1=[XY(:,1).^2, XY(:,1).*XY(:,2),XY(:,2).^2];  D2=[XY,ones(N,1)];
    %Design matrix D=[D1, D2]; %Scatter matrix S=D' *D=[S1, S2; S2', S3]
    S1=D1'*D1;  S2=D1'*D2;  S3=D2'*D2;
    Mat=Cinv*(S1-S2/S3*S2');
    [V,D]=eig(Mat,'nobalance');   %Eigenvector and  Eigenvalue
    %Find the smallest positive egienvalue
    D_diag=diag(D);  D_diag(D_diag<=0)=inf;  [lamda, idx]=min(D_diag);
    alpha1=V(:,idx); alpha2=-inv(S3)*S2'*alpha1;
    %The general equation of ellipse: Ax2+Bxy+Cy2+Dx+Ey+F=0 F=@(p,x)p(1)*x(:,1).^2+p(2)*x(:,1).*x(:,2)+p(3)*x(:,2).^2+p(4)*x(:,1)+p(5)*x(:,2)+p(6);
    pr1=[alpha1',alpha2']/alpha2(3); %Set F to 1
    coef(k,:)=pr1;
    xmin=0;  xmax=550;  ymin=0;  ymax=500;
    hold on;
    h=ezplot(@(x,y)F(pr1,[y,x]),[-1+xmin,1+xmax,-1+ymin,1+ymax]); set(h,'Color','r');
end
end

function XY=neighbour(x,y,n,f) % Return points whose value is 1 in n*n neighbour-hoods of f(x,y)
[M,N]=size(f);
n1=floor(n/2); % n*n neighbour-hood
if x-n1<1
    x1=1:n;
elseif x+n1>M
    x1=(M-n+1):M;
else
    x1=x-n1:x+n1;
end
if y-n1<1
    y1=1:n;
elseif y+n1>N
    y1=(N-n+1):N;
else
    y1=y-n1:y+n1;
end
fc=f(x1,y1);  [XY(:,1),XY(:,2)]=find(fc~=0);  % m*1
xc=find(x1==x);   yc=find(y1==y); %Find (x,y) in the new matrix
XY(:,1)=XY(:,1)+x-xc;  XY(:,2)=XY(:,2)+y-yc %Back to image coordinate.;
end
