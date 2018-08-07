% Demo of linear MCMC reconstruction for parallel beams.
% Copyright by @Neolus and @Mathwork Co.
% Data£º2018/06

%% Stimulate the projection
% You can skip this section if you could replace the right arguments of projection process.
r1=1;r2=4.5;r3=6.5; %FTO parameters.
u1=0;u2=0.0399*18.25;u3=0.0332*8.9;


dtc=0.04; % The length of a detector.
% nd=512; dtc=0.05; %The number of detectors.
H=dtc;
sz=round(2*(r3+0.5)/H);
nd=sz; % nd=350;  %The number of detectors.
ng=nd/2;

SM=CmpSysMatCrcPar(H,dtc,ng,nd); %Compute the system matrix of a horizoontal slice.
SM2D=CmpSysMatCrcDiag(SM,nd);% System matrix of parallel beams.

FTO=GetFTO(r1,r2,r3,u1,u2,u3,sz,H); %FTO density.
xt=FTO(:,ng+1:2*ng);
xt=reshape(xt',nd*ng,1);

OPL=SM2D*xt;
OPL=reshape(OPL,nd,nd)';
% figure(11),imshow(OPL,[]);

% sig=1.5;
% hz=2*ceil(2.6*sig)+1;   %  Truncation is less than 1%.
% GLPF=fspecial('gaussian',[hz,hz],sig); %Gaussian low pass filter.
% 
% OPLB=conv2(OPL,GLPF,'same'); %Add a Gaussian blur to the projections.
% % figure(12),imshow(OPLB,[]);
% 
% % varn=(0.01635*std(OPLB(:)))^2; %SNR 12.6
% % varn=(0.006145*std(OPLB(:)))^2; %SNR 32.3
% varn=(0.00410*std(OPLB(:)))^2; %SNR 46.7
% % varn=0;
% OPLBNG=imnoise(OPLB/max(OPLB(:)),'gaussian',0,varn);
% OPLBNG=max(OPLB(:))*OPLBNG;
% figure,imshow(OPLBNG,[])

% OPLBNG=load('FTOBNG467.mat');
% OPLBNG=OPLBNG.OPLBNG;

% figure,plot(OPL(nd/2,:));
% hold on,plot(OPLB(nd/2,:),'r')
% plot(OPLBNG(nd/2,:),'g')
% hold off
% OLPM=OPLBNG-OPL;
% SNR=sqrt(sum(sum(OPL.^2))/sum(sum(OLPM.^2)+eps))

%% MCMC sample
%If you have the system matrix and the OPL data, you could start at here.
tStart=tic;

OPLData=OPL;  %Choose the data for reconstruction
np=nd/2;      %The number of unknowns
SM=SM;        %the system matrix of a horizoontal slice

dg=8*ones(np,1);
dr=-2*ones(np,1);
Q=spdiags([dr dg dr], -1:1, np, np);  %Procesion matrix
%  Q(1,1)=2; Q(np,np)=2;
QI=Q;

Nmax=6000;  %Number of samples
NChain=3;  %Number of Markov chains.
ATA=SM'*SM; 

noisem=zeros(1,np);
NDsc=Nmax/5;  %Discard the fisrt 1/5 samples
ndcmax=nd/2;  %Compute a quadrant

xk2D=zeros(np,ndcmax); %Reconstructed image in each iteration

QALL=zeros(np,np,ndcmax);

for k=1:ndcmax
    QALL(:,:,k)=Q;
end

XCM=zeros(np,NChain);
XCALL=zeros(np,Nmax,NChain);  %Keep samples of the center line.

d=zeros(1,Nmax);
l=zeros(1,Nmax);

m=size(SM,1);
n=size(SM,2);
nfr=np+1;

XCS=zeros(nd/2,np,NChain);
for i=1:NChain
    prm=rand(1,6);  %Initialize arguments randomly.
    lambda= prm(1);
    deta= prm(2);
    alpdt=prm(3);
    alplm=prm(4);
    betadt=prm(5);
    betalm=prm(6);
    
    lamall=rand(nd,1);
    betall=rand(nd,1);
 
    i
    for k=1:Nmax  %Repeat sampling
        for j=1:ndcmax  %Compute image row by row
            y=OPLData(j,:)'; %OPL of j-th row
            lamda=lamall(j);
            deta=betall(j);
            Q=QALL(:,:,j);
            
            Q12=chol(Q,'lower');
            A=lambda*ATA+deta*Q;
            
            v1=randn(1,m)';
            v2=randn(1,n)';
            
            noise=sqrt(lambda)*SM'*v1+sqrt(deta)*Q12*v2;
            yk=lambda*SM'*y+noise;
            L=chol(A,'lower'); %Cholesky factorition.
            
            v=L\yk;
            xk=L'\v;  %Samples the j-th line image.
            xk2D(:,j)=xk; %Place samples in xk2D
        end
        if k<NDsc  %Median filter 
            xk2D=medfilt2(xk2D,[5,5],'symmetric');
%                     xkbr=xk2D(1:17,:);
%                     xkbr=medfilt2(xkbr,[9,9],'symmetric');
%                     xk2D(1:5,:)=xkbr(1:5,:);
        end
        
        if k>NDsc
        xkbr=xk2D(1:15,:);  %Filter only the vertical center area. 
        xkbr=medfilt2(xkbr,[5,5],'symmetric');
        xk2D(1:3,:)=xkbr(1:3,:);
        end
        

        
        for j=1:ndcmax
            xk=xk2D(:,j);
            xk=max(0,xk);
            
            deta=gamrnd(np/2+alpdt/2,0.5*norm((SM*xk-y),2)+betadt); %Sample DETA
            lambda=gamrnd(np/2+alplm/2,0.5*xk'*Q*xk+betalm);  %Sample LAMBDA
            
            WM=QI-deta/2*QI*xk*(xk')*QI/(1+deta/2*(xk')*QI*xk);
            
            WMC=chol(WM);
            wisha=randn(nfr,np)*WMC;            
            % wisha=mvnrnd(noisem,WM,np+1);
            Q=wisha'*wisha;
            QALL(:,:,j)=Q;
            lamall(j)=lamda;
            betall(j)=deta;

            
            if k>NDsc
                if(j==nd/2)  %Keep samples of the center line.
%                     XCALL(:,k,i)=xk;
                    XCALL(:,k,i)=xk;
                end
            XCS(j,:,i)= XCS(j,:,i)+xk';
                
             
%                 XCS(j,:,i)= XCS(j,:,i)+xk2';
            end
        end
    end
   
end

XMC= XCS/(Nmax-NDsc);
XMC=mean(XMC,3); %MCMC result.

%% Analysis of result.
%Flip a quadrant to an image, assuming that FTO is symmetry.
XCMPF=[fliplr(XMC),XMC];
XCMPF=[XCMPF;flipud(XCMPF)];

save('N0C3N6000.mat','XCMPF')
save('N0C3N600CSALL.mat','XCALL')

%Compute root-mean-square error
r=sqrt(sum(sum((FTO-XCMPF).^2))/sum(sum(FTO.^2)))

XALLD=XCALL(:,NDsc+1:Nmax,:);
% XALLD=max(0,XALLD);
XS=mean(XALLD,3);

XSM=mean(XS,2);
XSStd=std(XS,1,2);

CIL=XSM-1.96*XSStd; %Cridibilty interval
CIU=XSM+1.96*XSStd;
CI=[CIL,CIU];
save('N0C3N6000CICS.mat','CI');

tEnd=toc(tStart)
RMSE=[r,tEnd];
save('RMSEN0.mat','RMSE');

% figure(10),imshow(XCMPF,[]);
% 
% figure(11),plot(XMC(nd/2,:));
% hold on,plot(FTO(nd/2,nd/2+1:nd),'r');
% hold off
% 
% 
% figure(12),plot(XCMPF(120,:));
% hold on; plot(FTO(120,:));
% hold off
% 
% 
% figure(13),plot(XCMPF(80,:));
% hold on; plot(FTO(80,:));
% hold off
% 
% 
% 
% figure(14),plot(FTO(nd/2,nd/2+1:nd),'b--');
% hold on
% plot(XSM,'k');
% plot(CIL,'r');
% plot(CIU,'r');
% hold off
% % 
