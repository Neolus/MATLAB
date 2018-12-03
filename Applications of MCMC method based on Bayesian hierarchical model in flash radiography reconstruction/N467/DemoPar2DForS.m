% Demo of linear MCMC reconstruction of blurred and noised data with SNR=46.7 for parallel beams.
% Note: the PUBLIC folder must be added to the MATLAB path.
% Copyright by @Neolus and @Mathwork Co.
% Data£º2018/06


%% Stimulate the projection
% You can skip this section if you could replace the right arguments of projection process.

% r1=1;r2=4.5;r3=6.5; %FTO parameters.
% u1=0;u2=0.0399*18.25;u3=0.0332*8.9;


dtc=0.04; % The length of a detector.
H=dtc;
sz=round(2*(r3+0.5)/H);
nd=sz; % nd=350;  %The number of detectors.
ng=nd/2;

SM=CmpSysMatCrcPar(H,dtc,ng,nd); %Compute the system matrix of a horizoontal slice.

FTO=load('FTO.mat'); %Load density data of FTO
FTO=FTO.FTO;
xt=FTO(:,ng+1:2*ng);
xt=reshape(xt',nd*ng,1);

% OPL=SM2D*xt; %Ideal optical depth data
% OPL=reshape(OPL,nd,nd)';


OPLBNG=load('FTOBNG467.mat'); %Load Optical depth data with SNR=46.7
OPLBNG=OPLBNG.OPLBNG;

% OLPM=OPLBNG-OPL;
% SNR=sqrt(sum(sum(OPL.^2))/sum(sum(OLPM.^2)+eps))

%% MCMC sample
%If you have the system matrix and the OPL data, you could start at here.
tStart=tic;

OPLData=OPLBNG;  %Choose the data for reconstruction
np=nd/2;      %The number of unknowns
SM=SM;        %the system matrix of a horizoontal slice

dg=8*ones(np,1);
dr=-2*ones(np,1);
Q=spdiags([dr dg dr], -1:1, np, np);  %Procesion matrix
QI=Q;

Nmax=6000;  %Number of samples
NChain=3;  %Number of Markov chains.
ATA=SM'*SM; 

noisem=zeros(1,np);
NDsc=Nmax/5;  %Discard the fisrt 1/5 samples
ndcmax=nd/2;  %Compute a quadrant

xk2D=zeros(np,ndcmax);

QALL=zeros(np,np,ndcmax);

for k=1:ndcmax
    QALL(:,:,k)=Q;
end

XCM=zeros(np,NChain);
XCALL=zeros(np,Nmax,NChain);  %Save samples of the center line.

d=zeros(1,Nmax);
l=zeros(1,Nmax);

m=size(SM,1);
n=size(SM,2);
nfr=np+1;

XCS=zeros(nd/2,np,NChain);
for i=1:NChain  % For each Markov chains
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
    for k=1:Nmax  %For each sample
        for j=1:ndcmax  % For each pixel
            y=OPLData(j,:)';
            lamda=lamall(j);  %Load parameters of Markov chains
            deta=betall(j);
            Q=QALL(:,:,j);
            
            Q12=chol(Q,'lower');
            A=lambda*ATA+deta*Q;
            
            v1=randn(1,m)';
            v2=randn(1,n)';
            
            noise=sqrt(lambda)*SM'*v1+sqrt(deta)*Q12*v2; %A quick technique
            yk=lambda*SM'*y+noise;
            L=chol(A,'lower'); %Cholesky factorition.
            
            v=L\yk; %Compute the sample of X
            xk=L'\v;
            xk2D(:,j)=xk;
        end
%         if k<NDsc
            xk2D=medfilt2(xk2D,[5,5],'symmetric');
%         end
        
%         if k>NDsc
%         xkbr=xk2D(1:15,:);
%         xkbr=medfilt2(xkbr,[5,5],'symmetric');
%         xk2D(1:3,:)=xkbr(1:3,:);
%         end
        

        
        for j=1:ndcmax
            xk=xk2D(:,j);
            xk=max(0,xk);  %Non-negative constraint
            
            deta=gamrnd(np/2+alpdt/2,0.5*norm((SM*xk-y),2)+betadt);  % Sampling \delta
            lambda=gamrnd(np/2+alplm/2,0.5*xk'*Q*xk+betalm); %Sampling \lambda
            
            WM=QI-deta/2*QI*xk*(xk')*QI/(1+deta/2*(xk')*QI*xk);  %Smapling the covariance/procesion matrix            
            WMC=chol(WM);
            wisha=randn(nfr,np)*WMC;            
            % wisha=mvnrnd(noisem,WM,np+1);
            Q=wisha'*wisha;
            
            QALL(:,:,j)=Q;
            lamall(j)=lamda;
            betall(j)=deta;
            
            
            if k>NDsc
                if(j==nd/2)  %Save samples of the center line.
                    XCALL(:,k,i)=xk;
                end
            XCS(j,:,i)= XCS(j,:,i)+xk';               
            end
        end
    end
   
end

XMC= XCS/(Nmax-NDsc);
XMC=mean(XMC,3); %MCMC result.

%% Analysis of the MCMC result.

XCMPF=[fliplr(XMC),XMC]; %Flip a quadrant to four quadrants according to the symmetry of FTO
XCMPF=[XCMPF;flipud(XCMPF)];  %MCMC reconstruction 2D

save('N467C3N6000.mat','XCMPF') 
save('N467C3N600CSALL.mat','XCALL') 


r=sqrt(sum(sum((FTO-XCMPF).^2))/sum(sum(FTO.^2))) %RMSE

XSamC=XCALL(:,NDsc+1:Nmax,:);
XSamC=reshape(XSamC,[np,(Nmax-NDsc)*NChain]);
CI=quantile(XSamC,[0.025,0.975],2); %Quantile 95% the credible interval

XSM=mean(XSamC,2);
save('N467C3N6000CICS.mat','CI');

tEnd=toc(tStart);
RMSE=[r,tEnd];
save('RMSEN0.mat','RMSE');

figure(12),plot(FTO(nd/2,nd/2+1:nd),'b--');
hold on, plot(XSM,'r');
% plot(xrc(nd/2,nd/2+1:nd))
plot(CI(:,1),'g');
plot(CI(:,2),'g');
hold off
legend('real data','MCMC mean','95% credible interval')
xlabel('location/pixels','Fontname','Times New Roman');
ylabel('linear attenuation coefficient/cm^{-1}','Fontname','Times New Roman')
set(gca,'fontname','Times New Roman')
