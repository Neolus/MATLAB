function [SM,SMU]=CmpSysMatCrcPar(dg,dtc,nr,nd)
%CMPSYSMATSQR Compute the system matrix of parallel beams by circle grids.
% SM=CMPSYSMATSQR(DG,DTC,OS,OD,NG,ND) compute system matrix for one view. 
% DG is the length of a grid.
% DTC is the length of a detector.
% OS is the distance from the object center to the fan beam source.
% OD is the distance from the object center to the center of detectors.
% NR is the number of circles.
% ND is the number of detectors.
% SM is the system matrix.
% Note: The coordinate origin is the object center.
%       Consider only the view that the souce is perpendicular to detectors.
%       The source is a parallel beam which is big enough to cover the object.
% Copyright Neolus and The MathWorks, Inc
% Date£º2018/06


% Change to cooddinate system.
% rmax=nr*dg;
SMU=zeros(nd/2,nr);

for k=1:nd/2
    dk=(k-0.5)*dtc;
    
    m=ceil(dk/dg);
   
    dst0=0;
    x=zeros(1,nr);
    for j=m:nr
        rj=j*dg;
        x(j)=2*sqrt((rj^2-dk^2))-dst0;
        dst0=dst0+x(j);        
    end
    SMU(k,:)=x;
    
end

SM=[flipud(SMU);SMU];
SMU=flipud(SMU);
end