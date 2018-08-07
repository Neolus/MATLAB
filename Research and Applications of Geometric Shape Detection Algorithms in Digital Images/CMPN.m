function [samples,center,cnt]=CMPN(BW,num,ns)
%CMPN Uses Chord Middle Points method to find NS points randomly which are the same ellipse.
%   [SAMPLES,CENTER,CNT]=CMPN(BW,T) computes the center of an ellipse in binary image BW
%   and returns samples which are the same ellipse. 
%   SIZE(SAMPLES) is the number of ellipses. SAMPLES(K) are three points. 
%   CENTER store all centers of ellipses found by this algorithm.
%   CNT is a counter of found ellipes.
%   NUM is the max number of ellipses which need to find.
%   NS is the number of samples.

[M,N]=size(BW);
idx=find(BW);  % All edge points. A column vector
[X,Y]=idx2xy(idx,M); %Numbers of row and column.

n=length(idx);
done=true; %Loop control
cnt=1;  %Counter of picking points
fn=0; %Fail counter

samples_c=cell(num,1); %Samples of true ellipses.
center=[];  %Centers of ellipses.
while (done&&(cnt<n)&&fn<num)
    % Get ns points with RAND
    idxn=round(1+(n-1)*rand(ns,1));
    [Sx,Sy]=idx2xy(idx(idxn),M); % Get row and column indexes.
    
    % All middle points
    MS=zeros(n,2,ns);
    
    for k=1:ns  % Compute middle points.
        MS(:,:,k)=[X+Sx(k),Y+Sy(k)]/2;
    end
    A=MS(:,:,1); %The crossover points of inscribed ellipses
    for k=1:ns-1  %Find the intersecting pointof all inscribed ellipses.
        A=intersect(A,MS(:,:,k+1),'rows');
    end
    
    if ~isempty(A)  %A is the center of a ellipse.
        center_new=mean(A,1);
        if isempty(center)||(~isfind(center_new(1),center_new(2),center))
            %Whether center_new has been detcected
            fn=fn+1;   %This is a new center.
            samples_c{fn}=[Sx,Sy];  % Add points to SAMPLES.
            center=[center;center_new]; %Add this new center.
        end
    else
        cnt=cnt+1;  %Increase the counter
    end
end
samples=zeros(ns,2,fn); %Place samples.

for k=1:fn  %Cell to matix
    samples(:,:,k)=samples_c{k};
end

for k=1:fn  % Double samples
    sak=samples(:,:,k);
    a=round([2*center(k,1)-sak(:,1),2*center(k,2)-sak(:,2)]);
    %Another endpoints. Center is the middle points.
    done1=true;
    while(done1)   %Waive endpoints which are out of the boundary
        idx_out=find(a(:,1)<0|a(:,1)>M);
        if ~isempty(idx_out)
            a(idx_out,:)=[];
        end
        idy_out=find(a(:,2)<0|a(:,2)>N);
        if ~isempty(idy_out)
            a(idy_out,:)=[];
        end
        done1=length(idx_out)|length(idy_out);
    end
    samples_new(:,:,k)=[samples(:,:,k);a];
end

if fn>=1  % At least one ellipse
    samples=samples_new;
else samples=[];
end
end
% Convert an index of matrix to coordinate x,y. M is the row of this matrix
function [X,Y]=idx2xy(idx,M)
Y=ceil(idx/M);
X=idx-M*(Y-1);
end

function r=isfind(Xc,Yc,XY)
%ISFIND Finds (Xc,Yc) in matrix XY.
% IF (XC,YC) is found, R is true;
%Consider a 3*3 neighbourhood

D=(XY(:,1)-Xc).^2+(XY(:,2)-Yc).^2;
if find(D<100)
    r=true;
else r=false;
end
end

