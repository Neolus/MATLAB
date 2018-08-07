function [samples,center,H]=CMHTN(BW,num,ns)
%CMHTN Combines Chord Middle Method and Huogh Transform to get samples.
%    [SAMPLES,CENTER,H]=CMHT(BW,NUM,NS) computes the ceter of an ellipse in binary image BW
%    and returns samples which are the same ellipse. SIZE(SAMPLES) is the
%    number of ellipses. SAMPLES(K) are three points.
%    CENTER store all centers of ellipses found by this algorithm.
%    H is the Hough transformed matrix.
%    NUM is the max number of ellipses which need to find.
%    NS is the number of samples.

[M,N]=size(BW);
H=zeros(M,N,'uint16'); %Hough transformed matrix, as a counter of centers of ellipses
idx=find(BW);  % All edge points. A column vector
[X,Y]=idx2xy(idx,M); %Numbers of row and column.
n=length(idx); %The number of edge point.
nc=uint16(n/2); %The number of points which are used to count center
idxn=uint32(1+(n-1)*rand(nc,1)); % NC samples randomly
[Sx,Sy]=idx2xy(idx(idxn),M);

%Similar to Hough transform
for k=1:nc
    MXY=round([X+Sx(k),Y+Sy(k)]/2);
    idxm=sub2ind(size(H),MXY(:,1),MXY(:,2));
    H(idxm)=H(idxm)+1;
end
%Find peaks of Hough transformed space:H. The locations of peaks are centers of ellipses.
nhood=size(BW)/25;
nhood=max(2*round(nhood/2)+1,1);
[r,c]=hough_peaks(H,num,max(H(:)/6),nhood); %Find peaks of H
center=[r',c'];
Nd=length(r);  %The number of centers.
SE=zeros(nc,2,Nd); %SE is the another endpoint. The middle points of SE and(Sx,Sy) is CENTER.
for  k=1:Nd
    SE(:,:,k)=[2*r(k)-Sx,2*c(k)-Sy]; %Compute SE
end
samples=zeros(2*ns,2,Nd); % Container of samples
cnsk=zeros(Nd,1); %Counter for samples of each ellipse.
for k=1:nc
    d=zeros(num,1);
    for j=1:Nd
        d(j)=isEllipse(SE(k,1,j),SE(k,2,j),7,M,N); %Whether the point SE(k,:,j) is on the ellipse with the center(j)
    end
    d=uint8(d);
    id_first=find(d,1,'first');  %The first ellipse
    id_last=find(d,1,'last');  % The last ellipse
    if ~isempty(id_first)&id_first==id_last
        %(SE(x,y,id_fist)) is a unique sample of an ellipse.
        if cnsk(id_first)<2*ns;  %The numbner of samples of an ellipse.
            %Push back two sampes            samples(cnsk(id_first)+1:cnsk(id_first)+2,:,id_first)=[Sx(k),Sy(k);SE(k,1,id_first),SE(k,2,id_first)];
            %Increase counter of samples(:,:,id_first).
            cnsk(id_first)=cnsk(id_first)+2;
        end
    end
    if min(cnsk(:))>=2*ns
        break;  %Samlpes of all ellipses are enough
    end
end
for k=1:size(samples,3)
    if(find(samples(:,:,k)==0))  %If samples are not enough, discard them
        samples(:,:,k)=[];
        break;
    end
end
% A nested function To see if the point(x,y) is on an ellipse. [M,N]=SIZE(BW);
    function d=isEllipse(x,y,n,M,N)
        if x<1||x>M||y<1||y>N  %Not in the image.
            d=false;
        else
            n1=floor(n/2); %n*n neighbour-hood of (x,y)
            if x-n1<1
                X1=1:n;
            elseif x+n1>M
                X1=M-n1+1:M;
            else
                X1=x-n1:x+n1;
            end
            if y-n1<1
                Y1=1:n1;
            elseif y+n1>N
                Y1=N-n1+1:N;
            else
                Y1=y-n1:y+n1;
            end
            BC=BW(X1,Y1); %Neighbour-hood
            % Once there is 1 in BC, we consider it's on a ellipse.
            if max(BC(:))
                d=true;
            else
                d=false;
            end
        end
    end
end
