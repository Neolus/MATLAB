function   polygon(lines,dis)
% POLYGON Dectects polygons of an image.
%    POLYGON(LINES,DIS) finds polygons of deteced lines. LINES is a result
%    of HOUGH_LINES. DIS is the maximum distance of two adjacent points
%    which are on seperate lines. In this case, those two lines are
%    considered to connected.
%    Polygons:
%    Triangle: Equilateral triangle, Isoceles triangle, Right triangle.
%    Quadrilateral: Parallelogram, Rectangle, Square.
%    Pentagon...n-side polygon.
 
if nargin<2
    dis=20;   %Default max distance of two vertex.
end
 
n=size(lines,2);    
end_p=zeros(n,4);
deg_dis=zeros(n,2);
 
if n<3
    error('wrong imput')
else
    
%Initialize
for k=1:n
        end_p(k,:)=[lines(k).point1,lines(k).point2];  %Endpoints
        deg_dis(k,:)=[lines(k).theta,lines(k).length];  % Degrees and length of lines 
end
    
    adj=adjacent(end_p,dis); % Creat an adjacent matrix
    linen=[];
    
    for k=1:n % Every line
        adjk=adj;        
        adjk{2*k-1}(1)=[];
        adjk{2*k}(1)=[];
        path=findPath(2*k-1,2*k,adjk); % Find the path form an endpoint to another.
        path=sort(path); 
        
        i=1;
        while i<length(path)   %Delete some bad vertex.            
            if mod(path(i),2)==0 %Even
                if i==1  %The first point
                    path(i)=[];
                else
                    if path(i)-1==path(i-1) %Successive index is good
                        i=i+2;
                    else  %Un-successive.
                        path(i)=[];
                    end
                end
                       
            else   %Odd
                if path(i)+1==path(i+1)
                    i=i+2;
                else
                    path(i)=[];
                end  
            end
        end
        
        m=length(path);
        % The indexes of lines which are on the same polygon.
        line=round(path(1:2:m)/2); 
        line_num=m/2; %The number of sides of polygon.
        line_ind=0;
 % Index of lines in decimal system. eg.matrix [2 3 5 6]-> number '6532'        
        for j=1:line_num  %Code the line indexes.
            line_ind=line_ind+line(j)*10^(j-1);
        end
        
        tol=2e-2;             %Find if this polygon has been detected
        if (~isempty(line))&(isempty(linen)|(isempty(find(linen==line_ind))))
            linen=[linen;line_ind];  %Not been detected            
            if line_num==3  %Triangle                
                dis3=deg_dis(line,2);
                dis3=sort(dis3);     %Pythagorean           
             if (sqrt(abs(dis3(1)^2+dis3(2)^2))-dis3(3))/dis3(3)<(4*tol) 
                    fprintf('This is a Right triangle.\n')
                elseif (abs(dis3(1)-dis3(2))/dis3(2)<tol)|(abs(dis3(3)-dis3(2))/dis3(2)<tol) % Lsoceles triangle                    
                    if abs(mean(dis3)-dis3(2))/dis3(2)<tol  % Equilateral triangl
                        fprintf('This is an Equilateral triangle.\n')                        
                    else
                        fprintf('This is a Isoceles triangle.\n')
                    end
                else
                    fprintf('This is a Triangle\n');
                end
                
            elseif line_num==4 %Quadrilateral
                deg4=(deg_dis(line,1)); %Degrees of four lines
                dis4=deg_dis(line,2);   %Lengths of four lines 
                
                degm=abs(deg4(2:4)-deg4(1)); %Use the first as a standard.
                degm=sort(degm);
                if degm(1)<=5&&(degm(3)-degm(2))<=8  %Parallelogram
                    
                    if abs(degm(2)-90)<=5&&abs(degm(3)-90)<=5  %rectangle                        
                        if abs((dis4-mean(dis4)))/mean(dis4)<4*tol %Square
                            fprintf('This is a Square.\n');
                        else
                            fprintf('This is a Rectangle.\n');
                        end
                    else
                        fprintf('This is a Parallelogram.\n');
                    end
                else
                    fprintf('This is a general Quadrilateral.\n');
                end
               
            else  % Others  num>=5
                fprintf('This is a %d-sided polygon.\n',line_num);
            end
        end
    end
end
end 

function adj_mat=adjacent(endpoints,dis)
%ADJACENT Computes and connects the adjacent vertex.
%    ADJ_MAT=ADJACENT(ENDPOINTS,DIS) creats an adjacent matrix, ADJ_MAT,
%    which is a graph. DIS is the maximum threshold distance of two points.
 
 n=length(endpoints);
disq=dis^2;
adj_mat=cell(2*n,1); %2*n vertexes of n lines. vertex k and k+1 are endpoints of line k/2;
 
%Initialize ADJ_MAT. Connect endpoints of a line.
for k=1:2:2*n
    adj_mat{k}=k+1;
    adj_mat{k+1}=k;
end
 
%Compute adjacent matrix with  a distance.
for k=1:2*n
    ke=uint8(k/2);  %The index of vertex k in ENDPOINTS.
    
    if mod(k,2)==1
        Pk=endpoints(ke,1:2);% The current point.
    else
        Pk=endpoints(ke,3:4);% The current point.
    end
    
    %Connect two near vertexes.
    for j=ke+1:n  %Traverse
        ak=endpoints(j,:);
        %Two endpoints of one line.
        dis1=((ak(1)-Pk(1))^2+(ak(2)-Pk(2))^2);
        dis2=((ak(3)-Pk(1))^2+(ak(4)-Pk(2))^2);
         
        if dis1<disq  %Connect each other.
            adj_mat{k}=[adj_mat{k},2*j-1];
            adj_mat{2*j-1}=[adj_mat{2*j-1},k];
        end
        
        if dis2<disq
            adj_mat{k}=[adj_mat{k}, 2*j];
            adj_mat{2*j}=[adj_mat{2*j},k];
        end
    end
end
end
 
function path1=findPath(fs,ls,adj_mat)
%FINDPATH Finds a path from form vertex FS to vertex LS.
%    PATH=FINDPATH(FS,LS,ADJ_MAT), the driver function of DSF. ADJ_MAT is
%    the adjacent matrix of vertexes.  IF PATH1 is empty,  there is no such a circulate.

n=size(adj_mat,1);
path=1:n;
visited=false(n,1);

path=dfs(fs,ls,visited,adj_mat,path);  %Call DFS to get path form FS to LS.
if path(ls)==ls  %Not find such a path
    path1=[];
else  
    k=ls;
    res=[];
    while k~=fs
        res=[res,k];
        k=path(k);
    end
    res=[res,fs];
    path1=res(end:-1:1); %Reverse PATH
end 
end
 
function [path,visited]=dfs(fs,ls,visited,adj_mat,path)
%DFS: Depth-first Search method to traverse form vertex FS to vertex LS.
%    PATH=DFS(FS,LS,VISITED,ADJ_MAT,PATH) computes and return a route:
%    PATH contains previous vertex index. VISITED is a vector which stores
%    the access information of vertexes. DFS is a recusive function.
%    information of whether vertexes is vistited.
 
visited(fs)=true;
na=length(adj_mat{fs});
 
if find(adj_mat{fs}==ls) % Find the last vertext
    path(ls)=fs;
else
    for k=1:na
        if ~visited(adj_mat{fs}(k))
            path(adj_mat{fs}(k))=fs;  %Set PATH
            [path,visited]=dfs(adj_mat{fs}(k),ls,visited,adj_mat,path); %DFS
        end
    end
end
end
