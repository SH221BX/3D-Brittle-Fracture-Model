function [V,C]=voronoi3d_cuboid(Seed,Vcub,tol)
if nargin < 3
    tol=10;
end
eps=10^(-tol);

% ========================================================
%Verify if 'Vcub' forms a rectangular cuboid
CUB0=min(Vcub); CUB1=max(Vcub);
for x=[CUB0(1) CUB1(1)];for y=[CUB0(2 ...
            ) CUB1(2)];for z=[CUB0(3) CUB1(3)] %#ok<ALIGN>
if any([x y z]~=0)
    idx=find(Vcub(:,1)==x & Vcub(:,2)==y & Vcub(:,3)==z, 1);
    if isempty(idx)
        error('#1 your input "Vcub" does not form a rectangular cuboid')
    end
end
end;           end;           end

   
% ========================================================
% Add seeds to avoid infinite cells inside the borders 
% during the voronoi diagrams

dim=max([Seed; Vcub])-min([Seed; Vcub]);
mid=(max([Seed; Vcub])+min([Seed; Vcub]))/2; 
LIMITER=zeros(3*9-1,3); n=0;
for i=[-1 0 1];for j=[-1 0 1];for k=[-1 0 1] %#ok<ALIGN>
if any([i j k]~=0)
    n=n+1;
     LIMITER(n,:)=3*(dim.*[i j k])/2+mid;
%     scatter3(LIMITER(:,1),LIMITER(:,2),LIMITER(:,3),55, ...
%          'Marker','o','MarkerFaceColor',[0 0 1], 'MarkerEdgeColor','k');
%     hold on
end
end;           end;           end
% disp(LIMITER);

% dim = max([Seed; Vcub]) - min([Seed; Vcub]);
% mid = (max([Seed; Vcub]) + min([Seed; Vcub])) / 2;
% LIMITER = [];
% for i = [-1,0.5, 1]
%     for j = [-1, 1]
%         for k = [-1, 1]
%             LIMITER = [LIMITER; (dim .* [i j k]) / 2 + mid];
%         end
%     end
% end
% scatter3(LIMITER(:,1), LIMITER(:,2), LIMITER(:,3), 55, 'Marker', 'o', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'k');
% hold on;
bx = 100;
half = bx*2;
BX = [-half -half -half; half -half -half; -half half -half; half half -half; -half -half half; half -half half; -half half half; half half half];

% ========================================================
% Classical voronoi diagram without limits
[V0,c]=voronoin([Seed; BX]);
C=c(1:size(Seed,1));
V0=round(V0,tol);


% ========================================================
% Find the polyhedron that contains the corners and 
% add it to their vertices
for j=1:size(Vcub,1)
    idx=1;
    MIN=norm(Vcub(j,:)-Seed(1,:));
for i=2:size(Seed,1)
    tamp=norm(Vcub(j,:)-Seed(i,:));
    if tamp<MIN
        MIN=tamp;
        idx=i;
    elseif tamp==MIN
        idx=[idx i];
    end    
end
    M=Vcub(j,:);
    idxV=find(V0(:,1)==M(1) & V0(:,2)==M(2) & V0(:,3)==M(3));
    if isempty(idxV)
        V0=[V0; M]; idxV=size(V0,1); 
    elseif length(idxV)>1
        error('#2 There are duplicates in V')
    end
    for k=idx
        C{k}=[C{k} idxV];
    end
end


% ========================================================
% Find all the points outside the box (==1 for inside)
V_tst=zeros(size(V0,1),1);
for pt=1:size(V0,1)
if all(V0(pt,:)>=CUB0) && all(V0(pt,:)<=CUB1)
    V_tst(pt)=1;
end
end


% ========================================================
% ========================================================
% Cut all the polyhedron on the edges of the box
M=[0 0 0];
% List of box faces
      %i0 i1 i2 M(i0)
listF=[1 2 3 CUB0(1); 2 1 3 CUB0(2); 3 2 1 CUB0(3);
       1 2 3 CUB1(1); 2 1 3 CUB1(2); 3 2 1 CUB1(3);];
% List of box lines
      %i0 i0dim i1 i1dim i2
listL=[1 CUB0(1) 2 CUB0(2) 3; 1 CUB0(1) 3 CUB0(3) 2; 3 CUB0(3) 2 CUB0(2) 1;
       1 CUB1(1) 2 CUB0(2) 3; 1 CUB1(1) 3 CUB0(3) 2; 3 CUB1(3) 2 CUB0(2) 1;
       1 CUB0(1) 2 CUB1(2) 3; 1 CUB0(1) 3 CUB1(3) 2; 3 CUB0(3) 2 CUB1(2) 1; 
       1 CUB1(1) 2 CUB1(2) 3; 1 CUB1(1) 3 CUB1(3) 2; 3 CUB1(3) 2 CUB1(2) 1];
C_tst=zeros(length(C),1);
% for k=954
for k=1:length(C)
Ck=C{k};

if all(Ck~=1) && all(V_tst(Ck)==1) % All vertices are already inside the box
    C_tst(k)=1; %disp([1 k])    
    
elseif all(Ck~=1) && ... % The polyhedron is finite but has parts outside the box
       ( (any(V_tst(Ck)==1) && any(V_tst(Ck)==0)) || all(V_tst(Ck)==0) )
   
    C_tst(k)=2; %disp([2 k])    
    Vk = V0(Ck,:);
    C{k}=Ck(V_tst(Ck)==1);
    Fk=convhull(Vk);    

    % ========================================================
    % for each face of the polyhedron look for :
    for f=1:size(Fk,1)
        
        % Intersection between cube's edges and polyhedron's face f
        Ap=Vk(Fk(f,1),:); Bp=Vk(Fk(f,2),:); Cp=Vk(Fk(f,3),:);     
        for l=1:12
            i1=listL(l,1); i2=listL(l,3); i3=listL(l,5);
            Al=zeros(1,3); Al(i1)=listL(l,2); Al(i2)=listL(l,4); Al(i3)=CUB0(i3); 
            Bl=zeros(1,3); Bl(i1)=listL(l,2); Bl(i2)=listL(l,4); Bl(i3)=CUB1(i3); 
            [M,TST]=intersection_line_plan(Ap,Bp,Cp,Al,Bl,eps);
            if TST
                M=round(M,tol);
                idx=find(V0(:,1)==M(1) & V0(:,2)==M(2) & V0(:,3)==M(3));
                if isempty(idx)
                    V0=[V0; M]; idx=size(V0,1); 
                elseif length(idx)>1
                    error('#3 There are duplicates in V')
                end
                C{k}=[C{k} idx];
            end
        end

        % Intersection between cube's face p and poly's edge v
        for v=1:3
            Ai=Fk(f,v); Bi=Fk(f,mod(v+1-1,3)+1);
            A=Vk(Ai,:); B=Vk(Bi,:);
            for p=1:6
            i0=listF(p,1); i1=listF(p,2); i2=listF(p,3); 
            if (B(i0)-A(i0))~=0
                M(i0)=listF(p,4); t=(M(i0)-A(i0))/(B(i0)-A(i0)); 
                M(i1)=A(i1)+t*(B(i1)-A(i1));
                M(i2)=A(i2)+t*(B(i2)-A(i2));
                if (all(M>=(min([A; B])-eps)) && all(M<=(max([A; B])+eps)))
                if (all(M>=(CUB0-eps)) && all(M<=(CUB1+eps)))  
                    M=round(M,tol);              
                    idx=find(V0(:,1)==M(1) & V0(:,2)==M(2) & V0(:,3)==M(3));
                    if isempty(idx)
                        V0=[V0; M]; idx=size(V0,1); 
                    elseif length(idx)>1
                        error('#4 There are duplicates in V')
                    end
                    C{k}=[C{k} idx];                
                end
                end
            end
            end
        end
        
    end
    
    %Suppress all unecessary multiple points placed on a single same edge
    C{k}=unique(C{k});
    Vk=V0(C{k},:);
    LIST=[];
    for i0=1:size(Vk,1)
        l=[1:(i0-1) (i0+1):size(Vk,1)];
        x=(Vk(l,:)-Vk(i0,:))./vecnorm(Vk(l,:)-Vk(i0,:),2,2);
        t=acosd(round(x*x',tol));
        if all(t<180)
            LIST=[LIST i0];
        end
    end
%     disp(C{k})
%     disp(Vk)
%     disp(LIST)
    C{k}=C{k}(LIST);
    
else % All the remaining infinite cells
    C{k}=[]; C_tst(k)=0;
end

end

% % Round the vertices
% Vrnd=round(V0,tol);
Vrnd=V0;

% Eleminate vertices duplicates
iV=sort(unique([C{:}]));
V=sortrows(unique(Vrnd(iV,:),'rows'),[3 2 1]);
mergeTolerance = 1;

% Supress all the unused vertices and rearrange their numbers
for k=1:length(C)    
if ~isempty(C{k})
    C{k}=sort(unique(C{k}));
    Vk = Vrnd(C{k},:);
    
    % If all the vertices are outside, or all on the edges the box, the cell 
    % is none existent.
    if any(all(Vk<=CUB0)) || any(all(Vk>=CUB1))
        C{k}=[];
        C_tst(k)=0;
    else
        Fk = convhull(Vk); 
        %if exist('mergeNearlyCoplanarFaces.m','file')==2
      %  [Vk, ~] = mergeNearlyCoplanarFaces(Vk, Fk, mergeTolerance);
%         else
%             disp('It would be nice to install geom3d tool box')
%             disp('https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d')
       % end
        C{k}=zeros(1,size(Vk,1));
        for i=1:size(Vk,1)
            idx=find(V(:,1)==Vk(i,1) & V(:,2)==Vk(i,2) & V(:,3)==Vk(i,3));
            if length(idx)>1
                error('#5 There is an issue (it shouldn''t be empty)')
            end
            C{k}(i)=idx(1);
        end
    end

else
    C_tst(k)=0;
end
end

end