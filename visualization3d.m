T = readmatrix('./tets.txt')+1;
size(T)
X = readmatrix('./points3d.txt');
size(X) 



% Calculate the centroids of each of the simplices
C = (X(T(:,1),:) + X(T(:,2),:) + X(T(:,3),:) + X(T(:,4),:))/4.0;

xGtHalf = C(:, 1) > 0.5;
yGtHalf = C(:, 2) > 0.5;
zGtHalf = C(:, 3) > 0.5;

acceptMask = ~(~yGtHalf & zGtHalf);%~(xGtHalf & zGtHalf) & ~(yGtHalf;
acceptMask = ~(yGtHalf & zGtHalf);

%C(1,:)




off = 20;
i = 3;
%tetramesh(T(acceptMask,:),X,'FaceAlpha',1)
%camorbit(180,0)

n = 41;
%scatter3(X(:,1), X(:,2), X(:,3));
%plot3(X(:,1),X(:,2),X(:,3),'r+');
%plot3(C(:,1),C(:,2),C(:,3),'r+');
inds = off*(n)*(n):(off+1)*(n)*(n);
%scatter(X(inds,1), X(inds,2))

DT = delaunay(X(inds,1), X(inds,2));
triplot(DT,X(inds,1), X(inds,2));