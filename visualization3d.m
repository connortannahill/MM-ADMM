T = readmatrix('./tets.txt')+1;
size(T)
X = readmatrix('./points3d.txt');
size(X) 

% Calculate the centroids of each of the simplices
%Vp = zeros(size(T)(1), 3);

C = (X(T(:,1),:) + X(T(:,2),:) + X(T(:,3),:) + X(T(:,4),:))/4.0;

xGtHalf = C(:, 1) > 0.5;
yGtHalf = C(:, 2) > 0.5;
zGtHalf = C(:, 3) > 0.5;
camorbit(180,0)

acceptMask = ~(~yGtHalf & zGtHalf)%~(xGtHalf & zGtHalf) & ~(yGtHalf;

C(1,:)




off = 5;
i = 3;
tetramesh(T(acceptMask,:),X,'FaceAlpha',1)

%n = 11;
%scatter3(X(:,1), X(:,2), X(:,3));
%plot3(X(:,1),X(:,2),X(:,3),'r+');
%plot3(C(:,1),C(:,2),C(:,3),'r+');
%inds = off*(n)*(n):(off+1)*(n)*(n);
%scatter(X(inds,1), X(inds,2))