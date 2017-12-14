function X = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)
%% LinearTriangulation
% Find 3D positions of the point correspondences using the relative
% position of one camera from another
% Inputs:
%     C1 - size (3 x 1) translation of the first camera pose
%     R1 - size (3 x 3) rotation of the first camera pose
%     C2 - size (3 x 1) translation of the second camera
%     R2 - size (3 x 3) rotation of the second camera pose
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Outputs: 
%     X - size (N x 3) matrix whos rows represent the 3D triangulated
%       points


x1=[x1 ones(length(x1),1)];

x2=[x2 ones(length(x2),1)];

X=zeros(length(x2),3);

for i=1:length(x1)

  A=[Vec2Skew(x1(i,:))*(K*[R1 -R1*C1]); 
     Vec2Skew(x2(i,:))*(K*[R2 -R2*C2])];
 
  [U,D,V]=svd(A);
    
   V=V(:,end)./V(end,end);

   X(i,:)=V(1:end-1);
  
end

end

