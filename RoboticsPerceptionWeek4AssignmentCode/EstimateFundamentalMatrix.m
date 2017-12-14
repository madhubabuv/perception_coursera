function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2


A=zeros(length(x1),9);


for i=1:length(x1)
   
    x_1=x1(i,:);
    x_2=x2(i,:);
    
   A(i,:)=[x_1(1)*x_2(1) x_1(2)*x_2(1) x_2(1)...
           x_1(1)*x_2(2) x_1(2)*x_2(2) x_2(2)...
           x_1(1) x_1(2) 1];
   
   
 end
% 
 [U D V]=svd(A);
% 
 F=V(:,end);
% 
 F=reshape(F,3,3)';
% 
 [U D V]=svd(F);
% 
 D(end,end)=0.0;
% 
 F=U*D*V';
 
 F=F./norm(F);

end
