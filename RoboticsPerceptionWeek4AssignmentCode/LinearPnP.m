function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

A=[];
for i=1:length(X)
   
  X_i=X(i,:);
  x_i=x(i,:);
  
  a_x=[X_i(1) X_i(2) X_i(3) 1 0 0 0 0 -x_i(1)*X_i(1) -x_i(1)*X_i(2) -x_i(1)*X_i(3) -x_i(1)];
  a_y=[0 0 0 0 X_i(1) X_i(2) X_i(3) 1 -x_i(2)*X_i(1) -x_i(2)*X_i(2) -x_i(2)*X_i(3) -x_i(2)];
  A=[A;a_x;a_y];
end

[U D V]=svd(A);

P=reshape(V(:,end),4,3)';

temp=inv(K)*P;

R=temp(:,1:3);
t=temp(:,end);

[U D V]=svd(R);

R=U*V';

t=t./D(1,1);

if det(U*V')<0
  
  R=-1.0*R;
  t=-1*t;
  end
C=-1.0*R'*t;

end





