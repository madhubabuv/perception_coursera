function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 

X=[];
x1=[x1 ones(length(x1),1)];
x2=[x2 ones(length(x2),1)];
x3=[x3 ones(length(x2),1)];

for i=1:length(x1)

 x= Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), X0(i,:));  
 
 X=[X;x];
 
end
end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
 
 i_max=1;
 i=1;
while i<=i_max
  i=i+1;
  
    
  J1=Jacobian_Triangulation(C1,R1,K,X0);
 
  J2= Jacobian_Triangulation(C2,R2,K,X0); 
  
  J3= Jacobian_Triangulation(C3,R3,K,X0);
  
  b=[x1(1) x1(2) x2(1) x2(2) x3(1) x3(2)]';
  
  J=[J1 ;J2; J3];
  
  u1=K*R1*(X0'-C1);
  
  u2=K*R2*(X0'-C2);
  
  u3=K*R3*(X0'-C3);
  
  f_x=[u1(1)/u1(end) u1(2)/u1(end) u2(1)/u2(end) u2(2)/u2(end) u3(1)/u3(end) u3(2)/u3(end)]';
  
  del_x=inv(J' * J) * J' * (b-f_x);
 
  X0=X0+del_x';
 end
 X=X0;
end

function J = Jacobian_Triangulation(C, R, K, X)
  px=K(1,3);
  py=K(2,3);
  f=K(1,1); 

  u_x=[f*R(1,1)+px*R(3,1) f*R(1,2)+px*R(3,2) f*R(1,3)+px*R(3,3)];
  
 
  v_x=[f*R(2,1)+py*R(3,1) f*R(2,2)+py*R(3,2) f*R(2,3)+py*R(3,3)];
 
  w_x=[R(3,1) R(3,2) R(3,3)];
    
  u=K*R*(X'-C);
    
  a=(u(end).*u_x-u(1).*w_x)./u(end)^2;
  
  b=(u(end).*v_x-u(2).*w_x)./u(end)^2;
  
  J=[a;b];
     
end
