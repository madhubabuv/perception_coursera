function [proj_points, t, R] = ar_cube(H,render_points,K)
%% ar_cube
% Estimate your position and orientation with respect to a set of 4 points on the ground
% Inputs:
%    H - the computed homography from the corners in the image
%    render_points - size (N x 3) matrix of world points to project
%    K - size (3 x 3) calibration matrix for the camera
% Outputs: 
%    proj_points - size (N x 2) matrix of the projected points in pixel
%      coordinates
%    t - size (3 x 1) vector of the translation of the transformation
%    R - size (3 x 3) matrix of the rotation of the transformation
% Written by Stephen Phillips for the Coursera Robotics:Perception course

% YOUR CODE HERE: Extract the pose from the homography
    z=H(end,end);
    if z<0
        H=-1.0.*H;
    end
    R=[H(:,1:2) cross(H(:,1),H(:,2))];
    
    [U,D,V]=svd(R);
    D=eye(length(U));
    D(end,end) = det(U*V);
    R=U*D*V;
    t=H(:,3)./norm(H(:,1));
    
% YOUR CODE HERE: Project the points using the pose
    [m,n]=size(render_points);
%     render_points=[render_points ones(m,1)];
%     [m,n]=size(render_points);
%     n
%     homo_R=[R' zeros(n-1,1)]';
%     homo_t=[t;1];
%     new_H=[homo_R homo_t];
    proj_points=[];
    for i=1:m
        
        w=render_points(i,:);
        c=R*w'+t;
        i=K*c;
        i=i./i(end);
        proj_points=[proj_points i];
       
    end
     
    proj_points=proj_points';
    proj_points=proj_points(:,1:2);
   
end
