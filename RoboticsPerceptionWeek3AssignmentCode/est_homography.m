function [ H ] = est_homography(video_pts, logo_pts)
% est_homography estimates the homography to transform each of the
% video_pts into the logo_pts
% Inputs:
%     video_pts: a 4x2 matrix of corner points in the video
%     logo_pts: a 4x2 matrix of logo points that correspond to video_pts
% Outputs:
%     H: a 3x3 homography matrix such that logo_pts ~ H*video_pts
% Written for the University of Pennsylvania's Robotics:Perception course

% YOUR CODE HERE
H = [];
A=[];

for i=1:length(video_pts)
  
   x_1=video_pts(i,:);
   
   x_2=logo_pts(i,:);
   
   
   a_x=[-x_1(1) -x_1(2) -1 0 0 0 x_1(1)*x_2(1) x_1(2)*x_2(1) x_2(1)];
   
   a_y=[0 0 0 -x_1(1) -x_1(2) -1 x_1(1)*x_2(2) x_1(2)*x_2(2) x_2(2)];
   
   a=[a_x;a_y];
   
   A=[A;a];
    
  end
   
   [U S V]=svd(A);
     
   h=V(:,end);
   
   H=reshape(h,3,3)';
   

end

