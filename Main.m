%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Motion estimation using 8 points algortihm %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                 vs                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    5 points algorithm with IMU           %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name :    Mohit Ahuja
%            MSCV - 2
%      University of Burgundy
%        Le Creusot Campus
%

close all;
clear all;
clc;

%% Test 1 : example with a particular data
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    data generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lowerbound=-300;
upperbound=300;
nbpoints = 50;

P = [lowerbound+(upperbound-lowerbound).*rand(nbpoints,3),ones(nbpoints,1)];

% the data belong on a cube [-300*300]*[-300*300]*[-300*300]

% camera parameter (camera is calibrated)
f=1; u0 = 0; v0 = 0;
K=[f 0 u0;0 f v0;0 0 1];
K1=[f 0 u0 0;0 f v0 0;0 0 1 0]; 

%%%%%%%%%%%%%%%%%%%%%%%%%% Camera 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Camera position at time 1
tx1=10; ty1=-20; tz1=50;
% Rotation angles in degree
roll1=4; pitch1=-5; yaw1=0;

% Rotation of the camera 1       

Rp1=[cos(deg2rad(pitch1)) 0 sin(deg2rad(pitch1));...
     0 1 0;...
     -sin(deg2rad(pitch1)) 0 cos(deg2rad(pitch1))]; %Pitch 
Ry1=[cos(deg2rad(yaw1)) -sin(deg2rad(yaw1)) 0;...
    sin(deg2rad(yaw1)) cos(deg2rad(yaw1)) 0;...
    0 0 1]; % Yaw 
Rr1=[1 0 0;...
    0 cos(deg2rad(roll1)) -sin(deg2rad(roll1));...
    0 sin(deg2rad(roll1)) cos(deg2rad(roll1))]; % Roll
R1=Ry1*Rp1*Rr1;


% Translation of the camera 1  
T1 = [tx1,ty1,tz1]';

% camera image 1 :
for i =1 : nbpoints
P1(i,:) = K1*[R1' , -R1'*T1; 0 0 0 1]*P(i,:)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Camera 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Camera position at time 2
tx2=150; ty2=-60; tz2=100;
% Rotation angles in degree
roll2=3; pitch2=15; yaw2=4;

% Rotation of the camera 2       
Rp2=[cos(deg2rad(pitch2)) 0 sin(deg2rad(pitch2));...
     0 1 0;...
     -sin(deg2rad(pitch2)) 0 cos(deg2rad(pitch2))]; % Pitch 
Ry2=[cos(deg2rad(yaw2)) -sin(deg2rad(yaw2)) 0;...
     sin(deg2rad(yaw2)) cos(deg2rad(yaw2)) 0;...
     0 0 1];  %Yaw
Rr2=[1 0 0;...
     0 cos(deg2rad(roll2)) -sin(deg2rad(roll2));...
     0 sin(deg2rad(roll2)) cos(deg2rad(roll2))]; % Roll
R2 =Ry2*Rp2*Rr2;

% Translation of the camera 2
T2 = [tx2,ty2,tz2]';

% Camera image 2 :
for i =1 : nbpoints
P2(i,:) = K1*[R2' , -R2'*T2; 0 0 0 1]*P(i,:)';
end

%%%%%%%%%%%%%%%%%%%%%%%    display data     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
plot3(P(:,1),P(:,2),P(:,3),'*')
plot3(tx1,ty1,tz1,'g*')
Rx = R1*[100,0,0]';
line([tx1,tx1+Rx(1)],[ty1,ty1+Rx(2)],[tz1,tz1+Rx(3)],'Color','b')
Ry = R1*[0,100,0]';
line([tx1,tx1+Ry(1)],[ty1,ty1+Ry(2)],[tz1,tz1+Ry(3)],'Color','g')
Rz = R1*[0,0,100]';
line([tx1,tx1+Rz(1)],[ty1,ty1+Rz(2)],[tz1,tz1+Rz(3)],'Color','r')
plot3(tx2,ty2,tz2,'r*')
Rx = R2*[100,0,0]';
line([tx2,tx2+Rx(1)],[ty2,ty2+Rx(2)],[tz2,tz2+Rx(3)],'Color','b')
Ry = R2*[0,100,0]';
line([tx2,tx2+Ry(1)],[ty2,ty2+Ry(2)],[tz2,tz2+Ry(3)],'Color','g')
Rz = R2*[0,0,100]';
line([tx2,tx2+Rz(1)],[ty2,ty2+Rz(2)],[tz2,tz2+Rz(3)],'Color','r')
axis equal


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    Preliminary questions %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%% Preliminary Questions: %%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ');

% 1.1-> Relation between Rt, Tt and R1, R2, T1, T2
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('1.1-> Relation between Rt, Tt and R1, R2, T1, T2');
disp('Rt = R1*R2');
disp('Tt = (T1-T2)*R2');
Rt = R2'*R1;
Tt = R2'*(T1-T2);

% 1.2-> Theoretical Essential Matrix E:
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('1.2-> Essential Matrix Using Rt and Tt:');
Et = skew3(Tt)*Rt

% 1.3-> Verification for the First Point
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('1.3-> Verification of Epipolar Geometry for the First Point:');
disp('%%%%%    This should be equal to 0    %%%%%');
Verify = round(P2(1,:)*Et*P1(1,:)')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    8 pts algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%% 8 Point Algorithm: %%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ');

% 2.1.1-> Essential matrix estimation
[F,e1,e2] = fundmatrix(P1',P2');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('2.1.1-> Essential Matrix using 8 Point Algorithm:');
E = K'*F*K

% 2.1.2-> Verification for the first point
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('2.1.2-> Verification of E in 8 Point Algorithm:');
disp('%%%%%    This should be equal to 0    %%%%%');
Verify_For_First_Point = round(P2(1,:)*E*P1(1,:)')

% 2.2.1-> Essential matrix estimation with robust estimation 
[F, inliers] = ransacfitfundmatrix(P1', P2', 0.001);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('2.2.1-> Essential matrix estimation with Robust Estimation:');
 E_After_Ransac = K'*F*K
 
% 2.2.2-> Verification for the first point
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('2.2.2-> Verification after Robust Estimation:');
disp('%%%%%    This should be equal to 0    %%%%%');
Verify_for_First_Point_Ransac_E = round(P2(1,:)*E_After_Ransac*P1(1,:)')

% 2.3-> Essential Matrix decompostion on R and T
[U, R1E, R2E, t1E, t2E] = PoseEMat(E);

% Verification....
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('2.3-> Verification after Decomposition:');
R_1_computed = norm(R1E-Rt);   
R_2_computed = norm(R2E-Rt);  % R2E is the original rotational matrix.
Original_T = Tt/norm(Tt)
T_1_Computed = t1E/norm(t1E);
T_2_Computed = t2E/norm(t2E)  % As the norm is near to the Original_T, So this is our Translation Matrix
disp('We found our translation as T2_Computed = Original T');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   5+1 pts algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%% 5+1 Point Algorithm: %%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ');

% 3.1-> Essential matrix estimation with the vertical known
%       virtual points (derotation of roll and pitch angle)
R1derot = Rp1*Rr1;
R2derot = Rp2*Rr2;
PV1 = (R1derot*P1(:,1:3)')';
PV2 = (R2derot*P2(:,1:3)')';
for i =1 : nbpoints
    PV1(i,:)=PV1(i,:)/PV1(i,3);
    PV2(i,:)=PV2(i,:)/PV2(i,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.2-> Theoretical Essential matrix E :
T = Ry2'*(T1-T2);
R = Ry2'*Ry1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.3.1-> Yaw Angle Computation
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('3.3.1-> Yaw Angle Computation:');
Yaw_Angle = atan2(-R(1,2),R(1,1)) % in radian

% 3.3.2-> Essential Matrix using Yaw Angle
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('3.3.2-> Essential Matrix using Yaw Angle:');
Et_5Point = [-T(3)*sin(Yaw_Angle)  -T(3)*cos(Yaw_Angle)  T(2);...
       T(3)*cos(Yaw_Angle)  -T(3)*sin(Yaw_Angle) -T(1);...
      -T(2)*cos(Yaw_Angle)+T(1)*sin(Yaw_Angle) T(1)*cos(Yaw_Angle)+T(2)*sin(Yaw_Angle) 0]

% 3.3.3-> Verification by epipolar geometry 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('3.3.3-> Verification using 1st Point:');
disp('%%%%%    This should be equal to 0    %%%%%');
Verify_5_POint = round(PV2(1,:)*Et_5Point*PV1(1,:)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.4.1-> Essential matrix estimation using Essmatrix5pt_IMU
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('3.3.2-> Essential Matrix using Essmatrix5pt_IMU:');
[E_New_5Point,e1,e2] = Essmatrix5pt_IMU(PV1,PV2);

% 3.4.2-> Verification by epipolar geometry 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('3.4.2-> Verification using 1st Point:');
disp('%%%%%    This should be equal to 0    %%%%%');
Verify_5_POint_Prof_Func = round(PV2(1,:)*E_New_5Point*PV1(1,:)')

% 3.4.3-> Essential matrix estimation using My Function
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('3.4.3-> My Function to compute E for 5 Point Algorithm');
E_5Point_Algo = Compute_E( PV1,PV2 )

% 3.4.4-> Verification by epipolar geometry 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('3.4.4-> Verification of E using My Function:');
disp('%%%%%    This should be equal to 0    %%%%%');
Verify_5_POint_Mohit = round(PV2(1,:)*E_5Point_Algo*PV1(1,:)')


% 3.4.5-> Essential matrix estimation with robust estimation 
[E_5_point_Robust,inliers] = RansacEssmatrix5pt_IMU(PV1',PV2',0.001);

% 3.4.6-> Verification by epipolar geometry
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('3.4.6-> Verification after Robust Estimation in 5-Point Algo:');
disp('%%%%%    This should be equal to 0    %%%%%');
Verify_5_POint_Robust = round(PV2(1,:)*E_5_point_Robust*PV1(1,:)')


