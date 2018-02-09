function E = Compute_E( P1,P2 )
% This function is made by Mohit Kumar Ahuja under the guidancce of Prof.
% Cedric to complete the course Multi-Sensor Fusion and Tracking.
%
% This function can be used to compute the Essential Matrix using the 2D
% points of both the Images
% For Example: Essential_Matrix = Compute_E( Point_1, Point_2 );


A = [P2(1,1) (P1(1,1)*P2(1,2))-(P2(1,1)*P1(1,2)) (P1(1,2)*P2(1,2))+(P1(1,1)*P2(1,1)) P2(1,2) P1(1,1) P1(1,2)
     P2(2,1) (P1(2,1)*P2(2,2))-(P2(2,1)*P1(2,2)) (P1(2,2)*P2(2,2))+(P1(2,1)*P2(2,1)) P2(2,2) P1(2,1) P1(2,2)
     P2(3,1) (P1(3,1)*P2(3,2))-(P2(3,1)*P1(3,2)) (P1(3,2)*P2(3,2))+(P1(3,1)*P2(3,1)) P2(3,2) P1(3,1) P1(3,2)
     P2(4,1) (P1(4,1)*P2(4,2))-(P2(4,1)*P1(4,2)) (P1(4,2)*P2(4,2))+(P1(4,1)*P2(4,1)) P2(4,2) P1(4,1) P1(4,2)
     P2(5,1) (P1(5,1)*P2(5,2))-(P2(5,1)*P1(5,2)) (P1(5,2)*P2(5,2))+(P1(5,1)*P2(5,1)) P2(5,2) P1(5,1) P1(5,2)
    ];

% SVD of A
[~,~,V] = svd(A);

% Vectorising V
E = V(:,end);

% Essential Matrix
E = [E(3,1) -E(2,1) E(1,1) 
    E(2,1) E(3,1) E(4,1) 
    E(5,1) E(6,1) 0 ];

end

