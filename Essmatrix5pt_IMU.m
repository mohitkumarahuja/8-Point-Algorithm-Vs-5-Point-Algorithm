

function [E,e1,e2] = Essmatrix5pt_IMU(P1,P2)
    
%     [x1, x2, npts] = checkargs(varargin(:));
%     Octave = exist('OCTAVE_VERSION') ~= 0;  % Are we running under Octave?    
    
    % Normalise each set of points so that the origin 
    % is at centroid and mean distance from origin is sqrt(2). 
    % normalise2dpts also ensures the scale parameter is 1.
    %[x1, T1] = normalise2dpts(x1);
   % [x2, T2] = normalise2dpts(x2);
    
    % Build the constraint matrix
%     A = [P2(1,1)*P1(1,1) P2(1,1)*P1(1,2) P2(1,1) P2(1,2)*P1(1,1) P2(1,2)*P1(1,2) P2(1,2) P1(1,1) P1(1,2) 1
%          P2(2,1)*P1(2,1) P2(2,1)*P1(2,2) P2(2,1) P2(2,2)*P1(2,1) P2(2,2)*P1(2,2) P2(2,2) P1(2,1) P1(2,2) 1
%          P2(3,1)*P1(3,1) P2(3,1)*P1(3,2) P2(3,1) P2(3,2)*P1(3,1) P2(3,2)*P1(3,2) P2(3,2) P1(3,1) P1(3,2) 1
%          P2(4,1)*P1(4,1) P2(4,1)*P1(4,2) P2(4,1) P2(4,2)*P1(4,1) P2(4,2)*P1(4,2) P2(4,2) P1(4,1) P1(4,2) 1
%          P2(5,1)*P1(5,1) P2(5,1)*P1(5,2) P2(5,1) P2(5,2)*P1(5,1) P2(5,2)*P1(5,2) P2(5,2) P1(5,1) P1(5,2) 1];       

     
     A = [P2(1,1) (P1(1,1)*P2(1,2))-(P2(1,1)*P1(1,2)) (P1(1,2)*P2(1,2))+(P1(1,1)*P2(1,1)) P2(1,2) P1(1,1) P1(1,2)
          P2(2,1) (P1(2,1)*P2(2,2))-(P2(2,1)*P1(2,2)) (P1(2,2)*P2(2,2))+(P1(2,1)*P2(2,1)) P2(2,2) P1(2,1) P1(2,2)
          P2(3,1) (P1(3,1)*P2(3,2))-(P2(3,1)*P1(3,2)) (P1(3,2)*P2(3,2))+(P1(3,1)*P2(3,1)) P2(3,2) P1(3,1) P1(3,2)
          P2(4,1) (P1(4,1)*P2(4,2))-(P2(4,1)*P1(4,2)) (P1(4,2)*P2(4,2))+(P1(4,1)*P2(4,1)) P2(4,2) P1(4,1) P1(4,2)
          P2(5,1) (P1(5,1)*P2(5,2))-(P2(5,1)*P1(5,2)) (P1(5,2)*P2(5,2))+(P1(5,1)*P2(5,1)) P2(5,2) P1(5,1) P1(5,2)
          ];


%     if Octave
% 	[U,D,V] = svd(A);   % Don't seem to be able to use the economy
                        % decomposition under Octave here
%     else
	[U,D,V] = svd(A); % Under MATLAB use the economy decomposition
%     end

    % Extract essential matrix from the column of V corresponding to
    % smallest singular value.
    Er = V(:,6);
    E = [Er(3,1) -Er(2,1) Er(1,1) 
         Er(2,1) Er(3,1) Er(4,1) 
         Er(5,1) Er(6,1) 0 ];
    
    % Enforce constraint that fundamental matrix has rank 2 by performing
    % a svd and then reconstructing with the two largest singular values.
    %[U,D,V] = svd(E,0);
    %E = U*diag([D(1,1) D(2,2) 0])*V';
    
    % Denormalise
%     E = T2'*E*T1;
%     E = 
%     
    if nargout == 3  	% Solve for epipoles
	[U,D,V] = svd(E,0);
	e1 = hnormalise(V(:,3));
	e2 = hnormalise(U(:,3));
    end
    
%--------------------------------------------------------------------------
% Function to check argument values and set defaults

% function [x1, x2, npts] = checkargs(arg);
%     
%     if length(arg) == 2
%         x1 = arg{1};
%         x2 = arg{2};
%         if ~all(size(x1)==size(x2))
%             error('x1 and x2 must have the same size');
%         elseif size(x1,1) ~= 3
%             error('x1 and x2 must be 3xN');
%         end
%         
%     elseif length(arg) == 1
%         if size(arg{1},1) ~= 6
%             error('Single argument x must be 6xN');
%         else
%             x1 = arg{1}(1:3,:);
%             x2 = arg{1}(4:6,:);
%         end
%     else
%         error('Wrong number of arguments supplied');
%     end
%       
%     npts = size(x1,2);
%     if npts < 8
%         error('At least 8 points are needed to compute the fundamental matrix');
%     end
    
