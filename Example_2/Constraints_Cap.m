function [C, dc] = Constraints_Cap(N, Cap, Tilt, Shift, Movement)
%
% constructs the rows of the constraint matrix [C] for the 4 cap spheres,
% in the 4-vector "Cap".  The cap is in the x1-x2 plane.  The spheres in 
% the base are arranged counterclockwise when looking downward from above.
%
% Tilt = 0   No tilting of cap allowed
% Tilt = 1   Cap can tilt
%
  if Tilt
    C = zeros(2,6*N);
    dc = zeros(2,1);
%
%   average movement in the x3 direction
    C(1, 3*(Cap(1:4)-1) + 3) = 1;  dc(1) = 4*Movement;
%
%   no warping of the cap
    C(2, 3*(Cap([1,3])-1) + 3) =  1;  dc(2) = 0;
    C(2, 3*(Cap([2,4])-1) + 3) = -1;
  elseif ~Tilt
    C = zeros(4,6*N);
    dc = zeros(4,1);
%
%   all 4 spheres have equal vertical movement
    C(1, 3*(Cap(1)-1) + 3) = 1;  dc(1) = Movement;
    C(2, 3*(Cap(2)-1) + 3) = 1;  dc(2) = Movement;
    C(3, 3*(Cap(3)-1) + 3) = 1;  dc(3) = Movement;
    C(4, 3*(Cap(4)-1) + 3) = 1;  dc(4) = Movement;
  end
%
  if ~Shift
    C_sub = zeros(2,6*N);
    dc_sub = zeros(2,1);
%
% average movements in the x1 and x2 directions
    C_sub(1, 3*(Cap(1:4)-1) + 1) = 1;  dc_sub(1) = 0;
    C_sub(2, 3*(Cap(1:4)-1) + 2) = 1;  dc_sub(2) = 0;
%
    C = [C; C_sub];
    dc = [dc; dc_sub];
  end
