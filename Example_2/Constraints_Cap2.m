function [ C, dc, dp] = ...
         Constraints_Cap2( ...
           N, Cap, x, Tilt, Shift, Twist, dp, ...
           LoadControl, Cap_Movement, Cap_Load)
%
% constructs the rows of the constraint matrix [C] for the 4 cap spheres,
% in the 4-vector "Cap".  The cap is in the x1-x2 plane.  The spheres
% are arranged counterclockwise when looking downward from above.
%
% Tilt = 0   No tilting of cap allowed
% Tilt = 1   Cap can tilt
%
% initialize
  C = zeros(0, 6*N);
  dc = zeros(0,1);
%
%
% the cap is rigid. Prevent warping of cap
  iRows = 1;
  C(iRows, 3*(Cap([1,3])-1) + 3) =  1;  dc(iRows,1) = 0;
  C(iRows, 3*(Cap([2,4])-1) + 3) = -1;
%
  if ~Tilt
    iRows = iRows + 1;
    C(iRows, 3*(Cap([1,4])-1) + 3) =  1;  dc(iRows,1) = 0;
    C(iRows, 3*(Cap([2,3])-1) + 3) = -1;
%
    iRows = iRows + 1;
    C(iRows, 3*(Cap([1,2])-1) + 3) =  1;  dc(iRows,1) = 0;
    C(iRows, 3*(Cap([3,4])-1) + 3) = -1;
  end
%
  if ~Shift
%   average movements in the x1 and x2 directions
    iRows = iRows + 1;
    C(iRows, 3*(Cap(1:4)-1) + 1) = 1;  dc(iRows,1) = 0;
%
    iRows = iRows + 1;
    C(iRows, 3*(Cap(1:4)-1) + 2) = 1;  dc(iRows,1) = 0;
  end
%
  if ~Twist
%   average twirling of the cap spheres about the x3 axis
    iRows = iRows + 1;
    C(iRows, 3*(Cap-1) + 2) =  x(Cap,1);
    C(iRows, 3*(Cap-1) + 1) = -x(Cap,2);
  end
%
  if LoadControl
%   load control on the cap
    dp(3*(Cap([1,2,3,4])-1) + 3) = (1/4) *  Cap_Load;
  else
%   displacement control of the cap
    iRows = iRows + 1;
    C(iRows, 3*(Cap(1:4)-1) + 3) = 1;  dc(iRows,1) = 4*Cap_Movement;
  end
%
% remove extraneous rows
  C =   C(1:iRows,:);
  dc = dc(1:iRows,:);
