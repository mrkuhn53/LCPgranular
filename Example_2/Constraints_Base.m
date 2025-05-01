function [C, dc] = Constraints_Base(N, Base, x, Turn)
%
% constructs the rows of the constraint matrix [C] for the 4 base spheres,
% in the 4-vector "Base".  The base is in the x1-x2 plane.  The spheres in 
% the base are arranged counterclockwise when looking downward from above
%
  if Turn
    C = zeros(6,6*N);
    dc = zeros(6,1);
  else
    C = zeros(7,6*N);
    dc = zeros(7,1);
  end
%
% zero vertical movements
  C(1, 3*(Base(1)-1) + 3) = 1;  dc(1) = 0;
  C(2, 3*(Base(2)-1) + 3) = 1;  dc(2) = 0;
  C(3, 3*(Base(3)-1) + 3) = 1;  dc(3) = 0;
  C(4, 3*(Base(4)-1) + 3) = 1;  dc(4) = 0;
%
% average movements in the x1 and x2 directions
  C(5, 3*(Base-1) + 1) = 1;  dc(5) = 0;
  C(6, 3*(Base-1) + 2) = 1;  dc(6) = 0;
%
  if ~Turn
%   average twirling of the base spheres about the x3 axis
    C(7, 3*(Base-1) + 2) =  x(Base,1);
    C(7, 3*(Base-1) + 1) = -x(Base,2);
  end
