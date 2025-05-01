function [C, dc] = Constraints_NoRotations(N, Rotations, Cells)
%
% constructs the rows of the constraint matrix [C] to prevent rotations
% of each and every particle
%
  if Rotations==0
%   all particles are not allowed to rotate
    C = zeros(3*N,6*N);
    dc = zeros(3*N,1);
%
    for i = 1:3*N
      C(i,3*N+i) = 1;
    end
%
  elseif Rotations==1
%   all particles are allowed to rotate
    C = zeros(0,6*N);
    dc = zeros(0,1);
%
  elseif Rotations==2 || Rotations==3
%   the exterior spheres are not allowed to rotate, but the interior spheres
%   are allowed to rotate.
%
    if Rotations==2
      if Cells==1
        p_list = [1:4, 6:9];
      elseif Cells==2
        p_list = [1:4, 6:9, 11:14];
      elseif Cells==3
        p_list = [1:4, 6:9, 11:14, 16:19];
      end
    elseif Rotations==3
      if Cells==1
        p_list = [1:4, 6:9];
      elseif Cells==2
        p_list = [1:4, 11:14];
      elseif Cells==3
        p_list = [1:4, 16:19];
      end
    end
%
    nConstraints = 3 * length(p_list);
    C = zeros(nConstraints, 6*N);
    dc = zeros(nConstraints, 1);
%
    ip = 0;
    for p = p_list
      ip = ip + 1;
      C(ip,3*N + 3*(p-1) + 1) = 1;
      ip = ip + 1;
      C(ip,3*N + 3*(p-1) + 2) = 1;
      ip = ip + 1;
      C(ip,3*N + 3*(p-1) + 3) = 1;
    end
  end
