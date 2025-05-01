%
% The Cross product of the vector rows in two n X 3 matrices
%
%
function z = Cross(x, y)
%
  if isempty(x) || isempty(y)
    z = [0 0 0];
  else
    if size(x) == size(y)
      if size(x,2) == 3 && size(y,2) == 3
        z = [x(:,2).*y(:,3) - x(:,3).*y(:,2), ...
             x(:,3).*y(:,1) - x(:,1).*y(:,3), ...
             x(:,1).*y(:,2) - x(:,2).*y(:,1)];
      else
        disp('\n  Error 2 in function Cross!\n')
      end
    else
      disp('\n  Error 1 in function Cross!\n')
    end
  end
