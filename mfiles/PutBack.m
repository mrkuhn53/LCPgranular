function [x, dx] = PutBack(x, dx, xcell, dxcell, periodic_directions)
  for i = 1:3
     if periodic_directions(i)==1
       move = floor(x(:,i)./xcell(i,i));
       x =   x - move *  xcell(:,i)';
       dx = dx - move * dxcell(:,i)';
     end
  end
