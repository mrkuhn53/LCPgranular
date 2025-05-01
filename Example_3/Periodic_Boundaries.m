function [dxcell_pq] = ...
         Periodic_Boundaries(x, V_pq, xcell, idim, periodic_directions)
%
% with periodic boundaries, this function adds "pads" the N particles with
% additional particles on the top and right sides (2D) of the assembly.
% These "ghost" particles are added so that the Type 3 displacement
% restraints can be effected.  Note that this functions assumes that the
% particles have already been placed into the primary cell
%
%------------ Input --------------------------------------
% N =     number of particles, ghost particles will be added to the output N
% M =     number of contacts.  Will not be changed by the function
% x =     Nx3 array of particle positions
% V_pq =  Mx2 array of the p and q particles for each pq contact
% xcell = 3x3 size of the periodic cell
% idim = scalar, 2D or 3D
% periodic_directions = 1x3 vector of periodic directions
%
%------------ Output --------------------------------------
% N =        number of particles, including ghost particles
% x =        Nx3 array of particle positions, including ghost particles
% V_pq =     Mx2 array of the p and q particles for each pq contact, note that
%            the original p or q have been replaced by ghost particles
% p_Padded = nPadx1 list of the particles in the primary cell that are
%            are in contact with the added ghost particles 
% dim_Padded = the dimensions that are straddled
% nPad =     number of ghost particles that have been added, outside of 
%            the periodic cell
%
% number of contacts
  M = size(V_pq,1);
%
% branch vectors between the centers of the two particles.  Note that the
% "x" array of particle positions might have been culled of non-relevant
% particles
  rxcntr = x(V_pq(:,2),:) - x(V_pq(:,1),:);
%
  if    (idim==2 && all(periodic_directions(1:2))) ...
     || (idim==3 && all(periodic_directions(1:3)))
%   the number of periodic cell boundaries between the particle pair
%   (for example, when the pair straddles a periodic boundary)
    rcells = (inv(xcell) * rxcntr')';
%
%   the integer number of periodic boundaries crossed by the particle pair.
    nintx = round(rcells);
  end
%
% for each of the 6 deformation components, assign an integer that describes
% whether changes in the periodic cell dimensions, dxcell(:,:) should be
% added or subracted from differences, du(V_pq(:,2),:) - du(V_pq(:,1),:),
% in finding dudef
  dxcell_pq = zeros(M,6);
  dxcell_pq(:,1) = -nintx(:,1); % dxcell(1,1)  NOTE: the negative signs
  dxcell_pq(:,2) = -nintx(:,2); % dxcell(2,2)
  dxcell_pq(:,3) = -nintx(:,3); % dxcell(3,3)
  dxcell_pq(:,4) = -nintx(:,2); % dxcell(1,2)
  dxcell_pq(:,5) = -nintx(:,3); % dxcell(1,3)
  dxcell_pq(:,6) = -nintx(:,3); % dxcell(2,3)
%
