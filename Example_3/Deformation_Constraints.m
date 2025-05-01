function [ C, dc] = Deformation_Constraints(xcell, N, x, idim)
%
% this function creates the displacement constraints for assemblies
% with periodic boundaries, by building the constraint matrix C and the
% constraint vector dc
%
%--------------- INPUT -----------------------------
% icont = 1x6 vector of 0's and 1's.  The six directions correspond to 
%         deformation/stress components (1,1), (2,2), (3,3), (1,2), (1,3),
%         and (2,3). A "0" means deformation (strain) control in that
%         component.  A "1" means stress control
% def = the previous (accumulated) deformation gradient relative to the
%       original assembly
% ddefm = the increment in the deformation gradient since the previous
%         time step
% xcell = 3x3 matrix of previous periodic cell dimensions
% dxcell_unknown = 3x3 matrix of controlled stress.  For example, if
%         sigma(1,3) and sigma(2,2) are being controlled (items 3 and 4
%         in the stress vector), Then
%             dxcell_unknown = [0 0 1; 0 2 0; 0 0 0]  (i.e., first and second))
%             nStress = 2  (i.e., two stresses are being controlled)
%             iStress = [3 4]
% N =     scalar number of particles.  This number is reduced by the 
%         non-relevant particles, and the number is increased by the
%         ghost particles that were added so that periodic boundary
%         constraints can be implemented
% x =     Nx3 array that gives the locations of the x particles
% V_pq =  Mx2 incidence matrix for the M contacts.  Row pq contains
%         the p and q particles, as V_pq(pq,p) and V_pq(pq,q)
% p_Padded = list of particles in the primary periodic cell that have
%            a ghost particle in a neighboring cell
% dim_Padded = length(p_Padded)x3 array gives the offset directions 
%              of the ghost particles relative to the p_Padded particle
%              in the primary cell
% nPad =  number of ghost particles that have been added, outside of 
%         the periodic cell
% nStress = number of stresses being controlled
%
%--------------- OUTPUT ---------------------------
% def =      the updated (accumulated) deformation gradient relative to the
%            original assembly
% xcell =    3x3 matrix of updated periodic cell dimensions
% C =        constraint matrix
% dc =       constraint vector
% C_stress = constraints associated with the unkown dxcell 
%            of contorlled stresses
%
% the columns in "C" have the following meaning:
%   The displacement of particle "p" in direction "i" is located in this column
%         3*(p-1) + i
%   The rotation of particle "p" in direction "i" is located in this column
%         3*N + 3*(p-1) + i
%
% the previous cell dimensions
% xcello = xcell;
%
% advance the deformation gradient
% def = def + ddefm;
%
% Eulerian increment of deformation
% ddefme = ddefm*inv(def);
%
% increment in periodic cell dimensions
% dxcell = ddefme*xcell;
%
% advance the periodic cell dimensions
% xcell = xcell + dxcell;
%
% initialize number of displacement constraints
  nConstraints = 0;
%
% we must prevent drifting of the entire assembly
%
  if 1
%   in this approach, we "fix" one particle, to prevent its displacement,
%   so that the entire assembly cannot shift.  We apply this constraint
%   to the particle that is closest to the origin
%
%   total number of constraints: 3 to fix the corner particle, 3*nPad to
%   constrain relative rotations of ghost particles, sum(sum(dim_Padded))
%   to constrain relative displacements of ghost particles
    nTotal = idim;
%
%   initialize constraint matrix and vector
    C = zeros(nTotal,6*N+6); % add 6 columns for dxcell
    dc = zeros(nTotal,1);
%
%   direction of the vector from the origin to the far corner of the
%   periodic cell
    Corner_Direction = sum(xcell,2);
    Corner_Direction = (1/sqrt(sum(Corner_Direction.^2))) * Corner_Direction;
%
%   the particles' projected distances from the origin
    Distances_from_Origin = x * Corner_Direction;
%
%   particle pCorner is the corner particle
    [W, pCorner] = min(Distances_from_Origin);
%
%   fix the corner particle in the three displacement directions
    nConstraints = nConstraints + 1;
    C(nConstraints, 3*(pCorner-1)+1) = 1;
    dc(nConstraints) = 0;
%
    nConstraints = nConstraints + 1;
    C(nConstraints, 3*(pCorner-1)+2) = 1;
    dc(nConstraints) = 0;
%
    if idim==3
      nConstraints = nConstraints + 1;
      C(nConstraints, 3*(pCorner-1)+3) = 1;
      dc(nConstraints) = 0;
    end
  else
%   in this approach, we make sure that the mean translation and
%   rotation (due to particle displacements) of the assembly is 
%   conistent with dxcell
%
%   initialize constraint matrix and vector.  3 tranlations + 3 rotations
    C = zeros(6,6*N+6); % add 6 columns for dxcell
    dc = zeros(6,1);
%
%   mean translation constraint, x_1 direction
    nConstraints = nConstraints + 1;
    C(nConstraints,1:3:3*N) = 1/N;
    C(nConstraints,6*N+1) = -0.5;
    C(nConstraints,6*N+4) = -0.5;
    if idim==3
      C(nConstraints,6*N+5) = -0.5;
    end
    dc(nConstraints) = 0;
%
%   mean translation constraint, x_2 direction
    nConstraints = nConstraints + 1;
    C(nConstraints,2:3:3*N) = 1/N;
    C(nConstraints,6*N+2) = -0.5;
    if idim==3
      C(nConstraints,6*N+6) = -0.5;
    end
    dc(nConstraints) = 0;
%
%   mean translation constraint, x_3 direction
    if idim==3
      nConstraints = nConstraints + 1;
      C(nConstraints,3:3:3*N) = 1/N;
      C(nConstraints,6*N+3) = -0.5;
      dc(nConstraints) = 0;
    end
%
%   mean rotation constraint, about x_3 axis
    nConstraints = nConstraints + 1;
    C(nConstraints,1:3:3*N) = -(1/N)*x(:,2)';
    C(nConstraints,2:3:3*N) =  (1/N)*x(:,1)';
    C(nConstraints,6*N+6) = -(-0.5 / xcell(2,2));
    dc(nConstraints) = 0;
%
    if idim==3
%     mean rotation constraint, about x_2 axis
      nConstraints = nConstraints + 1;
      C(nConstraints,3:3:3*N) = -(1/N)*x(:,1)';
      C(nConstraints,1:3:3*N) =  (1/N)*x(:,3)';
      C(nConstraints,6*N+5) =  (-0.5 / xcell(3,3));
      dc(nConstraints) = 0;
%
%     mean rotation constraint, about x_1 axis
      nConstraints = nConstraints + 1;
      C(nConstraints,2:3:3*N) = -(1/N)*x(:,3)';
      C(nConstraints,3:3:3*N) =  (1/N)*x(:,2)';
      C(nConstraints,6*N+4) = -(-0.5 / xcell(3,3));
      dc(nConstraints) = 0;
    end
  end
