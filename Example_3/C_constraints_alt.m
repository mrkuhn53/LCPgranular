function [ C_w_dxcell, dc_w_dxcell, C_wo_dxcell, dc_wo_dxcell, ...
           C_wo_constr, dc_wo_constr] = ...
         C_constraints_alt(N_Rel, C, dc, dxcell_vec, iStress, idim);
%
% this function appends to the constraint matrix "C" the set of constraints
% on strain
%
% N_Rel is the number of relevant particles (those in the force-bearing
% network of particles).
  N = N_Rel;
%
  if idim==2
%   the 2D stress components: 11, 22, 12
    K = [1, 2, 4];
  elseif idim==3
    K = [1:6];
  end
%
% the known stress increments.  Do the following to avoid including 
% extraneous stress components
  iStress = intersect(iStress, K);
%
% the known strain increments
  [iStrain, kStrain] = setdiff(K, iStress);
%
% additional rows of the constraint matrix "C" from the known strain
% increments
  C_boundary = zeros(length(iStrain),6*N+6);
  for i = 1:length(iStrain)'
    C_boundary(i,6*N+iStrain(i)) = 1;
  end
%
% the constraint matrix "C" that appends the known strain increments
  C_w_dxcell =  [C; ...
                 C_boundary];
  dc_w_dxcell = [dc; ...
                 dxcell_vec(iStrain)];
%
% a C matrix with no dxcell columns
  C_wo_dxcell =  C(:,1:end-6);
  dc_wo_dxcell = dc;
%
% a C matrix with dxcell columns, but with no constraints on dxcell
  C_wo_constr = C;
  dc_wo_constr = dc;
%
  if idim==2
%   select the rows of the stacked vector [du / dtheta] that apply to 
%   the 2D case
    Columns2D = [ 1:3:3*N, 2:3:3*N, (3*N+3):3:6*N, 6*N+[1,2,4]];
    Columns2D = sort(Columns2D);
%
    C_w_dxcell = C_w_dxcell(:,Columns2D);
    C_wo_dxcell = C_wo_dxcell(:, [ 1:3:3*N, 2:3:3*N, (3*N+3):3:6*N]);
    C_wo_constr = C_wo_constr(:, Columns2D);
  end
%
