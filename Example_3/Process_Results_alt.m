function [ du, u, dtheta, Qp, df_pq, dm_pq, dxcell, xcell, ...
           dy, dy_moment, dy2, dycell2, dxcell_vec, Path_Stability] = ...
         Process_Results_alt( ...
           N, N_Rel, Map_p_Rel_to_All, u, Qp, ...
           H_x_x, H_x_lambda, H_lambda_x, H_lambda_lambda, ...
           H_xx_BD_inv, C_MP_inv, dp_dstress, dc_w_dxcell, ...
           lambda, S_pq_x_all, S_pq_lambda_all, M_lambda, ...
           map_lambda_to_pq_and_qp, ...
           V_pq, B_pq, FM_model_pq, idim, ...
           nStress, iStress, dxcell, xcell, dxcell_vec, I, J);
%
% this function does the following: 
%   1) computes the displacement/rotation vector from the values of dp, dc,
%      and lambda
%   2) re-shuffles the 6Nx1 displacement/rotation vector to N*3 matrices
%      of du and dtheta
%   3) revises the Qp quaternions of the particles' orientations
%   4) computes the incremental changes in the contact forces.  Although the
%      normal forces could be computed from contact overlaps, the tangential
%      force must be computed from increments.
%
%--------------- INPUT --------------------------
% N =        scalar number of particles.  Note that the input value includes 
%            ghost particles and excludes non-relevant particles.  These must be
%            corrected
% N_Rel =    number of particles inside the primary cells with a minimum
%            number of contacts
% nPad =     number of ghost particles that were created
% Map_p_Rel_to_All = vector mapping the relevant particles to all 
%                    original N particles
% u =        Nx3 array of particle positions
% Qp =       Nx4 array of particle quaternions
% H_x, H_lambda, C_MP_inv, H_x_BD_inv = matrices used to compute the lambdas
% dp =       6Nx1 vector of external force increments
% dc =       vector of contraints
% lambda =   solution of the increments in lambda
% S_pq_x_all, S_pq_lambda_all = arrays used in computing contact forces
% M_lambda = number of active contacts
% V_pq =     Mx2 array of p and q particles for each pq contact
% B_pq =     kinematics matrix for each pq contact
% idim =     2D or 3D
%
%--------------- OUTPUT -------------------------
% du =     Nx3 array of particle displacements
% u =      Nx3 array of update particle positions
% dtheta = Nx3 array of particle rotations
% Qp =     Mx4 ammended particle quaternions
% df_pq =  Mx3 array of contact force increments
% dm_pq =  Mx3 array of contact moment increments
%
% since we have solved lambda, find the particle movements. Later, "dx" is
% the movements of all particles; whereas, "dx4" is the movement of the
% relevant "Rel" particles
  dx4 = (C_MP_inv - H_xx_BD_inv * H_x_x * C_MP_inv) * dc_w_dxcell ...
        + H_xx_BD_inv * dp_dstress ...
        - H_xx_BD_inv * H_x_lambda * lambda;
%
% reaction forces
  dy = - dp_dstress ...
       + H_x_x * dx4 ...
       + H_x_lambda * lambda;
%
  if idim==2
%   particle and cell movements
    dx2 = dx4(1:(end-3));
    dxcell2 = dx4((end-3+1):end);
%
%   reactive forces and  stressses
    dy2 = dy(1:(end-3));
    dycell2 = dy((end-3+1):end);
  elseif idim==3
%   particle and cell movements
    dx2 = dx4(1:(end-6));
    dxcell2 = dx4((end-6+1):end);
%
%   reactive forces and  stressses
    dy2 = dy(1:(end-6));
    dycell2 = dy((end-6+1):end);
  end
%
% insert the solved dxcell components into dxcell
  for k = 1:nStress
    l = iStress(k);
    i = I(l);
    j = J(l);
%
%   insert dxcell
    dxcell(i,j) = dxcell2(l);
%
%   revise this component of xcell
    xcell(i,j) = xcell(i,j) + dxcell(i,j);
%
%   insert the solved dxcell into the vector
    dxcell_vec(l) = dxcell2(l);
  end
%
% Eulearian strain
  ddefme = dxcell * inv(xcell);
%
% compute a parameter to characterize the stability of the current path.
% The paramter is the second-order work in the direction of the current path
  if ~isempty(lambda)
    Path_Stability = dy'*dx4 - dp_dstress'*dx4;
  else
    Path_Stability = 0;
  end
%
% compute the incremental change in the contact forces
  M = size(V_pq,1);
%
% combined vector of incremental forces and moments at contacts
  df_dm_pq = zeros(M,6);
%
% compute the increments of the contact forces, df_dm_pq
  for pq = 1:M
    p = V_pq(pq,1);
    q = V_pq(pq,2);
%
%   columns of 6Nx6N global stiffness matrix that are multiplied by the
%   movements and rotations of p and q.  Remember that all of the displacements
%   are stacked on top of all of the rotations, as in KPD Eqs. 1_1 and 9_1
    if idim==2
      Columns = [          2*(p-1) + [1 2  ], ...  % movement of p
                           2*(q-1) + [1 2  ], ...  % movement of q
                 2*N_Rel + 1*(p-1) + [1    ], ...  % rotation of p
                 2*N_Rel + 1*(q-1) + [1    ]];     % rotation of q
      dx3 = zeros(12+6,1);
      dx3([1,2,4,5,9,12,13:18]) = [dx2(Columns); dxcell_vec];
    elseif idim==3
      Columns = [          3*(p-1) + [1 2 3], ...  % movement of p
                           3*(q-1) + [1 2 3], ...  % movement of q
                 3*N_Rel + 3*(p-1) + [1 2 3], ...  % rotation of p
                 3*N_Rel + 3*(q-1) + [1 2 3]];     % rotation of q
      dx3 = zeros(12+6,1);
      dx3 = [dx2(Columns); dxcell_vec];
    end
%
    df_dm_pq(pq,:) = (  squeeze(S_pq_x_all(pq,:,:)) ...
                      * squeeze(B_pq(pq,:,:)) ...
                      * dx3)';
  end % for pq = 1:M
%
  if all(FM_model_pq == 1)
    for i_lambda = 1:M_lambda
      pq = map_lambda_to_pq_and_qp(i_lambda,1);
      qp = map_lambda_to_pq_and_qp(i_lambda,2);
%
%     the "1" here only applies with simple linear-friction models,
%     FM_model_pq = 1
      df_dm_pq(pq,:) = df_dm_pq(pq,:) ...
                       + (lambda(i_lambda) * S_pq_lambda_all(pq,:,1));
      df_dm_pq(qp,:) = df_dm_pq(qp,:) ...
                       + (lambda(i_lambda) * S_pq_lambda_all(qp,:,1));
    end
  else
     disp(cstrcat('The values of FM_model_pq not accomodated in function', ...
                  ' Revise_Contact_Forces'))
     ERROR_Process_Results
  end
%
% increments of contact force and contact moment
  df_pq = df_dm_pq(:,1:3);
  dm_pq = df_dm_pq(:,4:6);
%
% input dx2 is a 6Nx1 (3D) or 3Nx1 (2D) vector, containing both displacements
% and rotations.  For 3D, the dx2 vector is arranged with the 3Nx1 displacements
% stacked on top of the 3Nx1 rotations. For 2D, the dx2 vector is arranged
% with the 2Nx1 displacements stacked on top of the Nx1 rotations
  if idim==2
    dx_Rel = [dx2(1:2:(2*N_Rel)), dx2(2:2:(2*N_Rel)), zeros(N_Rel,1)];
    dtheta_Rel = [zeros(N_Rel,2), dx2(2*N_Rel+1:end)];
  elseif idim==3
    dx_Rel = [dx2(1:3:(3*N_Rel)), dx2(2:3:(3*N_Rel)), dx2(3:3:(3*N_Rel))];
    dtheta_Rel = [dx2(3*N_Rel+1:3:end), ...
                  dx2(3*N_Rel+2:3:end), ...
                  dx2(3*N_Rel+3:3:end)];
  end
%
% similar for the reactive forces
  if idim==2
    dy_Rel = [dy2(1:2:(2*N_Rel)), dy2(2:2:(2*N_Rel)), zeros(N_Rel,1)];
    dy_moment_Rel = [zeros(N_Rel,2), dy2(2*N_Rel+1:end)];
  elseif idim==3
    dy_Rel = [dy2(1:3:(3*N_Rel)), dy2(2:3:(3*N_Rel)), dy2(3:3:(3*N_Rel))];
    dy_moment_Rel = [dy2(3*N_Rel+1:3:end), ...
                     dy2(3*N_Rel+2:3:end), ...
                     dy2(3*N_Rel+3:3:end)];
  end
%
% map the relevant displacements to the displacements of all particles
  du = zeros(N,3);
  dtheta = zeros(N,3);
%
% the displacements/rotations of relevant particles
  du(Map_p_Rel_to_All,:) = dx_Rel;
  dtheta(Map_p_Rel_to_All,:) = dtheta_Rel;
%
% the reactive forces/moments of relevant particles
  dy = zeros(N,3);
  dy_moment = zeros(N,3);
%
  dy(Map_p_Rel_to_All,:) = dy_Rel;
  dy_moment(Map_p_Rel_to_All,:) = dy_moment_Rel;
%
% displacements of rattler particles
  if 1
%   the set of non-relevant particles
    NonRel = setdiff(1:N, Map_p_Rel_to_All);
    du(NonRel,:) = (ddefme * u(NonRel,:)')';
%
    dtheta(NonRel,3) = -0.5*ddefme(1,2);
    dtheta(NonRel,2) =  0.5*ddefme(1,3);
    dtheta(NonRel,1) = -0.5*ddefme(2,3);
  end
%
% revise the particles' positions
  u = u + du;
%
% now, revise the quaternions, which give the particles' orientations
  [Qp] = Quat_Inc_Rotate(Qp, dtheta);
