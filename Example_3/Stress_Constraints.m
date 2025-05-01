function [ H_stress_x, H_stress_lambda, dstress_vec] = ...
         Stress_Constraints( ...
           N, icont, xcello, dxcell, nStress, iStress, ...
           Vol, stress, stress_target, V_pq, r_pq_p, n_pq, f_pq, ...
           K_pq, K_inv, B_pq, H_m_pq_x, ...
           M_lambda, H_pq_lambda, map_pq_to_lambdas, idim, I, J)
%
% this function does the following computes the rows that will later 
% be appended to H_x, H_lambda, and dp to attain the stress increment dstress
%
% initialize the appended rows
  H_stress_x = zeros(6,6*N+6);
  H_stress_lambda = zeros(6,M_lambda);
%
% number of contacts pq and qp
  M = size(f_pq,1);
%
% if stresses are being controlled...
  if nStress > 0
%
%   compute the increments of the controlled stresses
    dstress3 = zeros(3,3);
    for k = 1:nStress
      i = I(iStress(k));
      j = J(iStress(k));
      dstress3(i,j) = stress_target(i,j) - stress(i,j);
    end
%
%   the stress increments, arranged as a column vector. Note: we will actually
%   use the boundary force = stress*xcell, since the force is conjugate
%   with dxcell
    if idim==2
      dstress_vec = [dstress3(1,1)*xcello(2,2);
                     dstress3(2,2)*xcello(1,1);
                     dstress3(3,3)*xcello(1,1);
                     dstress3(1,2)*xcello(1,1);
                     dstress3(1,3)*xcello(1,1);
                     dstress3(2,3)*xcello(1,1)];
    elseif idim==3
      dstress_vec = [dstress3(1,1)*xcello(2,2)*xcello(3,3);
                     dstress3(2,2)*xcello(1,1)*xcello(3,3);
                     dstress3(3,3)*xcello(1,1)*xcello(2,2);
                     dstress3(1,2)*xcello(1,1)*xcello(3,3);
                     dstress3(1,3)*xcello(1,1)*xcello(2,2);
                     dstress3(2,3)*xcello(2,2)*xcello(1,1)];
    end
%
%   the vector is negated
    dstress_vec = -dstress_vec;
%
%   stress change due to changes in contact vectors dr and df_x (the change
%   for force due to du_def, not dlambda)
    dstress_dr = zeros(3,3,6*N+6);          % add 6 columns for dxcell
    dstress_df_x = zeros(3,3,6*N+6);        % add 6 columns for dxcell
    dstress_df_lambda = zeros(3,3,M_lambda);
%
    for pq = 1:M
%     the p and q particles that touch at contact pq
      p = V_pq(pq,1);
      q = V_pq(pq,2);
%
%     columns of 6Nx6N global stiffness matrix that are multiplied by the
%     movements and rotations of p and q.  Remember that all of the 
%     displacements are stacked on top of all of the rotations, as in 
%     KPD Eqs. 1_1 and 9_1
      Columns = [      3*(p-1) + [1 2 3], ...  % movement of p
                       3*(q-1) + [1 2 3], ...  % movement of q
                 3*N + 3*(p-1) + [1 2 3], ...  % rotation of p
                 3*N + 3*(q-1) + [1 2 3], ...  % rotation of q
                 (6*N+1):(6*N+6)];             % the 6 dxcell's
%
%     KPD Eq. 26
      del_r_n_pq = ...
        [0.5*n_pq(pq,:)'*n_pq(pq,:), zeros(3,3)] *squeeze(B_pq(pq,:,:));
%
%     KPD Eq. 27
      del_r_t_pq = squeeze(K_inv(pq,:,:)) ...
                   * ([zeros(3), [          0, n_pq(pq,3), -n_pq(pq,2);
                                  -n_pq(pq,3),          0,  n_pq(pq,1);
                                   n_pq(pq,2) -n_pq(pq,1),           0]] ...
                      - ...
                      squeeze(K_pq(pq,2,:,:)) ...
                        * ( ...
                           [eye(3), zeros(3,3)] ...
                           - ...
                           [n_pq(pq,:)'*n_pq(pq,:), zeros(3,3)] ...
                          ) ...
                     ) ...
                   * squeeze(B_pq(pq,:,:));
%
%     KPD Eq. 25
      del_r_pq = del_r_n_pq + del_r_t_pq;
%
%     KPD Eq. 11
      d_r_pq = del_r_pq + [zeros(3,6), ...
                           [            0,  r_pq_p(pq,3), -r_pq_p(pq,2); ...
                            -r_pq_p(pq,3),             0,  r_pq_p(pq,1); ...
                             r_pq_p(pq,2), -r_pq_p(pq,1),             0], ...
                           zeros(3,3), ...
                           zeros(3,6) ...  % 6 dxcells, not involved
                          ];
%
      dstress_dr_pq = zeros(3,3,12+6); % add 6 dxcells
      dstress_dr_pq(:,1,:) = f_pq(pq,1) * d_r_pq;
      dstress_dr_pq(:,2,:) = f_pq(pq,2) * d_r_pq;
      dstress_dr_pq(:,3,:) = f_pq(pq,3) * d_r_pq;
%
      dstress_dr(:,:,Columns) = dstress_dr(:,:,Columns) + dstress_dr_pq;
%
%     KPD Eq. 19, first term on right
      del_n_pq = -squeeze(K_pq(pq,1,:,:)) * del_r_pq;
      del_n_pq_cross_n_pq = ...
        [n_pq(pq,3)*del_n_pq(2,:) - n_pq(pq,2)*del_n_pq(3,:);
         n_pq(pq,1)*del_n_pq(3,:) - n_pq(pq,3)*del_n_pq(1,:);
         n_pq(pq,2)*del_n_pq(1,:) - n_pq(pq,1)*del_n_pq(2,:)];
%
      f_pq_cross_del_n_pq_cross_n_pq = ...
        [  f_pq(pq,2) * del_n_pq_cross_n_pq(3,:) ...
         - f_pq(pq,3) * del_n_pq_cross_n_pq(2,:); ...
           f_pq(pq,3) * del_n_pq_cross_n_pq(1,:) ...
         - f_pq(pq,1) * del_n_pq_cross_n_pq(3,:); ...
           f_pq(pq,1) * del_n_pq_cross_n_pq(2,:) ...
         - f_pq(pq,2) * del_n_pq_cross_n_pq(1,:)];
%
%     part of the second term on the right of KDP Eq. 19
      del_theta_dot_n = n_pq(pq,:) * squeeze(B_pq(pq,4:6,:));
%
%     KDP Eq. 11 with n_pq
      f_pq_cross_n_pq = [f_pq(pq,2)*n_pq(pq,3) - f_pq(pq,3)*n_pq(pq,2); ...
                         f_pq(pq,3)*n_pq(pq,1) - f_pq(pq,1)*n_pq(pq,3); ...
                         f_pq(pq,1)*n_pq(pq,2) - f_pq(pq,2)*n_pq(pq,1)];
%
%     KPD Eq. 30 (2nd term in KPD Eq. 13)
      del_f_hat_pq = f_pq_cross_del_n_pq_cross_n_pq ...
                     - 0.5 * f_pq_cross_n_pq * del_theta_dot_n;
%
%     KPD last term in Eq. 13
      dtheta_p_cross_f_pq = [zeros(3,6), ...
                             [          0,  f_pq(pq,3), -f_pq(pq,2); ...
                              -f_pq(pq,3),           0,  f_pq(pq,1); ...
                               f_pq(pq,2)  -f_pq(pq,1),           0], ...
                             zeros(3,3), ...
                             zeros(3,6) ...  % 6 dxcells, not involved
                            ];
%
%     KPD Eq. 13
      df_pq_x = squeeze(H_m_pq_x(pq,1:3,:)) + del_f_hat_pq + dtheta_p_cross_f_pq;
%
      dstress_df_x_pq = zeros(3,3,12+6); % add 6 dxcells
%
      dstress_df_x_pq_1 = zeros(1,3,18);
      dstress_df_x_pq_2 = zeros(1,3,18);
      dstress_df_x_pq_3 = zeros(1,3,18);
%
      dstress_df_x_pq_1(1,:,:) = r_pq_p(pq,1) * df_pq_x;
      dstress_df_x_pq_2(1,:,:) = r_pq_p(pq,2) * df_pq_x;
      dstress_df_x_pq_3(1,:,:) = r_pq_p(pq,3) * df_pq_x;
%
      dstress_df_x_pq(1,:,:) = dstress_df_x_pq_1;
      dstress_df_x_pq(2,:,:) = dstress_df_x_pq_2;
      dstress_df_x_pq(3,:,:) = dstress_df_x_pq_3;
%
%     dstress_df_x_pq(1,:,:) = r_pq_p(pq,1) * df_pq_x;
%     dstress_df_x_pq(2,:,:) = r_pq_p(pq,2) * df_pq_x;
%     dstress_df_x_pq(3,:,:) = r_pq_p(pq,3) * df_pq_x;
%
      dstress_df_x(:,:,Columns) = dstress_df_x(:,:,Columns) + dstress_df_x_pq;
%
%     find the effect of lambda on df_pq and on dStress
      for k = find(map_pq_to_lambdas(pq,:))
        df_pq_lambda = squeeze(H_pq_lambda(pq,1:3,k));
%
        dstress_df_lambda_pq = zeros(3,3,1);
        dstress_df_lambda_pq(:,:,1) = r_pq_p(pq,:)' * df_pq_lambda;
%
        dstress_df_lambda(:,:,map_pq_to_lambdas(pq,k)) = ...
          dstress_df_lambda(:,:,map_pq_to_lambdas(pq,k)) ...
          + dstress_df_lambda_pq(:,:,1);
      end % k = find(map_pq_to_lambdas(pq,:))
    end % for pq = 1:M
%
%   divide by the volume
    dstress_dr = (1/Vol) * dstress_dr;
    dstress_df_x = (1/Vol) * dstress_df_x;
    dstress_df_lambda = (1/Vol) * dstress_df_lambda;
%
%   now, we build 5 matrices that are appended to the matrix equation
%      H_x             H_lambda
%      H_stress_x      H_stress_lambda
%
%   times
%      [dx/dxcell lambda]'
%   equals
%      dp
%      dstress_vec = dStress
%
    dstress_dV = zeros(6,6*N+6);
    dstress_dV(1,6*N+1:6*N+6) = ...
      [-(1/xcello(1,1))*stress(1,1), ...
       -(1/xcello(2,2))*stress(1,1), ...
       -(1/xcello(3,3))*stress(1,1), 0, 0, 0];
    dstress_dV(2,6*N+1:6*N+6) = ...
      [-(1/xcello(1,1))*stress(2,2), ...
       -(1/xcello(2,2))*stress(2,2), ...
       -(1/xcello(3,3))*stress(2,2), 0, 0, 0];
    dstress_dV(3,6*N+1:6*N+6) = ...
      [-(1/xcello(1,1))*stress(3,3), ...
       -(1/xcello(2,2))*stress(3,3), ...
       -(1/xcello(3,3))*stress(3,3), 0, 0, 0];
    dstress_dV(4,6*N+1:6*N+6) = ...
      [-(1/xcello(1,1))*stress(1,2), ...
       -(1/xcello(2,2))*stress(1,2), ...
       -(1/xcello(3,3))*stress(1,2), 0, 0, 0];
    dstress_dV(5,6*N+1:6*N+6) = ...
      [-(1/xcello(1,1))*stress(1,3), ...
       -(1/xcello(2,2))*stress(1,3), ...
       -(1/xcello(3,3))*stress(1,3), 0, 0, 0];
    dstress_dV(6,6*N+1:6*N+6) = ...
      [-(1/xcello(1,1))*stress(2,3), ...
       -(1/xcello(2,2))*stress(2,3), ...
       -(1/xcello(3,3))*stress(2,3), 0, 0, 0];
%
    H_stress_x = zeros(6,6*N+6);
    H_stress_lambda = zeros(6,M_lambda);
%
%   add the three sources of stress-stiffness
    for k = 1:6
      i = I(k);
      j = J(k);
%
      H_stress_x(k,:) =   squeeze(dstress_dr(i,j,:))' ...
                        + squeeze(dstress_df_x(i,j,:))' ...
                        + dstress_dV(k,:);
%
      H_stress_lambda(k,:) = squeeze(dstress_df_lambda(i,j,:));
%
%     note: we will use the boundary force instead of the boundary stress,
%     since the force is conjugate with dxcell
%     if idim==2
%       H_stress_x(k,:) = H_stress_x(k,:) * xcello(i,i);
%       H_stress_lambda(k,:) = H_stress_lambda(k,:) * xcello(i,i);
%     elseif idim==3
%       H_stress_x(k,:) = H_stress_x(k,:) * xcello(i,i)*xcello(j,j);
%       H_stress_lambda(k,:) = H_stress_lambda(k,:) * xcello(i,i)*xcello(j,j);
%     end
    end % for k = 1:6
%
    if idim==2
      H_stress_x(1,:) = (xcello(2,2)) * H_stress_x(1,:);
      H_stress_x(2,:) = (xcello(1,1)) * H_stress_x(2,:);
      H_stress_x(3,:) = (xcello(1,1)) * H_stress_x(3,:);
      H_stress_x(4,:) = (xcello(1,1)) * H_stress_x(4,:);
      H_stress_x(5,:) = (xcello(1,1)) * H_stress_x(5,:);
      H_stress_x(6,:) = (xcello(1,1)) * H_stress_x(6,:);
%
      H_stress_lambda(1,:) = (xcello(2,2)) * H_stress_lambda(1,:);
      H_stress_lambda(2,:) = (xcello(1,1)) * H_stress_lambda(2,:);
      H_stress_lambda(3,:) = (xcello(1,1)) * H_stress_lambda(3,:);
      H_stress_lambda(4,:) = (xcello(1,1)) * H_stress_lambda(4,:);
      H_stress_lambda(5,:) = (xcello(1,1)) * H_stress_lambda(5,:);
      H_stress_lambda(6,:) = (xcello(1,1)) * H_stress_lambda(6,:);
    elseif idim==3
      H_stress_x(1,:) = (xcello(2,2)*xcello(3,3)) * H_stress_x(1,:);
      H_stress_x(2,:) = (xcello(1,1)*xcello(3,3)) * H_stress_x(2,:);
      H_stress_x(3,:) = (xcello(1,1)*xcello(2,2)) * H_stress_x(3,:);
      H_stress_x(4,:) = (xcello(1,1)*xcello(3,3)) * H_stress_x(4,:);
      H_stress_x(5,:) = (xcello(1,1)*xcello(2,2)) * H_stress_x(5,:);
      H_stress_x(6,:) = (xcello(2,2)*xcello(1,1)) * H_stress_x(6,:);
%
      H_stress_lambda(1,:) = (xcello(2,2)*xcello(3,3)) * H_stress_lambda(1,:);
      H_stress_lambda(2,:) = (xcello(1,1)*xcello(3,3)) * H_stress_lambda(2,:);
      H_stress_lambda(3,:) = (xcello(1,1)*xcello(2,2)) * H_stress_lambda(3,:);
      H_stress_lambda(4,:) = (xcello(1,1)*xcello(3,3)) * H_stress_lambda(4,:);
      H_stress_lambda(5,:) = (xcello(1,1)*xcello(2,2)) * H_stress_lambda(5,:);
      H_stress_lambda(6,:) = (xcello(2,2)*xcello(1,1)) * H_stress_lambda(6,:);
    end
  else % if ~(nStress > 0)
    dstress_vec = zeros(6,1);
    H_stress_x = zeros(6,6*N+6);
    H_stress_lambda = zeros(6,M_lambda);
  end % if nStress > 0
