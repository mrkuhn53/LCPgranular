%
% |+|+|+|+|+|+|+|+|+|  Revision 2   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
function [F_pq, M_pq, S_pq_x, S_pq_lambda, Q_pq_x, Q_pq_lambda, ...
          lambdas_pq, k_lambdas_pq] = ...
         LCP3_F_M_Q_pq(FM_model_pq, n_pq, f_pq, m_pq, ...
                      k_pq, alpha_pq, mu_pq, lambda_tolerance, ...
                      other_parameters)
%
% this function finds the stiffness and compatability matrices for a 
% single contact:
%
% In comments below, reference is made to equations in KPD, which is the
% Kuhn, Prunier, Daouadji paper: "Stiffness pathologies in discrete 
% granular systems: Bifurcation, neutral equilibrium, and instability 
% in the presence of kinematic constraints," Int. J. Num. Anal. Methods Geomech.,
% 2019
%
% In KPD Eq. 20-23, the stiffnesses are used for finding the contact force:
%   [df_pq / dm_pq] = [F_pq / M_pq] * [delta_u_pq_def / delta_theta_pq_def]
% where [df_pq / dm_pq] is the stacked vector of contact force and moment
% and [delta_u_pq_def / delta_theta_pq_def] is the stacked vector of
% the relative displacements and rotations at the contact (KPD Eqs. 17-18).
%
% In the linear complementarity formulation of the problem, yield parameters
% are used, with a compatability condition applied to the yield parmaeters:
%     [[[F_pq; M_pq], H_pq_lambda]; 
%      [Q_pq_x, Q_pq_lambda]] ...
%     * [du_pq_def; dlambda] ...
%     = [df_pq; dm_pq; dG_hat_pq]
%
% Note that both F_pq and M_pq are 3x6 matrices.  The stacked matrix 
% [F_pq / M_pq] is 6x6.
%
% Revision 2 is customized for treating the problem as a linear complementarity
% problem.
%
% ---------- INPUT ---------------------------
%   FM_model_pq =  the model for contact force (currently, only a single model,
%               F_model=1, for a linear-frictional contact, as in KPD Section 2.5,
%               with zero contact moments)
%   n_pq =     1x3 (will be converted to 3x1) unit normal vector, outward from p
%   f_pq =     1x3 (will be converted to 3x1) contact force vector, acting on p
%   m_pq =     1x3 (will be converted to 3x1) contact moment vector, acting on p
%   k_pq =     1x1 scalar normal spring stiffness
%   alpha_pq = 1x1 scalar tangential stiffness coefficient
%   mu_pq =    1x1 scalar friction coeficient
%   lambda_tolerance = an kx1 vector that gives the tolerances associated
%              with the ativation of the possible slip conditions (lambdas).
%              For a simple linear-frictional contact, k=1 and the 
%              lambda_tolerance is a scalar
%   other_parameters = currently [], but might be used in future models
%
% ---------- OUTPUT ---------------------------
%   F_pq =     the 3x6 contact force stiffness array of matrices.
%   M_pq =     the 3x6 contact moment stiffness array of matrices.
%              Together, the matrix [F_pq; M_pq] is the matrix S_pq_x
%   S_pq_x =   [F_pq; M_pq]
%   S_pq_lambda = the 6xk matrix that is multiplied by dlambda to adjust
%              the contact force
%   Q_pq_x =   the kx6 vectors associated with the consistency of the lambda
%              parameters of the contact. These vectors are multiplied by the
%              vector u_pq_def.
%   Q_pq_lambda = the kxk scalars associated with the consistency of the lambda
%              parameters of the contact.  These parameters are multiplied by 
%              the lambdas.  The vector Q_pq_lambda is used to build the
%              matrix Q_lambda
%   lambdas_pq = a vector of 0s and/or 1s, indicating whether the lambda
%              is possibly active.  For example, with a linear-frictional
%              contact, lambdas_pq is a 1x1 scalar (only one possible
%              slip mechanism) that is equal to 1 if the tangential force
%              is sufficiently large to possibly slip
%   k_lambdas_pq = the number of lambdas associated with the contact model.
%              For example, is equal to 1 for a simple linear-frictional contact;
%              and is equal to 2 for a linear-frictional contact with rolling
%              friction
%
%   In general, the five matrices are used as follows:
%     [[[F_pq; M_pq], H_pq_lambda]; [Q_pq_x, Q_pq_lambda]] ...
%      * [du_pq_def; dlambda] ...
%     = [df_pq; dm_pq; dG_hat_pq]
%
% ---------- OTHER ----------------------------
%   h_pq =     the 3x1 unit vector of the direction of the tangential force.
%              This vector is used to build the matrices H_pq_lambda and
%              Q_pq_lambda
%
% convert to column vectors
  n_pq = n_pq';
  f_pq = f_pq';
  m_pq = m_pq';
%
  if FM_model_pq==1 
%   a simple linear-frictiona model of contact force.  In KPD, this model 
%   has two stiffness branches: no-slip (branch 1), slip (branch 2)
%   However, in the LCP model, the concept of branches is not used, and
%   the complementarity condition is applied with a single yield parameter,
%   lambda
%
%   initialize the 3 x 6 stiffness matrix for all branches
    F_pq = zeros(3,6);
%
%   the first, elastic (no-slip) branch (KPD Eq. 82_1). The 3x3 matrix of zeros
%   at the end of the matrix nullifies the effect of delta_theta_def_pq 
%   of the relative contact rotation
    F_pq = [k_pq*(alpha_pq*eye(3) + (1-alpha_pq)*n_pq*n_pq'), zeros(3,3)];
%
%   zero contact moment with FM_model_pq==1
    M_pq = zeros(3,6);
%
%   stiffness with respect to [delta_u_def_pq / delta_theta_def_pq]]
    S_pq_x = [F_pq; M_pq];
%
%   the normal force magnitude
    f_dot_n_pq = n_pq' * f_pq;
%
%   tangential contact force.  Needed to find the sliding direction, h_pq,
%   which is needed to find S_pq_lambda and Q_pq_x
    ft_pq = f_pq - f_dot_n_pq * n_pq;
%
%   the sliding direction (KPD Eq. 81)
    ft_pq_mag = sqrt(sum(ft_pq.^2));
%
%   number of lambdas for this type of type of contact mechanism
    k_lambdas_pq = 1;
%
%   stiffness and consistency matrices associated with this contact pq and
%   for all of the lambdas (all k_lambdas_pq of them) that apply to the 
%   contact model of this contact
    S_pq_lambda = zeros(6,k_lambdas_pq);
    Q_pq_x =      zeros(k_lambdas_pq,6);
    Q_pq_lambda = zeros(k_lambdas_pq,k_lambdas_pq);
%
%   test whether lambda is active, with the contact possibly slipping.  Note
%   that f_dot_n_pq is negative when the particles are touching
    lambdas_pq(1) = -f_dot_n_pq > 0 ...
                    & ft_pq_mag > (1-lambda_tolerance(1))*mu_pq*(-f_dot_n_pq);
%
    if lambdas_pq(1)
%     the contact is possibly slipping, so compute the various quantities
%     that are associated with the linear complementarity condition that
%     applies to the slip lambda
%
%     the unit direction of the tangential force
      if ft_pq_mag ~= 0
        h_pq = (1 / ft_pq_mag) * ft_pq;
      else
%       no tangential force at this contact
        h_pq = zeros(3,1);
      end
%
%     used for computing df_pq. The effect of lambda on the contact force
      S_pq_lambda(:,1) = [-alpha_pq*k_pq*h_pq; zeros(3,1)];
%
%     the matrix components that are associated with the consistency conditions
%     for activation of the "lambda" parameters. For FM_model_pq==1, there is
%     only a single vector and a single scalar.  Row vector "-Q_pq_x" is 
%     multiplied by [du_pq_def / dtheta_pq_def]; and the scalar "Q_pq_lambda" 
%     is multiplied by dlambda.  Note that Q_pq_x corresponds to g_pq in
%     the KPD Eq. 80.  Note the negative "-", used for calculating
%     dG_hat = Q_pq_x * dx + Q_pq_lambda * dlambda >= 0
      Q_pq_x(1,:) = -[(k_pq*alpha_pq*h_pq + k_pq*mu_pq*n_pq); zeros(3,1)]';
%
%     Note that if the model includes both sliding and rolling friction
%     (say, with lambdas_pq = 2), then Q_pq_lambda is a k x k matrix.
%     For a simple linear-frictional with a single lambda, then Q_pq_lambda
%     is a 1x1 scalar
      Q_pq_lambda(1,1) = k_pq*alpha_pq;
    end
  else
    disp('Oops! Only FM_model_pq==1 is allowed in LCP_F_M_pq.m')
    ERROR_in_LCP_F_M_pq
  end
