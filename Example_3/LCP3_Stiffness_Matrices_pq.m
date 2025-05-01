%
% |+|+|+|+|+|+|+|+|  Revision 2   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
function [H_m_pq_x, H_g1_pq, H_g2_pq, H_g3_pq, B_pq, H_m_max, ...
          H_pq_lambda, Q_pq_x_all, Q_pq_lambda_all, ...
          k_lambdas_pq_all, lambdas_pq_all, S_pq_x_all, S_pq_lambda_all] = ...
  LCP3_Stiffness_Matrices_pq(r_pq_p, r_pq_q, n_pq, f_pq, m_pq, dxcell_pq, ...
          K_pq, K_inv, mu_pq, k_pq, alpha_pq, other_parameters, ...
          lambda_tolerance, FM_model_pq, m_lambdas)
%
% this function computes the stiffnesses (geometric and mechanical) that are
% associated with each contact.  These are local stiffness matrices, 
% not global stiffness matrices. That is, they are of the form
%
%   [db_p / db_q / dw_p / dw_q] = [H_pq] * [du_p / du_q / dtheta_p / dtheta_q]
%
% where [H_pq] is a local stiffness matrix (geometric or mechanical) that
% we mulitply by the movements (and rotations) of particles p and q
% to yield the change in the external force (and moment) on particles p and q.
% In general, [H_pq] has size 12 x 12, since each of the vectors delta_b_p,
% delta_m_p, du_p, etc. is a 3-vector.  
%
% Revision 2 is customized for treating the problem as a linear complementarity
% problem.
%
% In comments below, reference is made to equations in KPD, which is the
% Kuhn, Prunier, Daouadji paper: "Stiffness pathologies in discrete 
% granular systems: Bifurcation, neutral equilibrium, and instability 
% in the presence of kinematic constraints," Int. J. Num. Anal. Methods Geomech.,
% 2019
%
% Note with 2D problems, we will also compute a 12 x 12 matrix, even 
% though many of the elements are zero.  With 2D problems, some rows and
% columns can be ignored, but this neglect is done outside of this function.
%
% Note that the forces and moments on particle q (that is, db_q and dw_q)
% are likely ignored outside of this function, as these forces are computed
% as the qp variant of the contact between p and q.  That is, outside of
% this funciton we will likely ignore six of the twelve rows of the
% H_pq matrices.
%
% This function does not assemble the local stifffnesses to create global
% stiffness and consistency matrices.  This task is done outide of the function
%
%-------------- INPUT ------------------------
% r_pq_p =     2Mx3 matrix of contact vectors.  The 3-vector in row
%              r_pq_p(pq,:) is the vector from the center of particle p
%              to its contact point with particle q.
% r_pq_q =     2Mx3 matrix of contact vectors.  The 3-vector in row
%              r_pq_q(pq,:) is the vector from the center of particle q
%              to its contact point with particle p.
% n_pq =       2Mx3 matrix of unit contact normal vectors.  The 3-vector
%              in row n_pq(ps,:) is the normal vector of contact pq, outward
%              from particle p
% f_pq =       2Mx3 matrix of contact force vectors.  THe 3-vector
%              in row f_pq(pq,:) is the contact force that acts upon particle p
% K_pq =       2Mx2x3x3 array of curvature tensors.  The 3x3 matrix
%              K_pq(pq,1,:,:) is the curvature tensor of particle p at
%              its contact pq with particle q.  The 3x3 matrix K_pq(pq,2,:,:)
%              is the curvature tensor of particle q at its contact pq with
%              particle p.  See KDP Eqs/ 27 and 33.
% K_inv =      2Mx3x3 array of the pseudo-inverse of the curvature sums
%              (K_pq_p + K_pq_q)+ in KDP Eq. 27.  The 3x3 matrix K_inv(pq,:,:)
%              is the pseudo-inverse of the sum for contact pq between 
%              particles p and q
% mu_pq =      2Mx1 vector of contact friction coefficients.  The scalar
%              mu_pq(pq) is the coefficient of contact pq
% k_pq =       2Mx1 vector of contact normal stiffnesses.  The scalar
%              k_pq(pq) is the normal stiffness of contact pq
% alpha_pq =   2Mx1 vector of contact tangential stiffness coefficients.
%              The scalar alpha_pq(pq) is the coefficient of contact pq
% other_parameters = place-holder for other contact properties, as might
%              be needed in the future.  Use [] if not needed.
% FM_model_pq =2Mx1 vector of contact force models.  The model for contact
%              pq is in FM_model_pq(pq).  Only one model is currently
%              programmed: the linear-frictional model of Section 2.5, with
%              zero contact moment
% lambda_tolerance = an kx1 vector that gives the tolerances associated
%              with the ativation of the possible slip conditions (lambdas).
%              For a simple linear-frictional contact, k=1 and the
%              lambda_tolerance is a scalar
%
%------------------- OUTPUT --------------------------------
% H_m_pq_x = 2M x 12 x 12 array of mechanical contact stiffness
%            matrices.
% H_g1_pq =  2Mx12x12 array of the "g-1" geometric contact stiffnesses.
%            The 12x12 matrix H_g1_pq(pq,:,:) is the stiffness of contact pq
% H_g2_pq =  2Mx12x12 array of the "g-2" geometric contact stiffnesses.
%            The 12x12 matrix H_g1_pq(pq,:,:) is the stiffness of contact pq
% H_g3_pq =  2Mx12x12 array of the "g-3" geometric contact stiffnesses.
%            The 12x12 matrix H_g1_pq(pq,:,:) is the stiffness of contact pq
% B_pq =     2Mx6x12 kinematic matrix for each contact.  The 6x12 matrix
%            B_pq(pq,:,:) is the kinematics matrix for contact pq, such that
%              [delta_u_def_pq / delta_theta_def_pq]
%                = [B_pq] * [du_p / du_q / dtheta_p / dtheta_q]
% H_m_max =  the largest value of the H_m_pq_x matrices.  This will be used
%            for determining a tolerance for determinant, eigenvalues, etc.
% H_pq_lambda = 2M x 12 x m array of the lambda stiffnesses for the m lambdas
%            associated with the contact model for this contact
% Q_pq_x_all = the 2M x m x 6 array of consistency dependence on x
% Q_pq_lambda_all = the 2M x 6 x 6 array of consistency dependence on lambdas
% k_lambdas_pq_all = the 2M vector of the number of lamdas associated with
%            the contacts.  For fully elastic contacts, the value is usually 0
% lambdas_pq_all = the 2M x m boolean array of the active lambdas for the
%            contacts
%
% Note that with each pair of touching particles, p and q, we have 
% two variants: a pq variant and a qp variant.  The two variants are stored 
% in separate rows of r_pq_p, n_pq, etc.  Therefore, the number of rows is twice
% the number of contacts.
  if mod(size(r_pq_p,1),2) == 0
    M = size(r_pq_p,1) / 2;
  else
    disp('ERROR in Contact_Stiffness_Matrices');
    disp('r_pq_p should have an even number of rows');
    zzz
  end
%
% maximum number of lambdas among the various contact model
% m_lambdas = 5;
%
% initialize the arrays that contain the mechanical and geometric
% stiffnesses.
  H_m_pq_x =         zeros(2*M,12,12+6);
  H_pq_lambda =      zeros(2*M,12,m_lambdas);
  Q_pq_x_all =       zeros(2*M,m_lambdas,12+6);
  Q_pq_lambda_all =  zeros(2*M,m_lambdas,m_lambdas);
  k_lambdas_pq_all = zeros(2*M,1);
  lambdas_pq_all =   zeros(2*M,m_lambdas);
  H_g1_pq =          zeros(2*M,12,12+6);
  H_g2_pq =          zeros(2*M,12,12+6);
  H_g3_pq =          zeros(2*M,12,12+6);
% B_pq =             zeros(2*M,6,12);   % ZZZ
  B_pq =             zeros(2*M,6,12+6); % 6 additional columns for dxcell
  S_pq_x_all =       zeros(2*M,6,6);
  S_pq_lambda_all =  zeros(2*M,6,m_lambdas);
%
% Note that the input matrices (r_pq_p, n_pq, etc.) are 2*M x 3, with the
% vectors of each contact arranged as the rows of these matrices
% 
% for each contact ...
  for pq = 1:2*M
%
%========= create the mechanical stiffness E_m_pq, =========================
%          as in KDP Eqs. 20-24, for this single contact pq.
%
%   we begin by creating the 6x12 kinematics matrix B_pq for the single pq 
%   contact (see KDP Eqs. 17-19).  This matrix does the following:
%     [delta_u_def_pq / delta_theta_def_pq]
%       = [B_pq] * [du_p / du_q / dtheta_p / dtheta_q]
%
%   Note that the function "Cross_Product" effects the product a X b,
%   so a change in sign is required for the dtheta X r products
%   B_pq(pq,:,:) = [-eye(3),eye(3), ...  % ZZZ
%                   Cross_Product(r_pq_p(pq,:)), -Cross_Product(r_pq_q(pq,:));
%                  zeros(3),zeros(3),-eye(3),eye(3)];
    B_pq(pq,:,:) = ...
      [-eye(3),eye(3), ...
      Cross_Product(r_pq_p(pq,:)), -Cross_Product(r_pq_q(pq,:)), ...
      [dxcell_pq(pq,1), 0, 0, dxcell_pq(pq,4) dxcell_pq(pq,5), 0;
       0, dxcell_pq(pq,2), 0, 0, 0, dxcell_pq(pq,6);
       0, 0, dxcell_pq(pq,3), 0, 0, 0];
      zeros(3),zeros(3),-eye(3),eye(3),zeros(3),zeros(3)];
%
%   the 12x6 statics matrix for the single contact pq, as in KPD Eqs. 7 & 19.
%   However, A_pq is the static matrix that is associated with the single
%   contact, not the entire particle assembly
%   A_pq = squeeze(B_pq(pq,:,:))';  % ZZZ
    
    A_pq = squeeze(B_pq(pq,:,1:12))'; % NOTE: don't include dxcell columns
%
%   the mechanical stiffness matrices for the single contact pq,
%   as in KDP Eqs. 20-23.  The incremental contact force also depends
%   on the activation of slip parameters "lambda"
    [F_pq, M_pq, S_pq_x, S_pq_lambda, Q_pq_x, Q_pq_lambda, lambdas_pq, ...
     k_lambdas_pq] = ...
     LCP3_F_M_Q_pq(FM_model_pq(pq), n_pq(pq,:), f_pq(pq,:), m_pq(pq,:), ...
                  k_pq(pq), alpha_pq(pq), mu_pq(pq), ...
                  lambda_tolerance, other_parameters);
%
    S_pq_x_all(pq,:,:) = S_pq_x;
    S_pq_lambda_all(pq,:,1:k_lambdas_pq) = S_pq_lambda;
%
%========= create the first geometric stiffness E_g1_pq, ===================
%          as in KDP Eq. 28, for the single contact pq.  This will be 
%          done in 5 steps.
%
%   (1) the 3x6 matrix "Dr_pq_n_pre" creates the vector delta_r_n_pq
%   in KDP Eq. 26.  The matrix operates as follows:
%      [delta_r_n_pq] = [Dr_pq_n_pre] * [delta_u_def_pq / delta_theta_def_pq]
%   This matrix is 3 x 6.  Note that the n_pq vector in KDP Eq. 26 is stored
%   in the row(!) pq of the Mx3 matrix n_pq
    Dr_pq_n_pre = 0.5 * n_pq(pq,:)' * [n_pq(pq,:), 0 0 0];
%   Note that the zeros at the end nullify the effect of delta_theta_def_pq
%
%   (2) the 3x6 matrix "Dr_pq_t_pre" that creates the vector delta_r_t_pq
%   in KDP Eq. 27.  The matrix operates as follows:
%      [delta_r_t_pq] = [Dr_t_pre] * [delta_u_def_pq / delta_theta_def_pq]
%   This matrix is 3 x 6.  Note that the n_pq vector in KDP Eq. 27 is stored
%   in the row(!) pq of the Mx3 matrix n_pq. 
%
%   (2a) the 3x3 pseudo-inverse of (K_p + K_q)+ in KDP Eq. 27 is stored in
%   the Mx3x3 array K_inv(pq,:,:).  We must extract the pq
%   part of this array for the particular contact pq
    K_inv_pq = squeeze(K_inv(pq,:,:));
%
%   (2b) the 3x3 curvature tensor K_q in KDP Eq. 27.  This 3x3 matrix is stored
%   in the 4-dimenaional array K(pq,2,:,:).  We must extract the pq part
%   of this array for the "q" (2nd) particle of the pq contact
    K_pq_q = squeeze(K_pq(pq,2,:,:));
%
%   (2c) the 3x6 matrix that will produce the cross-product delta_def_pq X n_pq,
%   inside of KDP Eq. 27.  Note that the function "Cross_Product" effects 
%   the product a X b, so a change in sign is required for the 
%   delta_def_pq X n_pq product
    Cross_with_n_pq = [zeros(3), -Cross_Product(n_pq(pq,:))];
%   Note that the zeros at the start nullify the effect of the
%   delta_u_def_pq
%
%   (2d) finally, the 3x6 matrix "Dr_pq_t_pre" that creates the vector 
%   delta_r_t_pq, as in KDP Eq. 27.  The 3x6 matrix Dr_pq_t_pre, when
%   multiplied by [delta_u_def_pq / delta_theta_def_pq], gives the 
%   3x1 vector delta_r_t_pq
    Dr_pq_t_pre = - K_inv_pq ...
                 * (  Cross_with_n_pq ...
                    - K_pq_q ...
                      * ([eye(3), zeros(3)] ...
                         -  n_pq(pq,:)' * [n_pq(pq,:), 0 0 0]));
%
%   (3) add together Dr_pq_n_pre and Dr_pq_t_pre, as in KDP Eq. 25
    Dr_pq_pre = Dr_pq_n_pre + Dr_pq_t_pre;
%
%   (4) multiply Dr_pq_pre by [B_pq].  In the matrices above, we were
%   multiplying by [delta_u_def_pq / delta_theta_def_pq] instead of
%   multiplying by [du_p / du_q / dtheta_p / dtheta_q].  This will
%   create the 3x(12+6) matrix Dr_pq
    Dr_pq = Dr_pq_pre * squeeze(B_pq(pq,:,:));
%
%   We will need the following matrix later to solve E_g2_pq (see
%   steps 4 & 6)
    Dr_pq_t = Dr_pq_t_pre * squeeze(B_pq(pq,:,:));
%
%   (5) we do the cross product in KDP Eq. 28 to create the
%   6x12 matrix E_g1_pq.  Note that the first 3 rows of this matrix
%   are filled with zeros, since the g-1 stiffness applies only to moments
%   and not to forces (note that the g-1 term is absent in KDP Eq. 15)
    E_g1_pq = [ zeros(3,12+6); ...
               -Cross_Product(f_pq(pq,:)) * Dr_pq];
%   note that the function "Cross_Product" effects the product a X b,
%   so a change in sign is required for Dr_pq X f_pq
%
%========= create the second geometric stiffness E_g2_pq, =================
%          as in KDP Eq. 34, for the single contact pq.  This will be 
%          done in 7 steps.
%
%   in steps (1)-(3) we find the vectors delta_f_hat_pq and delta_m_hat_pq
%   that are defined in KDP Eqs. 31-32.  These two vectors will be constructed
%   with the Dr_pq and delta_theta_def_pq vectors defined above.  Both
%   of these vectors are products of matrices and the vector of the
%   particles' movements, [du_p / du_q / dtheta_p / dtheta_q]
%
%   (1) an intermediate 3x3 matrix Qf that will produce the double 
%       product in the first term of KDP Eq. 31:
%          [Qf]*[delta_n_pq] = f_pq X (delta_n_pq X n_pq)
%
    Qf = Cross_Product(f_pq(pq,:)) * (-Cross_Product(n_pq(pq,:)));
%   We will use Qf in step 4
%
%   likewise for the moment m_pq:
    Qm = Cross_Product(m_pq(pq,:)) * (-Cross_Product(n_pq(pq,:)));
%   We will use Qm in step 6
%
%   (2) an intermediate 3x3 matrix Sf that will produce the second part
%       in KDP Eq. 31:
%         [Sf]*[delta_theta_def_ps]
%            = (delta_theta_def_pq * n_pq) f_pq X n_pq
%   note that n_pq(pq,:) is a row vector
    Sf = (Cross_Product(f_pq(pq,:))*n_pq(pq,:)') * n_pq(pq,:);
%   We will use Sf in step 4
%
%   likewise for the moment m_pq:
    Sm = (Cross_Product(m_pq(pq,:))*n_pq(pq,:)') * n_pq(pq,:);
%   We will use Sm in step 6
%
%   (3) the 3x3 curvature tensor K_p in KDP Eq. 33.  This 3x3 matrix is stored
%   in the 4-dimenaional array K(pq,1,:,:).  We must extract the pq part
%   of this array for the "p" (2nd) particle
    K_pq_p = squeeze(K_pq(pq,1,:,:));
%
%   (4) An intermediate 3x3 matrix T that will produce the following
%   cross-product in KDP Eq. 34:
%     [T]*[delta_f_hat_pq] = r_pq X delta_f_hat_pq
    T = Cross_Product(r_pq_p(pq,:));
%
%   (5) the first 3 rows of the 6x12 matrix E_g2_pq.  Note that
%   delta_theta_def_pq is found with the last 3 rows of the 6x12 matrix B_pq
    E_g2_pq_A = Qf * (-K_pq_p * Dr_pq_t) - 0.5*Sf * squeeze(B_pq(pq,4:6,:));
%
%   (6) the last 3 rows of the 6x12 matrix E_g2_pq
    E_g2_pq_B =  Qm * (-K_pq_p * Dr_pq_t) ...
                      - 0.5*Sm * squeeze(B_pq(pq,4:6,:)) ;%...
%                     + T*E_g2_pq_A;
%
%   (7) the combined 6x12 matrix E_g2_pq in KDP Eq. 34
    E_g2_pq = [E_g2_pq_A; E_g2_pq_B];
%
%========= create the third geometric stiffness E_g3_pq, ===================
%          as in KDP Eq. 35, for the single contact pq
%
    E_g3_pq = ...
      [ zeros(3), zeros(3), -Cross_Product(f_pq(pq,:)), zeros(3,3+6); ...
        zeros(3), zeros(3), -Cross_Product(m_pq(pq,:)), zeros(3,3+6) ];
%
%========= now, we compute the 12x12 "H_pq" counterparts of the ============
%          6x12 "E_pq" matrices.  Later, we will assembly the 
%          H_pq matrices into the global matrix for the entire assembly
%
    Eye = [ -eye(3), zeros(3); ...
           zeros(3),   eye(3); ...
           zeros(3),  -eye(3); ...
             eye(3), zeros(3)];
%
%   as in KDP Eqs. 29 and the second parts of 34 and 35
    H_g1_pq(pq,:,:) =  Eye * E_g1_pq;
    H_g2_pq(pq,:,:) = A_pq * E_g2_pq;
    H_g3_pq(pq,:,:) =  Eye * E_g3_pq;
%
%   the 6x6 mechanical stiffness S_m_pq is composed of F_pq and M_pq parts,
%   which are stacked as S_pq_x = [F_pq / M_pq].  Apply KDP Eq. 24:
    H_m_pq_x(pq,:,:) = A_pq * S_pq_x * squeeze(B_pq(pq,:,:));
%
%   the 12*k contributions to stiffness from slip parameters "lambda".
%   Note that when the contact is fully elastic with no possibility of
%   slip (lambdas_pq = 1, the matrix H_pq_lambda(pq,:,:) reamins filled
%   with zeros.  Note that A_pq is 12x6, and S_pq_lambda is 6xk, where
%   k is the number of active lambdas.  When lambdas_pq=1, then
%   S_pq_lambda is a null matrix
    H_pq_lambda(pq,:,1:k_lambdas_pq) = A_pq * S_pq_lambda;
%
%   all of the k x 6 contributions to consistency
    Q_pq_x_all(pq,1:k_lambdas_pq,:) = Q_pq_x * squeeze(B_pq(pq,:,:));
%
%   all of the k x k contributions to consistency
    Q_pq_lambda_all(pq,1:k_lambdas_pq,1:k_lambdas_pq) = Q_pq_lambda;
%
%   number of lambdas for the contact model of this contact
    k_lambdas_pq_all(pq) = k_lambdas_pq;
%
%   Boolean values of the lambdas for this contact.  That is, whether 
%   the contact is on the verge of sliding with respect to the contact's 
%   various lambdas
    lambdas_pq_all(pq,1:k_lambdas_pq) = lambdas_pq;
%
  end  % end of loop that resolves each contact pq
%
% the largest value of the H_m_pq matrices.  This can be used
% for determining a tolerance for determinant, eigenvalues, etc.
  H_m_max = max(max(max(abs(H_m_pq_x))));
