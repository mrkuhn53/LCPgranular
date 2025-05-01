%
% |+|+|+|+|+|+|+|+|+|  Revision 2   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
function  [H_x, H_m_x, H_g1, H_g2, H_g3, H_g1234, B, ...
           H_lambda, Q_x, Q_lambda, M_lambda, ...
           map_lambda_to_pq_and_qp, map_pq_to_qp] = ...
          LCP_Assemble_H_Q_Matrices( ...
            V_pq, H_m_pq_x, H_g1_pq, H_g2_pq, H_g3_pq, H_g4, ...
            B_pq, ...
            H_pq_lambda, Q_pq_x_all, Q_pq_lambda_all, ...
            k_lambdas_pq_all, lambdas_pq_all, ...
            N, FM_model_pq, lambda_tolerance);
%
% Assemble the 6Nx6N mechanical stiffnesses H_m for the entire assembly,
% as with the "squiggles" in KPD Eqs. 22. Also, assemble the geometric stiffnesses
%  H_g1, H_g2, and H_g3 (KPD Eqs. 28, 29, 33, and 34).  Also, assemble the 
% three matrices H_lambda, Q_x, Q_lambda that are associated with the
% lambdas, in XXX Eqs. XXX.
%
% Revision 2 is customized for treating the problem as a linear complementarity
% problem.
%
% In these comments, reference is made to equations in KPD, which is the
% Kuhn, Prunier, Daouadji paper: "Stiffness pathologies in discrete 
% granular systems: Bifurcation, neutral equilibrium, and instability 
% in the presence of kinematic constraints," Int. J. Num. Anal. Methods Geomech.,
% 2019
%
% The global stiffness matrix H_x is assembled from the 12x12 stiffnesses
% of individual contacts.  These 12x12 matrices are of the form H_pq_x, 
% which effect the following products
%   [db_p / db_q / dw_p / dw_q] = [H_m_pq_x] * [du_p / du_q / dtheta_p / dtheta_q]
% Note that rows 4,5,6 and 10, 11, 12 of [H_m_pq_x] involve the equilibrium 
% of particle q of the pair pq.  Remember that we have also found stiffnesses
% [H_m_pq_x] for both pq and qp variants of a contact (2M, the number of H_m_pq_x
% matrices, is twice the number of contacts).  In assembling the matrices, 
% we will ignore rows 4, 5, 6 and 10, 11, and 12, as these are given in the 
% qp contact as the rows 1, 2, 3 and 7, 8, 9.
%
% The rows and columns of the global stiffnesses are arranged as in
% KPD Eq. 9_1.  All forces are stacked above all moments.  All displacements
% are stacked above all rotations
%
%------------------- INPUT --------------------------------
%  V_pq =    2*M x 2 matrix that gives the assembly topology, where M is
%            the number of contacts (2*M, since there is a pq and a qp
%            variant for each pair of contacting particles).  V_pq(pq,1) 
%            is the "p" particle of "pq" contact; whereas, V_pq(c,2) is 
%            the "q" particle of the "pq" contact.
%  H_m_pq_x =  2*M x 12 x 12 array of mechanical contact stiffness
%            matrices.  
%  H_g1_pq = 2*Mx12x12 array of the "g-1" geometric contact stiffnesses.
%            The 12x12 matrix H_g1_pq(pq,:,:) is the stiffness of contact pq
%  H_g2_pq = 2*Mx12x12 array of the "g-2" geometric contact stiffnesses.
%            The 12x12 matrix H_g1_pq(pq,:,:) is the stiffness of contact pq
%  H_g3_pq = 2*Mx12x12 array of the "g-3" geometric contact stiffnesses.
%            The 12x12 matrix H_g1_pq(pq,:,:) is the stiffness of contact pq
%  H_g4 =    6N x 6N geometric stiffness.  If H_g4 is empty, then it is
%            ignored.
%  B_pq =    2*Mx6x12 array of the kinematic matrices of each contact
%  H_pq_lambda = 2M x 12 x m array of the lambda stiffnesses for the m lambdas
%            associated with the contact model for this contact
%  Q_pq_x_all = the 2M x m x 6 array of consistency dependence on x
%  Q_pq_lambda_all = the 2M x 6 x 6 array of consistency dependence on lambdas
%  k_lambdas_pq_all = the 2M vector of the number of lamdas associated with
%            the contacts.
%  lambdas_pq_all = the 2M x m boolean array of the active lambdas for the
%            contacts. For fully elastic contacts, the value is usually 0
%  N =       scalar number of particles.
%  FM_model_pq =2Mx1 vector of contact force models.  The model for contact
%             pq is in FM_model_pq(pq).  Only one model is currently
%             programmed: the linear-frictional model of Section 2.5, with
%             zero contact moment
%  lambda_tolerance = an kx1 vector that gives the tolerances associated
%             with the ativation of the possible slip conditions (lambdas).
%             For a simple linear-frictional contact, k=1 and the
%             lambda_tolerance is a scalar
%
%------------------- OUTPUT --------------------------------
%  H_x =      6Nx6N global stiffness matrix, KPD Eq. 22, KPD Eqs. 37 & 38
%  H_m_x =    6Nx6N global mechanical stiffness matrix, KPD Eq. 22
%  H_g1 =     6Nx6N global geometric stiffness matrix, KPD Eq. 28
%  H_g2 =     6Nx6N global geometric stiffness matrix, KPD Eq. 33
%  H_g3 =     6Nx6N global geometric stiffness matrix, KPD Eq. 34
%  H_g1234 =  6Nx6N global geometric stiffness matrix, sum of H_g1, H_g2, H_g3,
%             and H_g4.  If H_g4 is empty, then it is ignored.
%  B =        6M*6N kinematics matrix, as shown in KPD Eq. 19.  The rows of B
%             correspond to a stacked vector [du_def/dtheta_def] of the 
%             3M contact relative displacement "du_def" on top of the
%             3M contact relative rotations "dtheta_def".  Matrix B will be
%             needed to test the contact branchs for a 6N vector of particle
%             movements, as in KPD Eq. 86
% H_lambda =  the 6N x M_lambda stiffness matrix associated with the lambdas
% Q_x =       the M_lambda x 6N consistency matrix associated with the dx
% Q_lambda =  the M_lambda x M_lambda consistency matrix associated with the
%             lambdas
% M_lambda =  the number of active lambdas
% map_pq_to_qp = 2M x 1 mapping from the pq variant to the qp variant of the same
%              contact
% map_lambda_to_pq_and_qp = M_lambda x 2 array that maps each lambda to its
%             pq and qp contact variants
% map_contact_to_pq_qp = a mapping from the M contacts to the two pq and
%             variants within the V_pq matrix (the latter gives the p and
%             q particles for a pq or qp contact variant). Note that with the
%             pq contact in the first column, p<q.  With the pq variant in
%             the second column, q<p
%
% Note that with each pair of touching particles, p and q, we have 
% two variants: a pq variant and a qp variant.  The two variants are stored 
% in separate rows of V_pq.  Therefore, the number of rows is twice
% the number of contacts M.
  if mod(size(V_pq,1),2) == 0
    M = size(V_pq,1) / 2;
  else
    disp('ERROR in Contact_Stiffness_Matrices');
    disp('r_pq_p should have an even number of rows');
    zzz
  end
%
% we will later assemble the consistency matrix Q.  Instead of having a
% lambda for both the pq and qp variants of a contact, we will have just one
% one lambda for the contact.  So we must identify the two variants in the
% V_pq matrix.  First, create matrices S1 and S2 with the V_pq and its flipped
% variant:
  [S1,I1] = sortrows(V_pq,[1,2]);
  [S2,I2] = sortrows(V_pq(:,[2 1]),[1,2]);
%
% the two matrices must be equal, otherwise there is an error in V_pq
  if sum(sum(S1~=S2)) ~= 0
    disp('OOPS! There is an error in the input matrix V_pq.')
    disp('There should be a qp variant for every pq variant of a contact')
    zzz
  end
%
% only a subset of all contacts will have active lambdas.  For each of these
% contacts with active lambdas, we will only include one of its variants
% (pq or qp) in the consistency Q matrix (we will include only the pq contact
% with p<q).  We must later identify the qp variant of each pq contact.
% Map the pq to the qp contacts
  map_pq_to_qp = zeros(size(V_pq,2),1); % a 2M x 1 matrix
  map_pq_to_qp(I1) = I2;
%
% the number of active lambdas, half of the combination of the pq and qp variants
  M_lambda = sum(sum(k_lambdas_pq_all)) / 2;
%
% initialize the stiffness matrices
  H_m_x =                   zeros(6*N, 6*N);
  H_lambda =                zeros(6*N,M_lambda);
  Q_x =                     zeros(M_lambda, 6*N);
  Q_lambda =                zeros(M_lambda, M_lambda);
  H_g1 =                    zeros(6*N, 6*N);
  H_g2 =                    zeros(6*N, 6*N);
  H_g3 =                    zeros(6*N, 6*N);
  H_g1234 =                 zeros(6*N, 6*N);
  map_lambda_to_pq_and_qp = zeros(M_lambda,2);
% B =                       zeros(2*M*6,6*N);
%
% reset the number of active lambdas.  Below, we use M_lambda as a counter
  M_lambda =  0;
%
% for each contact, add its 12x12 stiffness contents into the global matrix,
% as in KPD Eqs. 28, 33, and 34.  The contacts are mapped to the particles
% with matrix V_pq, as in KPD Eq. 1_3.
  for pq = 1:2*M
%   the p and q particles that touch at contact pq
    p = V_pq(pq,1);
    q = V_pq(pq,2);
%
%   rows of a 12x12 stiffness that contribute to equilibrium of p
    Rows_pq = [1,2,3, 7,8,9];
%
%   rows of 6Nx6N global stiffness matrix the involves the equilibrium of p.
%   Remember, all of the external forces are stacked on top of all of the
%   external moments, as in KPD Eqs. 9_2 and 9_1.
    Rows =    [      3*(p-1) + [1 2 3], ...  % forces on p
               3*N + 3*(p-1) + [1 2 3]];     % moments on p
%
%   columns of 6Nx6N global stiffness matrix that are multiplied by the
%   movements and rotations of p and q.  Remember that all of the displacements
%   are stacked on top of all of the rotations, as in KPD Eqs. 1_1 and 9_1
    Columns = [      3*(p-1) + [1 2 3], ...  % movement of p
                     3*(q-1) + [1 2 3], ...  % movement of q
               3*N + 3*(p-1) + [1 2 3], ...  % rotation of p
               3*N + 3*(q-1) + [1 2 3]];     % rotation of q
%
%   add the contact stiffness into the global mechanical stiffness
    H_m_x(Rows, Columns) = H_m_x(Rows, Columns) ...
                           + squeeze(H_m_pq_x(pq, Rows_pq, :));
%
%   identify those contacts with active lambdas, and create the H_m_lamda
%   Q_x, and Q_lambda matrices.  First, look at all of the lambdas that
%   apply to this contact
    for k = 1:k_lambdas_pq_all(pq)
%     consider only those lambdas that are active (for example, contacts at
%     the friction limit) and only the single variant (pq or qp) for which
%     p < q
      if lambdas_pq_all(pq,k) %if there is one or more lambdas for this contact
        if V_pq(pq,1) < V_pq(pq,2)
%         the two variants, pq and qp, share the same lambda, 
%         since lambda_pq = lambda_qp.  The total number of active 
%         lambdas within the assembly
          M_lambda = M_lambda + 1;
%
%         map the current M_lambda to the pq and qp contact variants
          map_lambda_to_pq_and_qp(M_lambda,:) = [pq, map_pq_to_qp(pq)];
        end
%
%       the stiffness and consistency matrices
        H_lambda(Rows,M_lambda) = H_lambda(Rows,M_lambda) ...
                                  + squeeze(H_pq_lambda(pq,Rows_pq,k))';
        Q_x(M_lambda,Columns) = squeeze(Q_pq_x_all(pq,k,:))';
        Q_lambda(M_lambda,M_lambda) = squeeze(Q_pq_lambda_all(pq,k,k));
      end
    end
%
%   add the contact stiffness into the global stiffness
    H_g1(Rows, Columns) = H_g1(Rows, Columns) ...
                          + squeeze(H_g1_pq(pq, Rows_pq, :));
    H_g2(Rows, Columns) = H_g2(Rows, Columns) ...
                          + squeeze(H_g2_pq(pq, Rows_pq, :));
    H_g3(Rows, Columns) = H_g3(Rows, Columns) ...
                          + squeeze(H_g3_pq(pq, Rows_pq, :));
%
%   rows and columns of the kinematic matrix "B".  We stack the relative
%   contact displacements above the relative contact rotations, as in KPD Eq. 19
    Rows_B = [      3*(pq-1) + [1 2 3], ...  % relative movements of pq
              6*M + 3*(pq-1) + [1 2 3]];     % relative rotations of pq
%
%   the kinematic matrix
%   B(Rows_B,Columns) = B(Rows_B,Columns) + squeeze(B_pq(pq,:,:));
    B = [];
  end % for pq = 1:2*M
%
% we had originally set M_lambda to the maximum possible number of lambdas
% (to avoid refactoring matrices within a loop).  We must rest Q_x, Q_lambda, 
% and H_lambda to the actual number of active lambdas
  H_lambda = H_lambda(:,1:M_lambda);
  Q_x =      Q_x(1:M_lambda,:);
  Q_lambda = Q_lambda(1:M_lambda,1:M_lambda);
  map_lambda_to_pq_and_qp = map_lambda_to_pq_and_qp(1:M_lambda,:);
%
% note that H_g4 is an input argument of this function
  if isempty(H_g4)
    H_g1234 = H_g1 + H_g2 + H_g3;
  else
    H_g1234 = H_g1 + H_g2 + H_g3 + H_g4;
  end
%
% the total H_x matrix
  H_x = H_m_x + H_g1234;
