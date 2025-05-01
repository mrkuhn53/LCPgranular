%
% |+|+|+|+|+|+|+|+|+|  Revision 3   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
function ...
          [Details, Headings, ...
           M_hat, q_hat, H_x_BD_inv, dp, RCond_BD, ...
           C, P_L, P_Lperp, C_MP_inv, ...
           Pathologies, Available_Branch_Combinations, Condition_numbers, ...
           ...
           ... % from function LCP_Assemble_H_Q_Matrices
           H_x, H_m_x, H_g1, H_g2, H_g3, H_g1234, B, ...
           H_lambda, Q_x, Q_lambda, M_lambda, ...
           map_lambda_to_pq_and_qp, map_pq_to_qp, ...
           ...
           ... % from function LCP_Stiffness_Matrices_pq
           H_m_pq_x, H_g1_pq, H_g2_pq, H_g3_pq, B_pq, H_m_max, ...
           H_m_pq_lambda, Q_pq_x_all, Q_pq_lambda_all, ...
           k_lambdas_pq_all, lambdas_pq_all ...
          ] = ...
         LCP2_Stiffness_Pathologies_BCC( ...
           idim, Problem, N, x, V_pq, r_pq_p, r_pq_q, ...
           n_pq, f_pq, m_pq, rho_pq, K_pq, K_inv, ...
           mu_pq, k_pq, alpha_pq, FM_model_pq, ...
           other_parameters, ...
           dp, Type_Constraint, c, f, dx_c, C, dc, ...
           H_g4, lambda_tolerance, Combinations, ...
           Approx_Singular_Matrix, Cond_numbers)
%
% this function finds all stiffness, consistency, projection, and inverse
% matrices to solve the LCP
%
% In these comments, reference is made to equations in KPD, which is the
% Kuhn, Prunier, Daouadji paper: "Stiffness pathologies in discrete 
% granular systems: Bifurcation, neutral equilibrium, and instability 
% in the presence of kinematic constraints," Int. J. Num. Anal. Methods Geomech.,
% 2019
%
% -------------- INPUT ------------------------
%  idim =        dimension of the problem:  2 or 3
%  Problem =     whether (1) all 2^M_lambda matrix equations are solved, or 
%                        (2) an LCP is solved
%  N =           number of particles, including both disks and platens
%  x =           locations of the centers of the particles, x(N,3).  Similar
%                to "u" in Fig. 1
%  V_pq =        2M x 2 matrix that gives the assembly topology, where M is 
%                the number of contacts.  V_pq(pq,1) is the "p" particle
%                of the "pq" contact; whereas, V_pq(c,2) is the "q" particle 
%                of the "pq" contact.
%  r_pq_p =      2M x 3 array. Row r_pq_p(pq,:) is the vector from particle p
%                to contact pq. Fig. 1.
%  r_pq_q =      2M x 3 array. Row r_pq_q(pq,:) is the vector from particle q
%                to contact pq. Fig. 1.
%  n_pq =        unit normal vector of contact pq, directed outward from 
%                particle p.  That is, outward from particle V_pq(pq,1).  
%                An 2M x 3 array.  Fig. 1
%  f_pq =        force at contact pq, acting upon particle p.  That is, acting
%                upon particle V_pq(pq,1).  An 2M x 3 array, where 2M is 
%                the number of contacts.  Fig. 1
%  m_pq =        moment at contact pq, acting upon particle p.  That is, acting
%                upon particle V_pq(pq,1).  An 2M x 3 array, where M is 
%                the number of contacts.  Fig. 1
%  rho_pq =      the radius of curvature of particles p and q at their
%                contact pq.  An 2M x 2 array. See preface of Section 4.
%  K_pq =        curvature tensors of the two particles, p and q, at their 
%                contact pq. It is an 2M x 2 x 3 x 3 array (a 4D array), 
%                where K_pq(pq,1,:,:) is the 3x3 tensor for particle p of 
%                the contact pq, and K_pq(pq,2,:,:) is the 3x3 tensor 
%                for particle q of contact pq. KPD Eqs. 27 & 32
%  K_inv =       the pseudo-inverse of [K_p + K_q].  An 2Mx3x3 array.
%  mu_pq =       2M x 1 array of the friction coefficients of the contacts.
%                Section 2.5.
%  k_pq =        2M x 1 array of the normal stiffness of the contacts.
%                Section 2.5.
%  alpha_pq =    2M x 1 array of the tangential stiffness coefficients.
%                Section 2.5.
%  FM_model_pq = 2M x 1 array.  The contact model for the force of each 
%                contact.  Only one model is coded in function F_M_pq.m: 
%                  = 1, the linear-frictional model of Section 2.5, with zero
%                       contact moment
%  other_parameters = other contact stiffness parameters.  Currently
%                     ignored
%  dp =          6Nx1 vector of the particle force & moment increments, 
%                with the 3Nx1 vector of force increments stacked on
%                top of the 3Nx1 vector of moment increments, as in
%                KPD Eqs. 1 and 2.
%  Type_Constraint = type of displacement constraints, as explained in
%                    Section 2.4:  either 1, 2, 3, or 4
%  dc =          vector of indices of [dx] that are constrained with
%                Type I constraint.  For other types of constraint, c is
%                ignored. Section 2.4.  Type 1.
%  f =           vector of indicies of [dx] that are not constrained
%                with Type I constraint. Ignored for other constraint types.
%                Section 2.4. Type 1.
%  dx_c =        with Type I constraint, The imposed displacements of the
%                "c" displacements.  Ignored for other constraint types.
%                Section 2.4. Type 1.
%  C =           r x 6N constraint matrix for types II and III constraints,
%                as in Section 2.4, KPD Eqs. 57 and 64.  Ignored for other 
%                constraint types.
%  dc =          rx1 vector of movements in KPD Eq. 64.  Ignored for other
%                constraint types.
%  H_g4 =        6N x 6N geometric stiffness. KPD Eqs. 35 and 44.
%  lambda_tolerance = an kx1 vector that gives the tolerances associated
%                with the ativation of the possible slip conditions (lambdas).
%                For a simple linear-frictional contact, k=1 and the
%                lambda_tolerance is a scalar
%  Combinations = with Problem==1, list of the combinations of contact 
%                stiffness branches that should be investigated.  A zero "0" 
%                or empty vector. "[]" means to investigate all possible
%                branches.  A negative value means to investigate no branches.
%  Approx_Singular_Matrix = 0 means using the determinant for assessing
%                          whether a matrix is singular.  Use null() to
%                          to find null space
%                       ~= 0 use the condition number to determine
%                          whether a matrix is singular.  Use sdv() to
%                          find the approximate null space.  The value
%                          of Approx_Singular_Matrix is the minimum
%                          cond() for assessing singularity
%  Cond_numbers = 0, do not record condition numbers; 1, record them
%                    for each combination
%
% ----------- OUTPUT ------------------------------------------
% NOTE: most 6N indices will be replaced with 3N for 2D problems
% M_hat =    k x k matrix in the LCP:  [q_hat] + [M_hat]*[dlambda] >=0,
%            [dlambda] >= 0, ([q_hat] + [M_hat]*[dlambda])'*[dlambda] = 0
% q_hat =    k x 1 vector in the LCP
% H_x_BD_inv =  6Nx6N Bott-Duffin inverse of H_m_x constrained to L, 
%            with L being the nullspace of C
% dp =       6Nx1 loading vector
% C =        r x 6N constraint matrix
% P_L =      6N x 6N matrix that projects onto the null space of C', KPD Eq. 58
% P_Lperp  = 6N x 6N matrix that projects onto the column space of C', KPD Eq. 59
% C_MP_inv = 6N x r Moore-Penrose inverse of C, near KPD Eqs. 67 & 68
%
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
% H_m_pq_lambda = 2M x 12 x m array of the lambda stiffnesses for the m lambdas
%            associated with the contact model for this contact
% Q_pq_x_all = the 2M x m x 6 array of consistency dependence on x
% Q_pq_lambda_all = the 2M x 6 x 6 array of consistency dependence on lambdas
% k_lambdas_pq_all = the 2M vector of the number of lamdas associated with
%            the contacts.  For fully elastic contacts, the value is usually 0
% lambdas_pq_all = the 2M x m boolean array of the active lambdas for the
%            contacts
%
%
%
%  REMOVE these output items?
%
%  Pathologies = n x 3 matrix of the pathologies that were discovered.
%                Column 1 contains the type of pathology:
%                  1) neutral equilibrium for a legitimate (consistent)
%                     combination of contact stiffness branches
%                  2) a legitimate solution of the stiffness equation.
%                     When more than one of these is reported, then
%                     bifurcation of displacments is possible.
%                  3) negative eigenvalue for a legitimate (consistent)
%                     combination of contact stiffness branches
%                Column 2 contains contains the particular "Combination"
%                  of the contact stiffness branches for which the 
%                  pathology applies.
%                Column 3 contains the "I_2" value of KPD Eq. 95 for
%                  possible bifurcated solutions
% Available_Branch_Combinations = MxC matrix of the combinations of
%                contact stiffness branches.  The "n" column gives the
%                index of the stiffness branch of each contact for 
%                the n'th combination.  The combination that produces
%                a pathology is identified in column 2 of "Pathologies",
%                as explained above.
% Condition_numbers = list of the condition numbers of all combinations
%
% declare these variables to have global scope.  These variables
% will be shared with the function "Domain_Search"
  global Subset s_pq GB SubSpace
%
% initiate the output
  Pathologies = zeros(0,13);
  Available_Branch_Combinations = [];
  RCond_BD = [];
%
% a cell array to store details of results, with these columns:
%  Column 1:  category of pathology
%  Column 2:  Combination
%  Column 3:  s_pq (0=elastic, 1=slip)
%  Column 4:  dx = displacements, for solutions, instabilities, and neu. equil.
%
%  Column 5:  du as a 3*N x 3 array
%  Column 6:  dtheta as a 3*N x 3 array
%
%
%  Solution (posible bifurcations, too)
%    Column 7:  I_2 value for solution
%    Column 12: dy
%
%  Neutral equilibrium
%    Column 8:  Rank_Null_Space (neutral equilibrium), empty=other pathologies
%
%  Instability
%    Column 9:  Number_Neg_Eigenvalues (instability), empty for other pathologies
%    Column 10: Lambda (neg. eigenvalue (instability)), empty=other pathologies
%  Other
%    Column 11:  Strain
%    Column 12:  dy (only for solutions)
%    Column 13:  db (an Nx3 arranfement of dy)
%    Column 14:  dw (an Nx3 arranfement of dy)
%
  Details = cell(500,14);
  nDetails = 0;
%
  Headings = cell(14,1);
  Headings(:) = ...
      {'Category', ...
       'Combination', ...
       's_pq', ...
       'dx', ...
       'du', ...
       'dtheta', ...
       'I_2', ...
       'Rank_Null_Space', ...
       'Number_Neg_Eigenvalues', ...
       'Lambda', ...
       'Strain', ...
       'dy', ...
       'db', ...
       'dw'};
%
  Max_cond_number = 0;
%
% a list of the condition numbers for each combination
  Condition_numbers = [];
%
% create projection matrices for constraint Types II, II, and IV.
% First, initialize them
  R = []; Prr = []; Pnrr = []; P_L = []; P_Lperp = []; C_MP_inv = [];
%
  if Type_Constraint==3
%   projection matrices for projecting movement onto the null space of "C",
%   and the onto the orthogonal complement, as in KPD Eqs. 58 and 59
    [P_L, P_Lperp, C_MP_inv] = LCP_Projections_Types_II_III(C);
  end
%
  if Problem==1
%   dummy values for the output arguments
    M_hat = []; q_hat = []; H_x_BD_inv = []; H_x = []; H_m_x = [];
    H_lambda = []; Q_x = []; Q_lambda = []; M_lambda = [];
    map_lambda_to_pq_and_qp = []; map_pq_to_qp = []; H_m_pq_x = [];
    H_m_pq_lambda = []; Q_pq_x_all = []; Q_pq_lambda_all = [];
    k_lambdas_pq_all = []; lambdas_pq_all = [];
%
%   create the 12x12 stiffness matrices of the individual contacts
    [H_m_pq, H_g1_pq, H_g2_pq, H_g3_pq, B_pq, H_m_max] ...
      = Contact_Stiffness_Matrices ...
          (r_pq_p, r_pq_q, n_pq, f_pq, m_pq, K_pq, K_inv, ...
           mu_pq, k_pq, alpha_pq, other_parameters, FM_model_pq);
%
%   determine the contact yield directions "g_pq" in KPD Eq. 80 for each contact.
%   The row vector g_pq(pq,:) is normal to the yield surface (in displacement-
%   space) for the contact pq.
    g_pq = Yield_Directions_pq(f_pq, m_pq, n_pq, mu_pq, k_pq, alpha_pq, ...
                               FM_model_pq, other_parameters);
%
%   determine which of the contacts have reached the frictional limit and
%   can enable alternative branches of a contact's stiffness matrix.
%   Note that "Slip_Tolerance" is the first element in "lambda_tolerance"
    [Available_Branches_pq, Contact_Can_Slip_pq] = ...
    Test_Available_Branches_pq(n_pq, f_pq, m_pq, mu_pq, FM_model_pq, ...
                               lambda_tolerance(1), other_parameters);
%
%   assemble the set of all combinations of contact branches
    [Available_Branch_Combinations] = ...
       Find_Available_Branch_Combinations(Available_Branches_pq, V_pq);
%
%   assemble the 6Nx6N global stiffnesses H_g1, H_g2, H_g3 for the entire 
%   assembly, as with the "squiggles" in KPD Eqs. 28, 33, and 34.  On the other 
%   hand, the mechanical stiffness H_m will depend upon the particular 
%   combination of branches among the contacts (H_m will be assembled 
%   elsewhere).  We also assemble the global kinematics matrix "B" for 
%   the entire assembly, as this will be needed to test the yielding 
%   of the contacts (as in KPD Eqs. 19 and 86)
    [H_g1, H_g2, H_g3, H_g1234, B] = Assemble_Hg(V_pq, H_g1_pq, H_g2_pq, ...
                                                 H_g3_pq, H_g4, B_pq, N);
%
%   we prepare the product [G][B] that appears in KPD Eq. 86.  This matrix
%   only applies to linear-frictional contacts and will be used to 
%   determine whether contact yielding (slip) is implied by a solution 
%   vector for each conact.  The matrix will be used to test the 
%   consistency of a solution vector with the the stiffness matrix that 
%   was used to find the vector (see Section 2.6 and the first enumerated 
%   item in Section 3)
    GB = GB_Test_Matrix(g_pq, Available_Branches_pq, B, FM_model_pq);
%
%   determine the matrix for finding the "best fit" of the particle
%   displacements to an affine displacement field (strain field)
    [Affine_Matrix] = Affine_Fit(x);
%
%   preparations for 2D problems
    [dp, H_g1, H_g2, H_g3, H_g4, H_g1234, GB, ...
     c, f, C, R, Prr, Pnrr, P_L, P_Lperp, C_MP_inv, Rows, HRows, ...
     Affine_Matrix] = ...
      Convert_2D_or_3D(dp, H_g1, H_g2, H_g3, H_g4, H_g1234, GB, ...
                       c, f, C, R, Prr, Pnrr, idim, P_L, P_Lperp, C_MP_inv, ...
                       Affine_Matrix, Type_Constraint);
%
  elseif Problem==2
%   dummy values for the output arguments
    Pathologies = []; Available_Branch_Combinations = []; Condition_numbers =[];
%
%   create the H and Q stiffness/consistency matrices of the individual contacts
    [H_m_pq_x, H_g1_pq, H_g2_pq, H_g3_pq, B_pq, H_m_max, ...
     H_m_pq_lambda, Q_pq_x_all, Q_pq_lambda_all, ...
     k_lambdas_pq_all, lambdas_pq_all] = ...
    LCP_Stiffness_Matrices_pq( ...
     r_pq_p, r_pq_q, n_pq, f_pq, m_pq, ...
     K_pq, K_inv, mu_pq, k_pq, alpha_pq, other_parameters, ...
     lambda_tolerance, FM_model_pq);
%
%   assemble the 6Nx6N global stiffnesses H_g1, H_g2, H_g3 for the entire 
%   assembly, as with the "squiggles" in KPD Eqs. 28, 33, and 34.  On the other 
%   hand, the mechanical stiffness H_m will depend upon the particular 
%   combination of branches among the contacts (H_m will be assembled 
%   elsewhere).  We also assemble the global kinematics matrix "B" for 
%   the entire assembly, as this will be needed to test the yielding 
%   of the contacts (as in KPD Eqs. 19 and 86)
    [H_x, H_m_x, H_g1, H_g2, H_g3, H_g1234, B, ...
     H_lambda, Q_x, Q_lambda, M_lambda, ...
     map_lambda_to_pq_and_qp, map_pq_to_qp] = ...
    LCP_Assemble_H_Q_Matrices( ...
     V_pq, H_m_pq_x, H_g1_pq, H_g2_pq, H_g3_pq, H_g4, ...
     B_pq, ...
     H_m_pq_lambda, Q_pq_x_all, Q_pq_lambda_all, ...
     k_lambdas_pq_all, lambdas_pq_all, ...
     N, FM_model_pq, lambda_tolerance);
%
%   determine the matrix for finding the "best fit" of the particle
%   displacements to an affine displacement field (strain field)
    [Affine_Matrix] = Affine_Fit(x);
%
%   preparations for 2D problems
    [H_x, H_m_x, H_g1, H_g2, H_g3, H_g4, H_g1234, H_lambda, Q_x, ...
     dp, c, f, C, R, Prr, Pnrr, P_L, P_Lperp, C_MP_inv, Rows, HRows, ...
     Affine_Matrix] = ...
    LCP_Convert_2D_or_3D( ...
      dp, H_x, H_m_x, H_g1, H_g2, H_g3, H_g4, H_g1234, H_lambda, Q_x, ...
      c, f, C, R, Prr, Pnrr, idim, P_L, P_Lperp, C_MP_inv, ...
      Affine_Matrix, Type_Constraint);
%   note that the vector "Rows" applies to both 2D and 3D systems
  end
%
% a tiny tolerance value that will be used for testing eigenvalues,
% determinants, etc.
  Tiny = 1e-9;
  Null_tol =  H_m_max * Tiny;
  S_tol = 10 * Tiny;
  Eig_tol = H_m_max * Tiny;
%
  if Problem==2
%   with the LCP, find the matrix M_hat and vector q_hat
%
    if 1
      [ M_hat, q_hat, H_x_BD_inv, RCond_BD] = ...
      Setup_LCP_BCC( ...
        H_x, H_lambda, dp, P_L, P_Lperp, Q_x, Q_lambda, C_MP_inv, dc);
    elseif 0
%     the Bott-Duffin inverse of H_m_x constrained to L, with L being 
%     the nullspace of C
      H_x_BD_inv = P_L * inv(H_x*P_L + P_Lperp);
%
%     the M and q matrices of the LCP problem
      M_hat = Q_lambda - Q_x * H_x_BD_inv * H_lambda;
%
      q_hat = Q_x * (C_MP_inv - H_x_BD_inv * H_x * C_MP_inv) * dc ...
              + Q_x * H_x_BD_inv * dp;
    end
%
  elseif Problem==1  % skip these lines for the LCP problem
%   whether to investigate Neutral Eqilibrium, as in KPD Section 3.1
    Determine_Neutral_Equilibrium = 1;
%  
%   whether to investigate Bifurcation, as in KPD Section 3.2
    Determine_Bifurcation = 1;
%  
%   whether to investigate Instability of Equilibrium, as in KPD Section 3.3
    Determine_Instability = 1;
%  
%   ========= Consider possible combinations of contact branches =========
%  
    if isempty(Combinations) || (length(Combinations==1) && Combinations==0)
      Combinations = 1:size(Available_Branch_Combinations,2);
    elseif (length(Combinations==1) && Combinations<0)
      Combinations = [];
    end
%  
    if Cond_numbers
      Condition_numbers = zeros(length(Combinations),1);
    end
%
    iCombination = 0;
%  
    for Combination = Combinations
%  
      iCombination = iCombination + 1;
%
%     track status during long runs, with screen output
      if mod(iCombination,5000)==0
        fprintf('Now solving Combination = %6i of %6i combinations\n', ...
                Combination, length(Combinations));
      end
%  
%     Mx1 vector that gives the branch for each contact, for the
%     particular Combination of branches
      Branch_pq = Available_Branch_Combinations(:,Combination);
%  
%     assemble the 6Nx6N global mechanical stiffness H_m_i for the particular
%     branch combination
      H_m_i = Assemble_Hm(V_pq, H_m_pq, Branch_pq, N);
%  
%     extract the rows that apply to the 2D or 3D system
      H_m_i = H_m_i(Rows, Rows);
%  
%     total stiffness matrix for this branch combination, as in KPD Eqs. 38 
%     and 87
      H_i = H_m_i + H_g1234;
%  
%     H_i = H_m_i;
%  
%     some intermediate matrices for investigating stiffness pathologies,
%     Section 3.  Depends on the type of constraint, as in Section 2.4
      if Type_Constraint==1
%       the H_ff sub-matrix in KPD Eq. 54. Used in determining neutral
%       equilibrium for Type I constraint, as in Table 1
        H_ff_i = H_i(f,f);
        H_neu_equil = H_ff_i;
%       minimum rank of null space (for determining neutral equilibrium)
        Min_Rank_Null_Space = 0;
%  
%       the U and V matrices in KPD Eq. 54.  Used in determining bifurcation
%       for Type I constraint, as in Table 2
        U = [-eye(length(c)), H_i(c,f); zeros(length(f),length(c)), H_i(f,f)];
        V = [-H_i(c,c), zeros(length(c),length(f)); -H_i(f,c), eye(length(f))];
%       stifffness matrix used in finding bifurcations (see KPD Eq. 54 and 
%       KPD Table 2)
        H_bifurcation = U;
%       right-hand-side used in finding bifurcations(see KPD Eq. 54 and Table 2)
        dp_bifurcation = V * [dx_c; dp(f)];
%  
%       the symmetric part of H_ff in Table 3 and its explanation.  Used in 
%       determining instability of equilibrium, as in KPD Table 3
        H_ff_i_sym = 0.5 * (H_ff_i + H_ff_i');
%       use this matrix to determine stability
        H_stability = H_ff_i_sym;
      elseif Type_Constraint==2 || Type_Constraint==3
%       stiffness matrix in KPD Eqs. 61 and 68, and in KPD Tables 1 and 2.  
%       Used for determining neutral equilibrium.  Also see KPD Table 1.
        X_i = H_i*P_L + P_Lperp;
        H_neu_equil = X_i;
%       minimum rank of null space (for determining neutral equilibrium)
        Min_Rank_Null_Space = 0;
%  
%       stifffness matrix used in finding bifurcations (see KPD Eqs. 61 and 68
%       and KPD Table 2)
        H_bifurcation = X_i;
%       right-hand-side used in finding bifurcations (KPD Table 2)
        if Type_Constraint==2
%         KPD Eq. 61 and KPD Table 2
          dp_bifurcation = dp;
        elseif Type_Constraint==3
%         KPD Eq. 68 and KPD Table 2
          dn_frak = H_i * C_MP_inv * dc; % C_MP_inv is the Moore-Penrose inverse
          dp_bifurcation = dp - dn_frak;
        end
%  
%       symmetric matrix in Table 3 for determining instability of equilibrium
        H_i_sym = P_L' * 0.5*(H_i + H_i') * P_L;
%
%       one more time, to ensure symmetry
        H_i_sym = 0.5 * (H_i_sym + H_i_sym');
%
%       use this matrix to determine stability
        H_stability = H_i_sym;
      elseif Type_Constraint==4
%       stiffness matrix in KPD Eq. 75 and in KPD Tables 1 and 2.  Used in 
%       determining neutral equilibrium
        Hbb_i = H_i * Pnrr;
        H_neu_equil = Hbb_i + Prr;
%       minimum rank of null space (for determining neutral equilibrium)
        Min_Rank_Null_Space = 0;
%  
%       stifffness matrix used in finding bifurcations (see KPD Eq. 75 
%       and KPD Table 2)
        H_bifurcation = Hbb_i + Prr;
%       right-hand-side used to find bifurcations (see KPD Eq. 75 and 
%       KPD Table 2)
        dp_bifurcation = dp;
%  
%       matrix in Table 3 for determining instability
        Hcal_i = Pnrr' * H_i * Pnrr;
%       symmetric part of the matrix, for determining instability of equilibrium
        Hcal_i_sym = 0.5 * (Hcal_i + Hcal_i');
%       use this matrix to determine stability
        H_stability = Hcal_i_sym;
      end
%  
%     initialize whether the stiffness matrix is singular, etc.
      is_Singular = 0;
%     null space of the stiffness matrix, as appears in Table 1
      Null_Space = [];
%  
%     we now address the three separate pathologies: neutral equilibrium
%     (controlability), bifurcation (path instability), and instability of
%     equilibrium.  Each pathology can involve a different stiffness
%     matrix and right-hand-side of the boundary-value problem (KPD Eq. 10).
%     To investigate a particular pathology, we must also consider the 
%     particular type of displacement constraints, and we must check whether
%     the solution vector [dx] is compatible with the slip condition of
%     each contact.
%  
%     now, the three pathology problems...
%  
%     determine the existance of neutral equilibrium for the particular
%     branch combination, KPD Section 3.1 and KPD Table 1
      if Determine_Neutral_Equilibrium
%       the conditions for neutral equilibrium depend upon the type of
%       displacement constraints (Section 2.4).  Of the three approaches
%       (listed as "a", "b", and "c" in Table 1), we use approach "b"
%  
%       with Types I, II, and III constraint, determining the possibility
%       of neutral equilibrium involves finding the determinant of the
%       stiffness matrix.  With Type IV, we must determine whether the
%       null space is of rank 7 or greater.  We will use a general approach.
%       Rather than finding the determinant of the stiffness matrix,
%       we will find the rank of the null space, which can be applied to
%       all Types of constraint
%  
%       whether to use an approximation of the null space
        if ~Approx_Singular_Matrix
%         use the determinant as an exact indicator of singular condition
%  
%         an orthonormal basis of the null space.  For Type I constraint,
%         the number of rows is the same as the length of "f"
          Null_Space = null(H_neu_equil, Null_tol);
%  
%         dimension of the null space
          Rank_Null_Space = size(Null_Space,2);
        else
%         use the condition number as an approximate indicatior of singular
%         condition
          cond_H_neu_equil = cond(H_neu_equil);
          if cond_H_neu_equil > Approx_Singular_Matrix
%           matrix is nearly singular
            [u,s,v] = svd(H_neu_equil);
%  
%           dimension of the null space
            diag_s = diag(s);
            s_small = find(diag_s < s(1)/Approx_Singular_Matrix);
            Rank_Null_Space = length(s_small);
%  
%           an orthonormal basis of the null space.  For Type I constraint,
%           the number of rows is the same as the length of "f"
            Null_Space = v(:,s_small);
          else
            Rank_Null_Space = 0;
          end
        end
%  
%       if singular, neutral equilibrium is possible, but later we must 
%       check the consistency of the null space with contact stiffness branches
        is_Singular = Rank_Null_Space > Min_Rank_Null_Space;
%  
        if Cond_numbers
          if ~Approx_Singular_Matrix
            Condition_numbers(Combination) = cond(H_neu_equil);
          else
            Condition_numbers(Combination) = cond_H_neu_equil;
          end
        end
%  
        if is_Singular
%         the stiffness matrix is singular, we must now determine whether
%         the null space of the matrix is consistent with the stiffness
%         branches of the contacts
%  
          if Type_Constraint==1
%           include zeros into the null space matrix, as with the "X" operator
%           in the Type I (b) part of Table 1. Note that with Type I constraint,
%           the matrix H_ff, with "f" rows and columns, is singular.  To test
%           the consistency of the null space with the contact branches, we
%           must insert zeros into the "c" rows of the null space matrix.
%           To accomplish this task, we use a temporary matrix, and then
%           we replace Null_Space with the temporary matrix
            temp = zeros(HRows, Rank_Null_Space);
            temp(f, :) = Null_Space;
            temp(c, :) = zeros(length(c), Rank_Null_Space);
            Null_Space = temp;
          elseif Type_Constraint==2 || Type_Constraint==3
%           apply the product [P_L][dz] as in KPD Eq. 61 to find the 
%           displacement vector [dx].  The [dz] that correspond to neutral
%           equilibrium are the columns of "Null_Space".  These columns 
%           will be replaced by their product with [P_L].  The columns 
%           will later be checked for consistency
            Null_Space = P_L * Null_Space;
          elseif Type_Constraint==4
            Null_Space = Pnrr * Null_Space;
          end
%  
%         determine whether the null space intersects the region of [dx]
%         that is compatible with the particular slip conditions of the
%         contacts.  This region is written as \Omega^{i} in the paper,
%         and it is pointed convex polygonal cone that is formed from the
%         half-spaces S_pq in KPD Eq. 79.  The vectors g_pq of KPD Eq. 80 were
%         computed above with the function "Yield_Directions_pq".  The
%         g_pq vectors were multiplied by kinematic matrices B and
%         assembled in the matrix GB in function "GB_Test_Matrix" above,
%         and GB corresponds to the product [G][B] in KPD Eq. 86.
%  
%         for linear-frictional contacts only (KPD Eq. 79). "S_pq" has size Mx1.
%         S_pq(pq) is 0 for a contact pq that has not reached the 
%         friction limit and, therefore, is elastic; S_pq(pq) is 1 for a
%         contact that has reached the friction limit and in which [dx] 
%         produces slip (S_pq >0 in KPD Eq. 79); and S_pq(pq) is -1 for a
%         contact that has reached the friction limit but in which [dx]
%         produces unloading and the contact is elastic.  In this case,
%         the vectors [dx] are the columns of the null space
%         Small rounding errors occur in the computation of the
%         product GB*dx, so replace tiny values  with 0
          S_pq = GB * Null_Space;
          S_pq(find(abs(S_pq) < S_tol)) = 0;
          S_pq = sign(S_pq);
%  
%         for linear-frictional contacts only.  There are two possible
%         branches: 1 and 2 (elastic and slip).  These two values are the 
%         index of the array Hm_i_pq(:,branch,:,:) and of the array
%         Available_Branches_pq{:,branches) and the entries in the columns 
%         of "Available_Branch_Combinations".  By subtracting 1 from 
%         a column, we have a Mx1 column filled with 0's and 1's. As such,
%         s_pq(pq) is 0 if contact pq is assumed to elastic with the
%         combination of elastic contacts and slipping contacts that
%         are being assumed.
          s_pq = Available_Branch_Combinations(:,Combination) - 1;
%  
%         now we check whether the columns of the "Null_Space" are
%         compatible with the assumed branches.  If the null space
%         is a single column (rank 1), checking compatability is 
%         straightforward: we must merely determine whether the
%         column lies inside of \Omega_i.  If the null space is 
%         of higher rank, the problem is non-trivial.  For example,
%         if the null space is of rank 2, then the "Null_Space" matrix
%         is composed of the two orthogonal column (basis) vectors
%         for the null space.  Each of the two vectors can lie outside 
%         of the \Omega_i space, but a linear combination of the two 
%         vectors can lie inside(!) of the space.  This presents a
%         much more difficult situation, which we do not yet consider.
%         We will simply consider the null vectors individually.
%  
%         the contact branches for direction [dx]. elastic = 0(or -1), slip = 1.
%         When S_pw(pq) is 0 or -1, the contact is elastic
          S_pq_forward = S_pq > 0;
%  
%         the contact branches for direction [-dx]. elastic = 0, slip = 1.
%         When S_pw(pq) is 0 or 1, the contact is elastic for the reverse
%         direction
          S_pq_reverse = S_pq < 0;
%  
%         subset of possibly slipping contacts
          Subset = find(Contact_Can_Slip_pq);
%  
          Strain = zeros(9,0);
%  
          if Rank_Null_Space==1
%           if S_pq (the "s" in KPD Eq. 86) matches s_pq (the "s" 
%           in KPD Eq. 88), then the null space lies within \Omega_i.  
%           The same applies to the negative of S_pq (the null space 
%           includes [dx] vectors in opposite directions)
            if all(S_pq_forward(Subset)==s_pq(Subset)) || ...
               all(S_pq_reverse(Subset)==s_pq(Subset))
%             a consistent single solution
              Affine = Affine_Matrix \ Null_Space(1:3*N,:);
              if all(S_pq_forward(Subset)==s_pq(Subset))
                Strain =  Affine(4:12);
              else
                Strain = -Affine(4:12);
              end
              Pathologies = [Pathologies; ...
                             [1, Combination, 0, Rank_Null_Space, Strain']];
              [ nDetails, Details] = ...
              SaveDetails( ...
                nDetails, Details, N, ...
                1, Combination, s_pq, Null_Space, [], Rank_Null_Space, ...
                [], [], Strain, []);
            end
          else
%           we will simply consider the null vectors individually, rather
%           than the entire null space
            Found_Consistent_Basis_Vectors = 0;
            for i = 1:Rank_Null_Space
              if all(S_pq_forward(Subset,i)==s_pq(Subset)) || ...
                 all(S_pq_reverse(Subset,i)==s_pq(Subset))
%  
                Found_Consistent_Basis_Vectors = ...
                  Found_Consistent_Basis_Vectors + 1;
                Affine = Affine_Matrix \ Null_Space(1:3*N,i);
                if all(S_pq_forward(Subset,i)==s_pq(Subset))
                  Strain = [Strain,  Affine(4:12)];
                else
                  Strain = [Strain, -Affine(4:12)];
                end
              end
            end
%  
            if Found_Consistent_Basis_Vectors==0
%             we have a situation in which the null space is of a dimension
%             greater than 1, but none of the basis vectors of the null
%             space is consistent with the slip / no-slip conditions of the
%             contacts.  We will use the function "Domain_Search" to search
%             for a combination of the basis vectors that is consistent.
%  
%             these parameters must be passed to the function "Domain_Search"
%             as global variables
              SubSpace = Null_Space;
%  
%             use fminsearch to search for a consistent solution
              [X, FVAL] = ...
                fminsearch(@Domain_Search, 0.5*ones(Rank_Null_Space,1));
%  
              if FVAL < -Tiny
%               if FVAL < 0, then we have found a combination of the basis
%               vectors of the null space that are consistent with the
%               slip / no-slip conditions fo the contacts
%  
%               normalize the combination vector X to unit length
                X = (1 / sqrt(sum(X.^2))) * X;
%  
                Affine = Affine_Matrix \ (Null_Space(1:3*N,:)*X);
                Strain = Affine(4:12);
%  
%               collect this information, to be returned as output
                Pathologies = [Pathologies; ...
                               [1.8, Combination, 0, Rank_Null_Space, Strain']];
                [ nDetails, Details] = ...
                SaveDetails( ...
                  nDetails, Details, N, ...
                  1.8, Combination, s_pq, Null_Space*X, [], Rank_Null_Space, ...
                  [], [], Strain, []);
%             else
%               no consistent solution was found
%  
%               collect this information, to be returned as output
%               Pathologies = [Pathologies; ...
%                              [1.5, Combination, 0, Rank_Null_Space, Strain']];
              end
            else
%             at least one of the basis vectors is consistent
              for i = 1:Found_Consistent_Basis_Vectors
                Pathologies = [Pathologies; ...
                               [1.9, Combination, 0, Rank_Null_Space, ...
                                Strain(:,i)']];
                [ nDetails, Details] = ...
                SaveDetails( ...
                  nDetails, Details, N, ...
                  1.9, Combination, s_pq, ...
                  Null_Space(:,i), [], Rank_Null_Space, ...
                  [], [], Strain(:,i), []);
              end
            end
          end % if Rank_Null_Space==1 / else
        end % if is_Singular
        [];
      end % if Determine_Neutral_Equilibrium
%  
%     determine the existence of multiple solutions for the given
%     boundary conditions (KPD Section 3.2).  Note that we only consider 
%     equilibrium paths (not neutral equilibrium states)
      if Determine_Bifurcation && ~is_Singular
%       yes, we will investigate the possible existence of multiple solutions
%  
%       ignore cases of zero loading
        if ~all(dp_bifurcation==0)
%         solve the matrix equation using the stiffness matrix and 
%         right-hand-side that was defined above
          H_bifurcation_inv = inv(H_bifurcation);
          dx =  H_bifurcation_inv * dp_bifurcation;
%  
          if Type_Constraint==1
%           determine the value of I_2, as in KPD Eq. 95 and Table 2.  Note that
%               dp_bifurcation = V * [dx_c / dp_f]
%           Also note that "H_bifurcation_inv" is inv(U), and "dp_bifurcation"
%           equals [V]*[dx_c / dp_f]
            I_2 = 0.5 * [dx_c; -dp(f)]' * H_bifurcation_inv * dp_bifurcation;
%  
%           now, create the full [dx] vector so that we can check consistency.
%           Combine the constraint vector dx_c with the solution vector dx_f, 
%           as in Table 2.  Note that with Type I constraint, we had checked
%           for the solution [dx] with the matrix H_bifurcation, which is 
%           the [U] matrix (KPD Eq. 55).  But the solution [dx] that is found
%           above actually represents the stacked vector [dp_c / dx_f] in 
%           KPD Eq. 55.  To test the consistency of the null space with the
%           contact branches, we must re-assembly the actual dx vector.  
%           We will do this by using a temporary vector, and then we will
%           replace [dx] with the temporary vector
            dx_temp = dx;  % save dx, as we will need it to find [dp] later
            temp = zeros(HRows, 1);
            temp(f) = dx(length(c)+1:end); % bottom part of the stacked vector
            temp(c) = dx_c;          % top part of vector on right of KPD Eq. 55
            dx = temp;
%  
          elseif Type_Constraint==2
%           now, extract the [dx] vector so that we can check consistency.
%           Above, we had solved for dx, which represents the vector [dz] in
%           KPD Eq. 60.  We will replace dx with the product [P_L][dz], which
%           will project [dz] onto the nullspace of the constraint matrix [C]
%           as in KPD Eq. 62 and in KPD Table 2
            dx = P_L * dx;
%  
%           determine the value of I_2, as in KPD Eq. 95 and KPD Table 2. 
%           First, the forces of constraint, KPD Eq. 62
            dy = -dp + H_i*dx;
%  
%           now apply KPD Eq. 95
            I_2 = 0.5 * (dy'*dx - dp'*dx);
%  
          elseif Type_Constraint==3
%           above, we had solved for dx, which represents the vector [dz] in
%           KPD Eq. 68.  We will replace dx (actually, dz) with the expression
%           in KPD Eq. 69 and in KPD Table 2
            dx = C_MP_inv*dc + P_L*dx;  % C_MP_inv is the Moore-Penrose inverse
%  
%           determine the value of I_2, as in KPD Eq. 95 and KPD Table 2. 
%           First, the forces of constraint, KPD Eq. 70
            dy = -dp + H_i*dx;
%  
%           now apply KPD Eq. 95
            I_2 = 0.5 * (dy'*dx - dp'*dx);
%  
          elseif Type_Constraint==4
            dx = Pnrr * dx;
            dy = -dp + H_i*dx;
            I_2 = 0.5 * (dy'*dx - dp'*dx);
          end
%  
%         see the explanation of "S_pq" given above
          S_pq = GB*dx;
          S_pq(find(abs(S_pq) < S_tol)) = 0;
          S_pq = sign(S_pq);
%  
%         see the explanation of "s_pq" given above
          s_pq = Available_Branch_Combinations(:,Combination) - 1;
%  
%         the contact branches for direction [dx]. elastic = 0, slip = 1.
%         When S_pw(pq) is 0 or -1, the contact is elastic
          S_pq_forward = S_pq > 0;
%  
%         subset of possibly slipping contacts
          Subset = find(Contact_Can_Slip_pq);
%  
%         a consistent solution
          Affine = Affine_Matrix \ dx(1:3*N);
          Strain = Affine(4:12);
%  
%         if S_pq (the "s" in KPD Eq. 86) matches s_pq (the "s" in KPD Eq. 88),
%         then the solution [dx] lies within \Omega_i.
          if all(S_pq_forward(Subset)==s_pq(Subset))
            Pathologies = [Pathologies;
                           [2, Combination, I_2, 0, Strain']];
            [ nDetails, Details] = ...
            SaveDetails( ...
              nDetails, Details, N, ...
              2, Combination, s_pq, dx, I_2, [], ...
              [], [], Strain, dy);
          end
        end
%
      elseif Determine_Bifurcation & is_Singular
%       we will check for a solution, but the matrix X_i is singular
%
%       ignore cases of zero loading
        if ~all(dp_bifurcation==0)
%         find the generalized Bott-Duffin inverse and solve dz
%         H_GBD_inv = H_bifurcation' * inv(H_bifurcation*H_bifurcation');
%
          [U, S, V] = svd(H_bifurcation);
%
          S_diag = diag(S);
          S_diag_zeros = find(abs(S_diag) < 10*eps);
          S_diag_nonzeros = setdiff([1:length(S_diag)], S_diag_zeros);
%
          S_diag_inv = zeros(size(S_diag));
          S_diag_inv(S_diag_zeros) = 0;
          S_diag_inv(S_diag_nonzeros) = 1 ./ S_diag(S_diag_nonzeros);
%
          S_inv = diag(S_diag_inv);
%
          H_GBD_inv = V * S_inv * U';
%
          dz = H_GBD_inv * dp_bifurcation;
%
          if Type_Constraint==3
%           above, we had solved for dx, which represents the vector [dz] in
%           KPD Eq. 68.  We will replace dx (actually, dz) with the expression
%           in KPD Eq. 69 and in KPD Table 2
            dx = C_MP_inv*dc + P_L*dz;  % C_MP_inv is the Moore-Penrose inverse
%
%           determine the value of I_2, as in KPD Eq. 95 and Table 2. First, 
%           the forces of constraint, KPD Eq. 70
            dy = -dp + H_i*dx;
%
%           now apply KPD Eq. 95
            I_2 = 0.5 * (dy'*dx - dp'*dx);
          end
%
%         see the explanation of "S_pq" given above
          S_pq = GB*dx;
          S_pq(find(abs(S_pq) < S_tol)) = 0;
          S_pq = sign(S_pq);
%
%         see the explanation of "s_pq" given above
          s_pq = Available_Branch_Combinations(:,Combination) - 1;
%
%         the contact branches for direction [dx]. elastic = 0, slip = 1.
%         When S_pw(pq) is 0 or -1, the contact is elastic
          S_pq_forward = S_pq > 0;
%
%         subset of possibly slipping contacts
          Subset = find(Contact_Can_Slip_pq);
%
%         a consistent solution
          Affine = Affine_Matrix \ dx(1:3*N);
          Strain = Affine(4:12);
%
%         if S_pq (the "s" in KPD Eq. 86) matches s_pq (the "s" in KPD Eq. 88),
%         then the solution [dx] lies within \Omega_i.
          if all(S_pq_forward(Subset)==s_pq(Subset))
            Pathologies = [Pathologies;
                           [2.1, Combination, I_2, 0, Strain']];
            [ nDetails, Details] = ...
            SaveDetails( ...
              nDetails, Details, N, ...
              2.1, Combination, s_pq, dx, I_2, [], ...
              [], [], Strain, dy);
          end
        end
        [];
      end % Determine_Bifurcation
%  
%     determine whether instability of equilibrium is possible for the given
%     boundary conditions (Section 3.3)
      if Determine_Instability
%       yes, we will investigate instability of equilibrium
%  
%       find the eigenvectors "Eta" and their eigenvalues "Lambda" for the
%       stiffness matrix that was defined above
        [Eta, Lambda] = eig(balance(H_stability,'noperm'));
%  
%       convert Lambda from a diagonal matrix to a vector
        Lambda = diag(Lambda);
%  
%       save the original Eta
        Eta_save = Eta;
%  
%       we are concerned with finding negative eigenvalues, which correspond
%       to possible instability modes.  Ignore the positive eigenvalues and
%       their eigenvectors
        if any(Lambda < -Eig_tol)
          if Type_Constraint==1
%           include zeros into the eigenvector matrix, as with the "X" operator
%           in the Type I(b) part of Table 3.  Note that with Type I constraint,
%           the matrix H_ff, is used for finding the eigenspace.  To test
%           the consistency of the eigenspace with the contact branches, we
%           must insert zeros into the "c" rows of the eigenvalue matrix.
%           To accomplish this task, we use a temporary matrix, and then
%           we replace "Eta" with the temporary matrix
            temp = zeros(HRows);
            temp(f, f) = Eta;
            Eta = temp;
%  
%           place the eigenvalues Lambda into a vector, so that the rows of Eta
%           correspond to the rows of Lambda
            temp = zeros(HRows, 1);
            temp(f) = Lambda;  % convert Lambda from a matrix to a vector
            Lambda = temp;
          elseif Type_Constraint==2 || Type_Constraint==3
%           the eigenvectors that were found above are the [dz] vectors in
%           Table 3.  To test the consistency of the eigenspace with the 
%           contact branches, we compute [dx] as the product [P_L][dz]
            Eta = P_L * Eta;
          elseif Type_Constraint==4
            Eta = Pnrr * Eta;
          end
%  
%         the steps above have re-constructed the eigenspace and eigenvalues.
%         Now, we will isolate the eigenvectors that have negative eigenvalues.
          Neg_Eigenvalues = find(Lambda < -Eig_tol);
          Neg_Eigenvectors = Eta(:,Neg_Eigenvalues);
%  
          Number_Neg_Eigenvalues = length(Neg_Eigenvalues);
%  
%         see the explanation of "S_pq" given above
          S_pq = GB * Neg_Eigenvectors;
          S_pq(find(abs(S_pq) < S_tol)) = 0;
          S_pq = sign(S_pq);
%  
%         see the explanation of "s_pq" given above
          s_pq = Available_Branch_Combinations(:,Combination) - 1;
%  
%         the contact branches for direction [dx]. elastic = 0, slip = 1.
%         When S_pq(pq) is 0 or -1, the contact is elastic
          S_pq_forward = S_pq > 0;
%  
%         the contact branches for direction [-dx]. elastic = 0, slip = 1.
%         When S_pq(pq) is 0 or 1, the contact is elastic for the reverse
%         direction
          S_pq_reverse = S_pq < 0;
%  
%         subset of possibly slipping contacts
          Subset = find(Contact_Can_Slip_pq);
%  
          Strain = zeros(9,0);
%  
          if Number_Neg_Eigenvalues==1
%           if S_pq (the "s" in KPD Eq. 86) matches s_pq (the "s" in KPD Eq.88),
%           then the null space lies within \Omega_i.  The same applies
%           to the negative of S_pq (the null space includes [dx] 
%           vectors in opposite directions)
            if all(S_pq_forward(Subset)==s_pq(Subset)) || ...
               all(S_pq_reverse(Subset)==s_pq(Subset))
%  
%             a consistent single solution
              Affine = Affine_Matrix \ Neg_Eigenvectors(1:3*N,:);
              if all(S_pq_forward(Subset)==s_pq(Subset))
                Strain =  Affine(4:12);
              else
                Strain = -Affine(4:12);
              end
%  
              Pathologies = [Pathologies; ...
                            [3, Combination, Lambda(Neg_Eigenvalues(1)), 1, ...
                             Strain']];
              [ nDetails, Details] = ...
              SaveDetails( ...
                nDetails, Details, N, ...
                3, Combination, s_pq, Neg_Eigenvectors, [], [], ...
                1, Lambda(Neg_Eigenvalues(1)), Strain, []);
            end
          elseif Number_Neg_Eigenvalues > 1
%           we will simply consider the null vectors individually, rather
%           than the entire null space
            Found_Consistent_Basis_Vectors = 0;
            Consistent_Eigenvectors = [];
            for i = 1:Number_Neg_Eigenvalues
              if all(S_pq_forward(Subset,i)==s_pq(Subset)) || ...
                 all(S_pq_reverse(Subset,i)==s_pq(Subset))
%  
                Found_Consistent_Basis_Vectors = ...
                  Found_Consistent_Basis_Vectors + 1;
                Affine = Affine_Matrix \ Neg_Eigenvectors(1:3*N,i);
                if all(S_pq_forward(Subset,i)==s_pq(Subset))
                  Strain = [Strain,  Affine(4:12)];
                else
                  Strain = [Strain, -Affine(4:12)];
                end
                Consistent_Eigenvectors = ...
                  [Consistent_Eigenvectors,Neg_Eigenvalues(i)];
              end
            end
%  
            if Found_Consistent_Basis_Vectors==0
%             we have a situation in which the negative eigenspace is of 
%             a dimension greater than 1, but none of the basis vectors of 
%             the negative eigenspace is consistent with the slip / no-slip 
%             conditions of the contacts.  We will use the function 
%             "Domain_Search" to search for a combination of the basis 
%             vectors that is consistent.
%  
%             these parameters must be passed to the function "Domain_Search"
%             as global variables
              SubSpace = Neg_Eigenvectors;
%  
%             use fminsearch to search for a consistent solution
              [X, FVAL] = ...
                fminsearch(@Domain_Search, 0.5*ones(Number_Neg_Eigenvalues,1));
%  
              if FVAL < -Tiny
%               if FVAL < 0, then we have found a combination of the basis
%               vectors of the negativ eigenspace that are consistent with the
%               slip / no-slip conditions fo the contacts
%  
%               normalize the combination vector X to unit length
                X = (1 / sqrt(sum(X.^2))) * X;
%  
                Affine = Affine_Matrix \ (SubSpace(1:3*N,:)*X);
                Strain = Affine(4:12);
%  
%               combine the negative eigenvectors
                Eta_neg = Eta_save(:,Neg_Eigenvalues) * X;
%  
%               normalize the eigenvector
                Eta_neg = (1/sum(Eta_neg.^2)) * Eta_neg;
%  
%               collect this information, to be returned as output
                Pathologies = [Pathologies; ...
                              [3.8, Combination, ...
                               Eta_neg'*balance(H_stability)*Eta_neg, ...
                               Number_Neg_Eigenvalues, ...
                               Strain']];
                [ nDetails, Details] = ...
                SaveDetails( ...
                  nDetails, Details, N, ...
                  3.8, Combination, s_pq, SubSpace*X, [], [], ...
                  Number_Neg_Eigenvalues, ...
                  Eta_neg'*balance(H_stability)*Eta_neg, Strain, []);
              end
            else
%             at least one of the basis vectors is consistent
              for i = 1:Found_Consistent_Basis_Vectors
                Pathologies = [Pathologies; ...
                              [3.9, Combination, ...
                               Lambda(Consistent_Eigenvectors(i)), ...
                               Number_Neg_Eigenvalues, ...
                               Strain(:,i)']];
                [ nDetails, Details] = ...
                SaveDetails( ...
                  nDetails, Details, N, ...
                  3.9, Combination, s_pq, ...
                  Eta(:,Consistent_Eigenvectors(i)), [], [], ...
                  Number_Neg_Eigenvalues, ...
                  Lambda(Consistent_Eigenvectors(i)), Strain(:,i), []);
              end
            end
          end
        end % any(Lambda<0)
        [];
      end % Determine_Instability
      [];
    end % Combination
  end % skip these lines
%
% disp(Max_cond_number)
% disp(Max_cond_number_combination)
% disp(size(H_neu_equil))
Details = Details(1:nDetails,:);
