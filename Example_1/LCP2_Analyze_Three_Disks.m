%
% |+|+|+|+|+|+|+|+|+|  Revision 3   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
  clear
% Example with a system of three disks in marginal equilibrium.
%
% The three disks touch at two contacts.  The disks are arranged so that 
% the disk in the middle can be marginally stable: if the top and bottom 
% particles are pressed closer while preventing their rotation, 
% the middle particle will tend to "squirt" away from the top and 
% bottom particles.  The stability depends upon the angle "Beta" 
% in KDP Fig. 7.
%
% The system of particles lies with the x_1 - x_2 plane, x_1 being horizontal.
%
% In comments below, reference is made to equations in KPD, which is the
% Kuhn, Prunier, Daouadji paper: "Stiffness pathologies in discrete 
% granular systems: Bifurcation, neutral equilibrium, and instability 
% in the presence of kinematic constraints," Int. J. Num. Anal. Methods Geomech.,
% 2019
%
% this code works with both Octave and Matlab.  Create the following 
% two constants. The two softwares have subtle differences.
  Octave = 1; Matlab = 2;
%
  addpath('../mfiles2/');
%
% select the software system here
  System = Octave;
%
% problem types, as described in the paper XXX.  Choose 1 or 2. (Problem 3
% is solved in the same manner as Problem 2.)
  Problem = 1; % Try solving the matrix equations for all stiffness matrices
  Problem = 2; % Treat the problem as a linear complementarity problem (LCP)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START providing data for the 3-disk problem
%
% disk radius of top and bottom particles
  Radius1 = 1.;
%
% disk radius of middle particle
% Radius2 = 1;
% for Radius2 = [0.5 1 2]
  for Radius2 = [1]
%
% increase the disks' curvatures at their contacts.  If the factor is
% [], 0, negative, or 1, use a factor of 1
  CurveFactor = [];
% CurveFactor = 2;
%
% friction coefficient "mu", the normal stiffness "k", and the 
% tangential stiffness coefficient "alpha".  These are described in
% KPD Section 2.5
  mu = 0.50;
  k = 1.0;     % normal stiffness = k
  alpha = 1.0; % tangential stiffness = k * alpha
%
% Condition = scalar choise of initial conditions, as follows:
%               = 1, the conditions of KPD Section 4.1, in which the two
%                    contacts are at the sliding limit, with ft = mu*fn
%               = 2, same as Condition=1, except that the initial forces
%                    b_1 and b_3 (see KPD Fig. 7) are vertical. Note that
%                    this condition requires that Beta >= acot(mu)
  Condition = 1;
%
% Note: in the KPD paper, the beta angle is 90 degrees minus the following
% Beta angle
%
% for Beta = [70]           % Case 1 in the paper
% for Beta = [atand(1/mu)]  % Case 2 in the paper
% for Beta = [50]           % Case 3 in the paper
 for Beta = [70, atand(1/mu), 50]
%
% future contact parameters
  other_parameters = [];
%
% normal contact force.  Because the indentations of the particles are
% typically much smaller than the particle size, the normal force is much 
% smaller than than the product Radius*k.  The ratio Fn / (Radius*k) is
% roughly the overlap ratio (overlap / radius)
  Fn = 0.001;
%
% numbers for the various cases in the paper
  Case_A = 1; Case_B = 2; Case_C = 3; Case_D = 4; Case_E = 5; 
  Case_F = 6; Case_other = 7;
%
% the radii of curvature for the two contacts in the 3-particle problem
  Case = Case_D;
%
  switch (Case)
    case (Case_A) % in the paper
      rho_pq = [Radius1, Radius2; ... % p,q = 1,2
                Radius2, Radius1; ... % p,q = 2,1
                Radius2, Radius1; ... % p,q = 2,3
                Radius1, Radius2];    % p,q = 3,2
    case (Case_B) % in the paper
      rho_pq = [  Radius1,   Radius2; ... % p,q = 1,2
                  Radius2,   Radius1; ... % p,q = 2,1
                  Radius2, 2*Radius1; ... % p,q = 2,3
                2*Radius1,   Radius2];    % p,q = 3,2
    case (Case_C) % in the paper
      rho_pq = [    Radius1,     Radius2; ... % p,q = 1,2
                    Radius2,     Radius1; ... % p,q = 2,1
                    Radius2, 0.5*Radius1; ... % p,q = 2,3
                0.5*Radius1,     Radius2];    % p,q = 3,2
    case (Case_D) % in the paper
      rho_pq = [  Radius1,   Radius2; ...     % p,q = 1,2
                  Radius2,   Radius1; ...     % p,q = 2,1
                0.5*Radius2, 0.5*Radius1; ... % p,q = 2,3
                0.5*Radius1, 0.5*Radius2];    % p,q = 3,2
    case (Case_E) % NOT in the paper
      rho_pq = [  Radius1,   Radius2; ... % p,q = 1,2
                  Radius2,   Radius1; ... % p,q = 2,1
                2*Radius2, 2*Radius1; ... % p,q = 2,3
                2*Radius1, 2*Radius2];    % p,q = 3,2
    case (Case_F) % NOT in the paper
      rho_pq = [1.0*Radius1, 0.1*Radius2; ... % p,q = 1,2
                0.1*Radius2, 1.0*Radius1; ... % p,q = 2,1
                0.5*Radius2, 0.5*Radius1; ... % p,q = 2,3
                0.5*Radius1, 0.5*Radius2];    % p,q = 3,2
    case Case_other
      'Your Case here'; ERROR
    otherwise
      'This case is not a proper choice', ERROR
  endswitch
%
% rho_pq = [Radius1, Radius2; ...
%           Radius2, Radius1; ...
%           Radius2, Radius1; ...
%           Radius1, Radius2];
%
% now, we create the "input" information for 3-disk system
  [N, M, ...
   x, V_pq, r_pq_p, r_pq_q, n_pq, f_pq, m_pq, rho_pq, ...
   K_pq, K_inv, ...
   mu_pq, k_pq, alpha_pq, ...
   FM_model_pq, idim] = ...
  LCP2_Input_for_Three_Disks( ...
    Radius1, Radius2, Fn, mu, k, alpha, Beta, CurveFactor, rho_pq, Condition);
%
%mu_pq(3) = 0.501;
%
% Case = 0;
%
% change the curvature of particle 3 at its contact with particle 2
  if 0
  if 0
%   Case A, the base case, in paper
    Case = 1;
  elseif 0
%   Case C, in paper
    Case = 2;
    K_pq(3,2,:,:) = 2*K_pq(3,2,:,:);
    K_pq(4,1,:,:) = 2*K_pq(4,1,:,:);
  elseif 1
%   Case B, in paper
    Case = 3;
    K_pq(3,2,:,:) = 0.5*K_pq(3,2,:,:);
    K_pq(4,1,:,:) = 0.5*K_pq(4,1,:,:);
  elseif 0
%   Case IV
    Case = 4;
    K_pq(3,2,:,:) = 2*K_pq(3,2,:,:);
    K_pq(4,1,:,:) = 2*K_pq(4,1,:,:);
    K_pq(3,1,:,:) = 2*K_pq(3,1,:,:);
    K_pq(4,2,:,:) = 2*K_pq(4,2,:,:);
  elseif 1
%   Case D, in paper
    Case = 5;
    K_pq(3,2,:,:) = 0.5*K_pq(3,2,:,:);
    K_pq(4,1,:,:) = 0.5*K_pq(4,1,:,:);
    K_pq(3,1,:,:) = 0.5*K_pq(3,1,:,:);
    K_pq(4,2,:,:) = 0.5*K_pq(4,2,:,:);
  elseif 0
%   Case VI
    Case = 6;
    K_pa = 2 * K_pq;
  elseif 1
%   Case VII
    Case = 7;
    K_pa = 0.5 * K_pq;
  end
  end
%
% initialize the increments of force and moment on the particles.  The
% vector [dp] is the stacked vector of force increments [db] and moment
% increments [dw], as in KPD Eq. 2_2.  Alternatively, with external follower
% forces, [dp] can represent the lower-dimensional increments in KPD Eq. 5_2
  dp = zeros(6*N,1);
%
% type of constraints, as in KPD Section 2.4.  Note that in the current
% version, only Type 3 is accomodated
  Type_Constraint = 3;
%
% initialize some vectors associated with constraints
  c = []; f = []; dx_c = []; C = []; dc = []; R = []; Prr = []; Pnrr = [];
%
% label the three disks
  Bottom = 1;
  Middle = 2;
  Top = 3;
%
  if Type_Constraint==3
%   with Type III constraint, we specify an r x 6*N constraint matrix C and
%   an r x 1 vector, as in KPD Eq. 64, such that the movements [dx] satisfy 
%   the non-homogeneous condition [C] * [dx] = [dc]
%
%   in this particular example, we specify the rate at which the top
%   and bottom approach each other; zero rotation of the bottom platen; 
%   and zero horizontal movement of the left platen
%
    C = zeros(100,6*N);
    BC = 1;
%
%   matrix of the six constraints
    C(1,      3*(Bottom - 1) + 1) =  1; %  bottom particle, horizontal movement
    C(2,      3*(Bottom - 1) + 2) =  1; %  bottom particle, vertical movement
    C(3,3*N + 3*(Bottom - 1) + 3) =  1; %  bottom particle, theta_3 rotation
    C(4,      3*(Top - 1) + 1) =     1; %  top particle, horizontal movement
    C(5,      3*(Bottom - 1) + 2) = -1; %  top particle, vertical movement
    C(5,      3*(Top - 1) + 2) =     1; %  top particle, vertical movement
    r = 5;
%
    if 1
      C(6,3*N + 3*(Top - 1) + 3) =   1; %  top particle, theta_3 rotation
      r = 6;
%
%     the rate at which the top and bottom particles approach each other
      dc = [ 0; 0; 0; 0; -1; 0]; % movements
    elseif 0
      BC = 2;
      r = 5;
      dc = [ 0; 0; 0; 0; -1]; % movements
    end
  else
    disp('NO CONSTRAINTS HAVE BEEN SPECIFIED')
  end
%
  C = C(1:r,:);
%
% external follower forces
  H_g4 = [];
%
% END providing data for the 3-disk problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% list of the combinations of contact stiffness branches
% that should be investigated.  A zero "0" or empty vector
% [] means to investigate all possible branches.  A negative
% value means to investigate no branches.
  Combinations = [];
%
% do not wait until the end of run to report screen information
  more off; page_output_immediately(1)
%
% lambda_tolerance = an kx1 vector that gives the tolerances associated
% with the activation of possible slip conditions (lambdas).  For a 
% simple linear-frictional contact, k=1 and the lambda_tolerance is a scalar.
% For contact models with multiple slip mechanisms (e.g., both sliding and
% rotational friction), lambda_tolerance will be a vector
  lambda_tolerance = [0.0001];
%
% XXXXX - What are these numbers?
  Approx_Singular_Matrix = 0;
  Approx_Singular_Matrix = 1e6;
  Cond_numbers = 0;
%
   [Details, Headings, ...
    M_hat, q_hat, H_x_BD_inv, dp, RCond_BD, ...
    C, P_L, P_Lperp, C_MP_inv, ...
    Pathologies, Available_Branch_Combinations, Condition_numbers, ...
    ...
    H_x, H_m_x, H_g1, H_g2, H_g3, H_g1234, B, ...
    H_lambda, Q_x, Q_lambda, M_lambda, ...
    map_lambda_to_pq_and_qp, map_pq_to_qp, ...
    ...
    H_m_pq_x, H_g1_pq, H_g2_pq, H_g3_pq, B_pq, H_m_max, ...
    H_m_pq_lambda, Q_pq_x_all, Q_pq_lambda_all, ...
    k_lambdas_pq_all, lambdas_pq_all ...
   ] = ...
  LCP2_Stiffness_Pathologies_Three_Disks( ...
    idim, Problem, N, x, V_pq, r_pq_p, r_pq_q, ...
    n_pq, f_pq, m_pq, rho_pq, K_pq, K_inv, ...
    mu_pq, k_pq, alpha_pq, FM_model_pq, ...
    other_parameters, ...
    dp, Type_Constraint, c, f, dx_c, C, dc, ...
    H_g4, lambda_tolerance, Combinations, ...
    Approx_Singular_Matrix, Cond_numbers);
%
  if Problem==2
%   solve the linear complementarity problem (LCP) for the lambda
%   values.  These lambda values will later be used for finding the dx
%   displacements.
%
%   Path to the directory that holds all of the libraries
    LibrariesPath = '~/LCP/Libraries/';
    addpath(cstrcat(LibrariesPath, 'Murty'));
%
    EPS1 = 1e-11;
%
    if 1
      [ Lambda, Lambda_r, DegenSoln, IsolatedSoln, NonIsolated, ...
        NonIsolatedGrp, NonIsolatedDim, ...
        ColCompetent, DegenMatrix, n_DegenMatrix] = ...
      LCP_enumerate6(M_hat, q_hat, EPS1, 10*eps);
    elseif 0
      [ Lambda, Lambda_r, DegenSoln, IsolatedSoln, NonIsolated,  ...
        NonIsolatedGrp, ...
        ColCompetent, DegenMatrix, n_DegenMatrix] = ...
      LCP_enumerate3(M_hat, q_hat, 10*eps);
    elseif 0
      [ Lambda, DegenSoln, IsolatedSoln, NonIsolated, NonIsolatedGrp, ...
        ColCompetent, DegenMatrix, n_DegenMatrix] = ...
      LCP_enumerate2(M_hat, q_hat, 10*eps);
    elseif 0
      [ Lambda, DegenSoln, IsolatedSoln, NonIsolated, ...
        ColCompetent, DegenMatrix, n_DegenMatrix] = ...
      LCP_enumerate(M_hat,q_hat);
    end
%  
%   number of active contacts
    M_lambda = length(q_hat);
%
%   number of solutions
    nLambda = size(Lambda, 2);
%
%   process the nLambda solutions, as many are duplicates.  Some of the
%   solutions are isolated, others are associated with multi-solution
%   non-isolated solutions.  We remove all of the duplicates, to arrive
%   at the list "Bif_all_list" of the solutions in "Lambda" that are 
%   legitimate
    [ Bif_I, Bif_II, Bif_III, Bif_II_Dim, Bif_III_Dim, ...
      Bif_I_Grp, Bif_II_Grp, Bif_III_Grp, Bif_I_X, ...
      Bif_II_Corners, Bif_III_Corners, ...
      Bif_I_list, Bif_II_III_list, Bif_all_list, Index_Type] = ...
    Count_Bifurcations( ...
      M_hat, q_hat, EPS1, ...
      Lambda, Lambda_r, DegenSoln, IsolatedSoln, NonIsolated, NonIsolatedGrp, ...
      NonIsolatedDim, ColCompetent, DegenMatrix, n_DegenMatrix);
%
%   the number of legitimate solutions
    nBif_all_list = length(Bif_all_list);
%
%   now that we have solved lambda, find the particle movements
    if nLambda > 0
      dx = (  (C_MP_inv - H_x_BD_inv * H_x * C_MP_inv) * dc ...
            + H_x_BD_inv * dp ...
           ) * ones(1,nLambda)  ...
           - H_x_BD_inv * H_lambda * Lambda;
%  
%     reaction forces associated with constrained forces
      dy = - dp * ones(1,nLambda) ...
           + H_x*dx ...
           + H_lambda*Lambda;
%
      if 0
        nLambda2 = 1;
%       nonIs = [4,5];
        nonIs = [3,4];
        Lambda2 = 0.5 * (Lambda(:,nonIs(1)) + Lambda(:,nonIs(2)));
        dx2 = (  (C_MP_inv - H_x_BD_inv * H_x * C_MP_inv) * dc ...
               + H_x_BD_inv * dp ...
              ) * ones(1,nLambda2)  ...
              - H_x_BD_inv * H_lambda * Lambda2;
        dy2 = - dp * ones(1,nLambda2) ...
              + H_x*dx2 ...
              + H_lambda*Lambda2;
        Path_Stability2 = 0.5 * diag(dy2'*dx2 - dp'*dx2)';
      end
%
      if 0
%       the types of complementarity condition for each element if lambda
        Limit1 = 1e-9;
        w = M_hat*Lambda + q_hat*ones(1,nLambda);
        Index_Type = 1000*ones(M_lambda,nLambda);
%
        Index_Type(find((Lambda > Limit1 & w > Limit1) ...
                        | ...
                        (Lambda< -Limit1& w< -Limit1))) = 0; % not a solution
        Index_Type(find(Lambda < Limit1 & w > Limit1)) =  1;
        Index_Type(find(Lambda > Limit1 & w < Limit1)) =  2;
        Index_Type(find(Lambda < Limit1 & w < Limit1)) =  3; % degenerate solution
      end
%
%     du and dtheta
      du =     zeros(nLambda, N, 2);
      dtheta = zeros(nLambda, N, 1);
      db =     zeros(nLambda, N, 2);
      dw =     zeros(nLambda, N, 1);
%
      for i = 1:nLambda
        du(i,:,:) = [dx(1:2:2*N,i), dx(2:2:2*N,i)];
        dtheta(i,:,:) = [dx(2*N+1:1:3*N,i)];
%
        db(i,:,:) = [dy(1:2:2*N,i),     dy(2:2:2*N,i)];
        dw(i,:,:) = [dy(2*N+1:1:3*N,i)];
      end
%
%     compute a parameter to characterize the stability of the current path.
%     The paramter is the second-order work in the direction of the current path
      if ~isempty(Lambda)
%       Path_Stability = dx' * H_x * dx ...
%                        + dx' * H_lambda * Lambda;
%       Path_Stability = 0.5 * diag(dy'*dx - dp'*dx)';
        Path_Stability = 0.5 * sum(dy.*dx - dp.*dx, 1);
      else
        Path_Stability = [];
      end
    else
%     no solutions
      dx = [];
      dy = [];
      du = [];
      db = [];
      dw = [];
      Index_Type = [];
      dtheta = [];
      Path_Stability = [];
    end
%
    disp('Finished running LCP_enumerate3');
%
    if 1
      [ Lambda_Instab, NegEigenGrp, NegEigenValue, NegEigenVector] = ...
        LCP_Instability2(P_L, H_x, H_lambda, Q_x, Q_lambda, 10*eps);
%
%     [ Lambda_Instab, NegEigenGrp, NegEigenValue, NegEigenVector] = ...
%       LCP_Instability(P_L, H_x, H_lambda, Q_x, Q_lambda, 10*eps);
%
      nLambda_Instab = size(Lambda_Instab,2);
%
%     du and dtheta
      du_Instab =     zeros(nLambda_Instab, N, 2);
      dtheta_Instab = zeros(nLambda_Instab, N, 1);
%
      for i = 1:nLambda_Instab
        du_Instab(i,:,:) = [NegEigenVector(1:2:2*N,i), ...
                            NegEigenVector(2:2:2*N,i)];
        dtheta_Instab(i,:,:) = [NegEigenVector(2*N+1:1:3*N,i)];
      end
%
      disp('Finished running LCP_Instability');
    else
      nLambda_Instab = [];
      Lambda_Instab = [];
      NegEigenGrp = [];
      NegEigenValue = [];
      NegEigenVector = [];
      du_Instab = [];
      dtheta_Instab = [];
    end
%
    if 1
      [NeutralVec, NeutralGrp, NeutralRank] = LCP_Neutral2(M_hat, EPS1, 10*eps);
%     [NeutralVec, NeutralGrp, NeutralRank] = LCP_Neutral(M_hat, 10*eps);
      nNeutralVec = size(NeutralVec,2);
%
    disp('Finished running LCP_Neutral');
    else
      nNeutralVec = [];
      NeutralVec = [];
      NeutralGrp = [];
      NeutralRank = [];
    end
%
    if 1
      [Sensitive] = LCP_Sensitive(M_hat, q_hat, Lambda, 10*eps);
%
      disp('Finished running LCP_Sensitive');
    else
      Sensitive = [];
    end
%
    if 1
      [ Lambda_, Lambda_r_, DegenSoln_, IsolatedSoln_, ...
        NonIsolated_, NonIsolatedGrp_, NonIsolatedDim_, ...
        ColCompetent_, DegenMatrix_, n_DegenMatrix_] = ...
      LCP_enumerate6(M_hat, 1e-4*ones(size(M_hat,1),1), EPS1, 10*eps);
%
      is_R_matrix = isempty(Lambda_) || all(all(Lambda_ < 1000*eps));
      is_P_matrix = isPmatrix(M_hat);
      is_N_matrix = isNmatrix(M_hat);
      is_Degen_matrix = n_DegenMatrix>0;
%     is_PD_matrix = all(real(eig(M_hat)) > 10*eps);
      is_PD_matrix = isdefinite(0.5*(M_hat + M_hat'), 10*eps);
      is_ND_matrix = isdefinite(-0.5*(M_hat + M_hat'), 10*eps);
      is_R0_matrix = ~nNeutralVec;
    end
%
    if 0
%     an alternative search of neutral equilibrium
      [ N_Lambda, N_DegenSoln, N_IsolatedSoln, N_NonIsolated, ...
        N_NonIsolatedGrp, N_ColCompetent, N_DegenMatrix, N_n_DegenMatrix] = ...
      LCP_enumerate3(M_hat, zeros(size(q_hat)), 10*eps);
    end
%
    nSensitive_Solutions_I = sum(Sensitive(Bif_I_list));
    nSensitive_Solutions_II_III_set = find(Sensitive(Bif_II_III_list));
    nSensitive_Solutions_II_III_grp = ...
      NonIsolatedGrp( Bif_II_III_list(nSensitive_Solutions_II_III_set));
    nSensitive_Solutions_II_III = ...
      length(unique(nSensitive_Solutions_II_III_grp));
    nSensitive_Solutions = nSensitive_Solutions_I + nSensitive_Solutions_II_III;
%
    [Min_Path_Stability, Stable_Lambda] = min(Path_Stability);
    [Min_NegEigenValue, iMin_NegEigenValue] = min(NegEigenValue);
    Stable_Solns = find(abs(Path_Stability - Min_Path_Stability) ...
                         < 1e-5*abs(Min_Path_Stability));
    mStable = [size(unique(Index_Type(:,Stable_Solns)', 'rows'), 1)];
    Stable_DegenSoln = [DegenSoln(Stable_Lambda)];
%
    Stable_Isolated = [IsolatedSoln(Stable_Lambda)];
    Stable_Sensitive = [Sensitive(:,Stable_Lambda)];
%
    Details =  cell(nLambda,36);
    Headings = cell(36,1);
%
    Details{ 1} = M_lambda;        %  1) M_lambda
    Details{ 2} = M_hat;           %  2) M_hat
    Details{ 3} = q_hat;           %  3) q_hat
    Details{ 4} = nLambda;         %  4) nLambda = number of solutions
    Details{ 5} = Index_Type;      %  5) Index_Type
    Details{ 6} = RCond_BD;        %  6) condition no. of Bott-Duffin inverse
    Details{ 7} = DegenSoln;       %  7) whether a degenerate solution
    Details{ 8} = IsolatedSoln;    %  8) whether an isolated solution
    Details{ 9} = NonIsolated;     %  9) whether a nonisolated solution
    Details(10) = NonIsolatedGrp;  % 10) r for NonIsolatedGrp
    Details{11} = ColCompetent;    % 11) whether a column competent matrix
    Details{12} = DegenMatrix;     % 12) whether a degenerate matrix
    Details{13} = n_DegenMatrix;   % 13) type of degenerate matrix
    Details{14} = dx;              % 14) vectors of dx of solutions
    Details{15} = du;              % 15) vectors of du, 3N x 3 of solutions
    Details{16} = dtheta;          % 16) vectors of dtheta, 3N x 3 of solutions
    Details{17} = dx;              % 17) vectors of dy of solutions
    Details{18} = du;              % 18) vectors of db, 3N x 3 of solutions
    Details{19} = dw;              % 19) vectors of dw, 3N x 3 of solutions
    Details{20} = Path_Stability;  % 20) stability parameter of solutions
    Details{21} = C_MP_inv;        % 21) C_MP_inv
    Details{22} = H_x_BD_inv;      % 22) H_x_BD_inv
    Details{23} = H_x;             % 23) H_x
    Details{24} = H_lambda;        % 24) H_lambda
    Details{25} = nLambda_Instab;  % 25) number of instability modes
    Details{26} = Lambda_Instab;   % 26) Lambda of instability modes
    Details{27} = NegEigenGrp;     % 27) enumeration group, r, of inst. mode
    Details{28} = NegEigenValue;   % 28) neg. eigenvalues of inst. modes
    Details{29} = NegEigenVector;  % 29) neg. eigenvectors of inst. modes
    Details{30} = du_Instab;       % 30) k x N x 3 du's of inst. modes
    Details{31} = dtheta_Instab;   % 31) k x N x 3 dtheta's of inst. modes
    Details{32} = nNeutralVec;     % 32) no. of neutral equilibrium (NE) modes
    Details{33} = NeutralVec;      % 33) neutral equilibrium (NE) vectors
    Details{34} = NeutralGrp;      % 34) enumeration group, r, of (NE) vector
    Details{35} = NeutralRank;     % 35) rank of null space in this r group
    Details{36} = Sensitive;       % 36) whether soln. lambda is path-sensitive
%
    Headings(:) = ...
      {'M_lambda', ...
       'M_hat', ...
       'q_hat', ...
       'nLambda', ...
       'Index_Type', ...
       'RCond_BD', ...
       'DegenSoln', ...
       'IsolatedSoln', ...
       'NonIsolated', ...
       'NonIsolatedGrp', ...
       'ColCompetent', ...
       'DegenMatrix', ...
       'n_DegenMatrix', ...
       'dx', ...
       'du', ...
       'dtheta', ...
       'dx', ...
       'du', ...
       'dtheta', ...
       'Path_Stability', ...
       'C_MP_inv', ...
       'H_x_BD_inv', ...
       'H_x', ...
       'H_lambda', ...
       'nLambda_Instab', ...
       'Lambda_Instab', ...
       'NegEigenGrp', ...
       'NegEigenValue', ...
       'NegEigenVector', ...
       'du_Instab', ...
       'dtheta_Instab', ...
       'nNeutralVec', ...
       'NeutralVec', ...
       'NeutralGrp', ...
       'NeutralRank', ...
       'Sensitive'};
%  
    disp(' ')
    printf('Beta = %5.2f\n', Beta)
    disp(' ')
    nLambda_Instab
    nSensitive_Solutions
%
    disp(' ')
    Stable_DegenSoln
    Stable_Isolated
    mStable
    Stable_Sensitive
%
    disp(' ')
    is_R_matrix
    is_P_matrix
    is_N_matrix
    is_Degen_matrix
    is_PD_matrix
    is_ND_matrix
    eig(0.5*(M_hat + M_hat'))
    is_R0_matrix
%
    disp(' ')
    Bif_I
    Bif_II
    Bif_III
    max_dim_Bif_II = max(Bif_II_Dim)
%
    disp(' ')
%
  elseif Problem==1
%
    disp(' ')
    disp(Radius2)
    disp(Beta)
%
%   if the list of pathologies is small, print to the screen
    if size(Pathologies,1) < 40
      printf(cstrcat(' %0.1f %6i %+12.3e %2i  ', ...
                     '%+0.6f %+0.6f %+0.6f %+0.6f\n'), Pathologies')
    end
%
%   page_output_immediately(0); more on
%
    if 0
%     file name
      Output_File_Name = ...
                         'Results_Analyze_Three_Disks_I';
%                        'Results_Analyze_Three_Disks_III';
%                        'Results_Analyze_Three_Disks_IV';
%
%     save results to a text file
      Output_File_Number = ...
        fopen(cstrcat(Output_File_Name,'.txt'),'w');
      fprintf(Output_File_Number, ...
              cstrcat(' %0.1f %6i %+12.3e %2i  ', ...
                      '%+0.6f %+0.6f %+0.6f %+0.6f\n'), Pathologies')
      fclose(Output_File_Number);
%
%     save to a Matlab-style "mat" file
      save('-mat7-binary', cstrcat(Output_File_Name,'.mat'),'Pathologies');
   end
%
  end  % Problem==1
%
  File_Name = sprintf(cstrcat( ...
                      'Prob_%1i_Rad1_%1.2f_Rad2_%1.2f_Beta_%7.1e_CF_%1.2f',
                      '_Case_%1i_BC_%1i'), ... 
                      Problem, Radius1, Radius2, Beta, CurveFactor, Case, BC);
%
  File_Path = cstrcat('Results/', File_Name, '.mat');
  save -binary Temp_Details ...
       Problem Radius1 Radius2 Beta CurveFactor ...
       C dc K_pq ...
       Fn mu k alpha Fn ...
       Details Headings
  system(cstrcat('mv Temp_Details ',File_Path));


  end  % Beta
  end  % Radius2
