%
% |+|+|+|+|+|+|+|+|+|  Revision 1   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
function [N, M, ...
          x, V_pq, r_pq_p, r_pq_q, n_pq, f_pq, m_pq, rho_pq, ...
          K_pq, K_inv, ...
          mu_pq, k_pq, alpha_pq, ...
          FM_model_pq, idim] = ...
  LCP2_Input_for_Three_Disks( ...
          Radius1, Radius2, Fn, mu, k, alpha, Beta, CurveFactor, rho_pq, ...
          Condition)
%
% This function creates the input for a system of three disks that touch
% at two contacts.  The conditions include those of Section 4.1 of KPD.
%
% In comments, reference is made to equations in KPD, which is the
% Kuhn, Prunier, Daouadji paper: "Stiffness pathologies in discrete 
% granular systems: Bifurcation, neutral equilibrium, and instability 
% in the presence of kinematic constraints," Int. J. Num. Anal. Methods Geomech.,
% 2019
%
% The system of particles lies with the x_1 - x_2 plane
%
% ------------ INPUT ----------------------------------------------------
% Radius1 =   radius of top and bottom disks
% Radius2 =   radius of middle disk
% Fn =        current normal contact force at the disk-disk contacts
% mu =        friction coefficient, disk-disk contacts
% k =         normal stiffness, disk-disk contacts
% alpha =     coefficient for tangential stiffness, disk-disk contacts
% Beta =      angle of the lower pair from horizontal. Note: in the KPD paper,
%             the beta angle is 90 degrees minus the following Beta angle
%             in this function
% CurveFactor=increase the cuvatures of the disk-disk contacts by this
%             factor.  If factor is [], 0, 1, or negative, use 1
% Condition = scalar choise of initial conditions, as follows:
%               = 1, the conditions of KPD, Section 4.1 in which the two
%                    contacts are at the sliding limit, with ft = mu*fn
%               = 2, same as Condition=1, except that the initial forces
%                    b_1 and b_3 (see KPD Fig. 7) are vertical. Note that
%                    this condition requires that Beta >= acot(mu)
%
% -------------- OUTPUT -------------------------------------------------
% N =         number of particles
% M =         number of contacts.  Note that a pair of touching particles has
%             two contacts(!): a pq variant and a qp variant.  Therefore,
%             M is twice the number of actual contacts (i.e. twice the number
%             of touching pairs of particles)
% x =         locations of the centers of the particles, x(N,3)
% V_pq =      M x 2 matrix that gives the assembly topology, where M is twice
%             the number of contacts.  V_pq(pq,1) is the "p" particle of 
%             "pq" contact; whereas, V_pq(c,2) is the "q" particle of the 
%             "pq" contact.
% r_pq_p =    M x 3 array. Row r_pq_p(pq,:) is the vector from particle p
%             to contact pq
% r_pq_q =    M x 3 array. Row r_pq_q(pq,:) is the vector from particle q
%             to contact pq
% n_pq =      unit normal vector of contact pq, directed outward from particle 
%             p.  That is, outward from particle V_pq(pq,1).  An M x 3 array.
% f_pq =      force at contact pq, acting upon particle p.  That is, acting
%             upon particle V_pq(pq,1).  An M x 3 array, where M is twice
%             the number of contacts.
% rho_pq =    the radius of curvature of particles p and q at their contact pq.
%             An M x 2 array
% K_pq =      curvature tensors of the two particles, p and q, at their 
%             contact pq. It is an M x 2 x 3 x 3 array (a 4D array), where 
%             K_pq(pq,1,:,:) is the 3x3 tensor for particle p of contact pq,
%             and K_pq(pq,2,:,:) is the 3x3 tensor for particle q of 
%             contact pq
% K_inv =     the pseudo-inverse of [K_p + K_q].  An Mx3x3 array.
% mu_pq =     M x 1 array of the friction coefficients of the contacts
% k_pq =      M x 1 array of the normal stiffness of the contacts
% alpha_pq =  M x 1 array of the tangential stiffness coefficients
% FM_model_pq =M x 1 array.  The contact model for the force of each contact.
%             Only one model is coded in function F_M_pq.m: 
%               = 1, the linear-frictional model of Section 2.5, with zero
%                 contact moments
%
% a 2D problem
  idim = 2;
%
%------- create the array of N disks with centers in matrix x(N,3) --------
%
% vector "x", which is the position vector of the three particles' centers.
% We number the particle as follows: 1=bottom, 2 = middle, 3= top
  x = zeros(3,3);
  x(1,:) = [0 0 0];
  x(2,:) = (Radius1 + Radius2)*[cosd(Beta), sind(Beta), 0];
  x(3,:) = 2*(Radius1 + Radius2)*[0, sind(Beta), 0];
%
% number of disks
  N = 3;
%
%
%-----identify the contacts among the disks and build the "v" matrix --------
%
% number of contacts. There are 2 contacts per pair of touching particles,
% since each contact has a pq and qp variant
  M = 4;
% the 4 contacts are as follows:
%   contact 1 = between particles 1 and 2
%   contact 2 = between particles 2 and 1
%   contact 3 = between particles 2 and 3
%   contact 4 = between particles 3 and 2
%
% initialize matrix "V_pq", which is a collapsed incidence matrix that gives
% the topology of the particle assembly.  V_pq(pq,1) is the "p" particle of
% "pq" contact; whereas, V_pq(c,2) is the "q" particle of the "pq" contact.
% If M is the total number of contacts, then "V_pq" will have size 2Mx2,
% since pq and qp are treated as separate entries.  Now, initialize the
% matrix "V_pq"
  V_pq = [1, 2; ...
          2, 1; ...
          2, 3; ...
          3, 2];
%
%----- contact properties of the disks: contact vectors, --------
%      normal vectors, and contact force vectors.
%
% the branch vectors for the contacts, from p to q
  rx = zeros(M,3);
  rx = x(V_pq(:,2),:) - x(V_pq(:,1),:);
%
% the normal vectors of the contacts, outward from p
  n_pq = rx ./ (sqrt(sum(rx.^2, 2)) * [1 1 1]);
%
% the vector normal contact forces of the contacts, stored in the rows of
% an Mx3 matrix "Fn_pq". Note that it acts in the direction opposite 
% the contact normal.
  Fn_pq = -Fn * n_pq;
%
  if Condition==1
%   the contact is ready to slide.  Magnitude of tangential force
    Ft_mag = mu * Fn;
  elseif Condition==2
    Ft_mag = cotd(Beta) * Fn;
  else
    disp('Oops! An invalid value of Condition.')
    ERROR
  end
%
% generic horizontal and vertical components of the tangential contact forces
  ft_x = Ft_mag * sind(Beta);
  ft_y = Ft_mag * cosd(Beta);
%
% the vector tangential contact forces of the contacts, stored in the rows of
% an Mx3 matrix "Ft_pq".
  Ft_pq = [ ft_x, -ft_y,  0; ...
           -ft_x,  ft_y,  0; ...
           -ft_x, -ft_y,  0; ...
            ft_x,  ft_y,  0];
%
% the contact forces at the contacts
  f_pq = Fn_pq + Ft_pq;
%
% no contact moments at the contacts
  m_pq = zeros(size(f_pq));
%
%------curvatures and curvature tensors of particles p and q ---------
%      at their contact pq
%
% radii of curvature of particles p and q.  These radii will be stored in
% the 2M x 2 matrix rho.  rho_pq(pq,1) is the radius of p; rho_pq(pq,2) is the
% radius of q. Thus far, we have only identified contacts between disks.
% The radius of curvature is the disk radius.
% rho_pq = Radius*ones(M,2);
%
  if isempty(CurveFactor) || CurveFactor<=0
    CurveFactor = 1;
  end
%
  radii_pq = [Radius1, Radius2; ...
              Radius2, Radius1; ...
              Radius2, Radius1; ...
              Radius1, Radius2];
%
% the contact vectors for each contact: from p to the contact pq, and
% from q to the contact pq
  r_pq_p =  ((radii_pq(:,1)./(radii_pq(:,1) + radii_pq(:,2)))*[1 1 1]).*rx;
  r_pq_q = -((radii_pq(:,2)./(radii_pq(:,1) + radii_pq(:,2)))*[1 1 1]).*rx;
%
% we store the curvature tensors in a 2M x 2 x 3 x 3 array (a 4D array), K_pq.
% Note that M is twice the number of contacts.  See KPD Eqs. 27 and 32
  K_pq = zeros(M,2,3,3);
%
% the pseud0-inverse of [k_p + K_q], a 2M x 3 x 3 array.  See KPD Eq. 27.
  K_inv = zeros(M,3,3);
%
% Each pq contact has a K_p curvature tensor and a K_q curvature tensor.
% Each of these two tensors is 3 x 3.  Refer to Kuhn & Bagi (2004),
% J. Engrg. Mech., vol. 130, no 7, 826-835.
  p = 1; q = 2;
  for pq = 1:M
%   we will use this matrix in our calculations:
    Matrix = [ n_pq(pq,2)*n_pq(pq,2), -n_pq(pq,1)*n_pq(pq,2), 0; ...
              -n_pq(pq,1)*n_pq(pq,2),  n_pq(pq,1)*n_pq(pq,1), 0; ...
                         0,                      0,           0];
%
%   the 3x3 K_p tensor for contact pq
    K_pq(pq,p,:,:) = -(1/rho_pq(pq,p)) * Matrix;
%
%   the 3x3 K_q tensor for contact pq
    K_pq(pq,q,:,:) = -(1/rho_pq(pq,q)) * Matrix;
%
%   also compute the 3x3 pseudo-inverse of [K_p + K_q], which involves
%   the quotient (rho_p * rho_q) / (rho_p + rho_q), as in eq. 41 of
%   Kuhn and Bagi (2004).
    if      isinf(rho_pq(pq,1)) && ~isinf(rho_pq(pq,2))
      Quotient = rho_pq(pq,2);
    elseif ~isinf(rho_pq(pq,1)) &&  isinf(rho_pq(pq,2))
      Quotient = rho_pq(pq,1);
    elseif ~isinf(rho_pq(pq,1)) && ~isinf(rho_pq(pq,2))
      Quotient = rho_pq(pq,1)*rho_pq(pq,2) / (rho_pq(pq,1) + rho_pq(pq,2));
    else
      Quotient = Inf;
    end
%
%   now, compute the 3x3 pseudo-inverse of [K_p + K_q]
    K_inv(pq,:,:) = - Quotient * Matrix;
  end
%
%-----contact stiffness and friction properties ---------
%
% contact frictional and stiffness properties. We allow for different
% properties for the disk-disk contacts
  mu_pq = mu * ones(M,1);
  k_pq = k * ones(M,1);
  alpha_pq = alpha * ones(M,1);
%
% the models for the contact force; and contact moments;
  FM_model_pq = 1*ones(M,1);
