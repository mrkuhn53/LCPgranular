function [ Contact_List, M, x, V_pq, C_pq, r_pq_p, r_pq_q, n_pq, ...
           Overlap_pq, rho_pq, V_pq_old, C_pq_old, Map_p_Rel_to_All_old, ...
           K_pq, K_inv, mu_pq, k_pq, alpha_pq, FM_model_pq, ...
           Coord_Number, idim] = ...
         Contacts_Find( ...
           N, Shape, Size, u, Qp, s_rad, b_rad, xlocal, ...
           xcell, periodic_directions, Near_list_2, Contact_List, ...
           mu, k, alpha, CurveFactor, V_pq, C_pq, Map_p_Rel_to_All)
%
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
% N =         scalar number of particles
% Shape =     scalar particle shape
% Size =      Nx1 vector of particles sizes (e.g., sphere radius)
% u =         Nx3 matrix of particle locations
% QP =        Nx4 matrix of particle orientation quaternions
% s_rad =     for multi-component sphere-clusters, vector of the radii of the 
%             disks/spheres
% b_rad =     for multi-component sphere-clusters, vector of the distance 
%             of centers of spheres from particle center
% xlocal =    for multi-component sphere-clusters, the coordinates of
%             centers of spheres relative to particle center
% xcell =     3x3 matrix of periodic cell size
% periodic_directions = 3x1 vector.  Whether perodic boundaries exist in
%             the three directions
% Near_list_2 = the M_near x 4 matrix of the results from the function
%             Near_Neighbors() of the near-neighbor particle-pairs
% Contact_List = the M x 4 matrix of contacts.  Empty [] for initial assembly
% mu =        friction coefficient, disk-disk contacts
% k =         normal stiffness, disk-disk contacts
% alpha =     coefficient for tangential stiffness, disk-disk contacts
% CurveFactor=increase the radii of cuvatures of the disk-disk contacts by this
%             factor.  If factor is [], 0, 1, or negative, use 1
% V_pq =      old assembly topology (see OUTPUT)
%
% -------------- OUTPUT -------------------------------------------------
% Contact_List = the M x 4 matrix of contacts.  Empty [] for initial assembly
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
% Overlap_pq = Mx1 vector of contact overlaps
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
% Coord_Number = N x 1 vector of the number of contacts of each particle
% idim =      scalar problem dimension, 2 or 3
%
% save previous values
  V_pq_old = V_pq;
  C_pq_old = C_pq;
  Map_p_Rel_to_All_old = Map_p_Rel_to_All;
%
% initialize
  M = []; x = []; V_pq = []; 
  r_pq_p = []; r_pq_q = []; n_pq = []; f_pq = []; m_pq = []; rho_pq = [];
  K_pq = []; K_inv = []; mu_pq = []; k_pq = []; alpha_pq = []; FM_model_pq = [];
%
  [Circle, Oval, Ellipse, Sphere, Ovoid, Nobby, Bumpy] = Shapes3();
%
% whether this is an initial assembly
  Initial_Assembly = isempty(Contact_List);
%
% basic particle shapes
  [Circle, Oval, Ellipse, Sphere, Ovoid, Nobby, Bumpy] = Shapes3();
  if Shape==Circle || Shape==Nobby
    idim = 2;
  elseif Shape==Sphere || Shape==Bumpy
    idim = 3;
  else
    disp('Only disks, spheres, nobbies, or bumpies allowed.')
    ERROR_Create_Assembly_D_1
  end
%
% the positions of the particles, Nx3
  x = u;
%
% now, we have a list of candidates for near-neighbors.  We will cull
% the list, so that the list only includes particles that are within a
% threshold distance of separation
%
% branch vectors between the centers of the two particles
  rxcntr = u(Near_list_2(:,2),:) - u(Near_list_2(:,1),:);
%
  if    (idim==2 && all(periodic_directions(1:2))) ...
     || (idim==3 && all(periodic_directions(1:3)))
%   the number of periodic cell boundaries between the particle pair
%   (for example, when the pair straddles a periodic boundary)
    rcells = (inv(xcell) * rxcntr')';
%
%   the integer number of periodic boundaries crossed by the particle pair.
    nintx = round(rcells);
%
%   adjust the distance between the two particles, adjusting for the 
%   number of periodic boundaries crossed by the particle pair. The
%   branch vector between the centers of the two particles
    rxcntr = rxcntr - (xcell * nintx')';
  end
%
  if Shape==Circle || Shape==Sphere
%   separation between particles in the near-neighbor list
    Separation = sqrt(sum(rxcntr.^2,2)) ...
                 - (Size(Near_list_2(:,1)) + Size(Near_list_2(:,2)));
%
%   list of rows in Near_list_2 in which particles are contacting
    Contacting = find(Separation <= 0);
%
    Contact_List = Near_list_2(Contacting,:);
%
%   number of contacts (half of M)
    nContacts = size(Contacting,1);
%
%   with disks and spheres, there is no offset from the center of the 
%   particle to the center of the disk/sphere
    Offsets_p = zeros(nContacts,3);
    Offsets_q = zeros(nContacts,3);
%
%   separations of all contacts
    Separations = Separation(Contacting);
%
  elseif Shape==Nobby
%   convert quaternions to orientation angles
    Theta = 2*atan2(Qp(:,4), Qp(:,1));
    cosTheta = cos(Theta);
    sinTheta = sin(Theta);
%
%   we must create a number of arrays to hold the data of confirmed contacts.
%   To save memory, reduce the length of these arrays below that of
%   Near_list_2
    Array_Factor = 0.90;
    Array_Length = floor(Array_Factor * size(Near_list_2,1));
%
%   create a shadow list of contacts
    Contact_List = zeros(Array_Length,4);
%
%   initialize list of rows in Near_list_2 in which particles are contacting
    Contacting = zeros(Array_Length,1);
%
%   offset vectors from particle-center to center of the componenent sphere
    Offsets_p = zeros(Array_Length,3);
    Offsets_q = zeros(Array_Length,3);
%
%   total number of contacts
    nContacts = 0;
%
    for pq = 1:size(Near_list_2,1)
%
%     the p and q particles
      p = Near_list_2(pq,1); q = Near_list_2(pq,2);
%
%     the component spheres (nobs) of the p and q sphere-clusters
      pc = Near_list_2(pq,3); qc = Near_list_2(pq,4);
%
%     offset vectors from the centers of p and q to the centers of the
%     component spheres
      Offset_p = (Size(p)*b_rad(pc))*[1 1 1] ...
                 .* ...
                 [  cosTheta(p)*xlocal(1,pc)' ...
                  - sinTheta(p)*xlocal(2,pc)', ...
                    sinTheta(p)*xlocal(1,pc)' ...
                  + cosTheta(p)*xlocal(2,pc)', ...
                  0];
%
      Offset_q = (Size(q)*b_rad(qc))*[1 1 1] ...
                 .* ...
                 [  cosTheta(q)*xlocal(1,qc) ...
                  - sinTheta(q)*xlocal(2,qc), ...
                    sinTheta(q)*xlocal(1,qc) ...
                  + cosTheta(q)*xlocal(2,qc), ...
                  0];
%
%     equal to rxcntr vector from the center of p to the center of q
%     plus vector from center of q to the q-nob minus vector from center
%     of p to p-nob
      rxcntr_2 = rxcntr(pq,:) + Offset_q - Offset_p;
%
%     separation between particles
      Separation = sqrt(sum(rxcntr_2.^2)) ...
                   - (  Size(p)*s_rad(pc) ...
                      + Size(q)*s_rad(qc));
%
%     if the pair is contacting, append to the list of contacts
      if Separation <=0
        nContacts = nContacts + 1;
        Contact_List(nContacts,:) = [p, q, pc, qc];
        Contacting(nContacts) = pq;
        Offsets_p(nContacts,:) = Offset_p;
        Offsets_q(nContacts,:) = Offset_q;
        Separations(nContacts,:) = Separation;
      end
%
    end
%
    Contact_List = Contact_List(1:nContacts,:);
    Contacting = Contacting(1:nContacts,:);
    Offsets_p = Offsets_p(1:nContacts,:);
    Offsets_q = Offsets_q(1:nContacts,:);
    Separations = Separations(1:nContacts,:);
%
  elseif Shape==Bumpy
%
%   rotation matrices of N particles in a 3x3xN array
    [QRot] = RotMatrix_from_Quat(Qp);
%
%   we must create a number of arrays to hold the data of confirmed contacts.
%   To save memory, reduce the length of these arrays below that of
%   Near_list_2
    Array_Factor = 0.50;
    Array_Length = floor(Array_Factor * size(Near_list_2,1));
%
%   create a shadow list of contacts
    Contact_List = zeros(Array_Length,4);
%
%   initialize list of rows in Near_list_2 in which particles are contacting
    Contacting = zeros(Array_Length,1);
%
%   offset vectors from particle-center to center of the componenent sphere
    Offsets_p = zeros(Array_Length,3);
    Offsets_q = zeros(Array_Length,3);
%
%   total number of contacts
    nContacts = 0;
%
    for pq = 1:size(Near_list_2,1)
%     the p and q particles
      p = Near_list_2(pq,1); q = Near_list_2(pq,2);
%
%     the component spheres (nobs) of the p and q sphere-clusters
      pc = Near_list_2(pq,3); qc = Near_list_2(pq,4);



%
%     offset vectors from the centers of p and q to the centers of the
%     component spheres
      Offset_p = (Size(p)*b_rad(pc))*[1 1 1] ...
                 .* ...
                 (QRot(:,:,p)' * xlocal(:,pc))';
%
      Offset_q = (Size(q)*b_rad(qc))*[1 1 1] ...
                 .* ...
                 (QRot(:,:,q)' * xlocal(:,qc))';
%
%     equal to rxcntr vector from the center of p to the center of q
%     plus vector from center of q to the q-nob minus vector from center
%     of p to p-nob
      rxcntr_2 = rxcntr(pq,:) + Offset_q - Offset_p;
%
%     separation between particles
      Separation = sqrt(sum(rxcntr_2.^2)) ...
                   - (  Size(p)*s_rad(pc) ...
                      + Size(q)*s_rad(qc));
%
%     if the pair is contacting, append to the list of contacts
      if Separation <=0
        nContacts = nContacts + 1;
        Contact_List(nContacts,:) = [p, q, pc, qc];
        Contacting(nContacts) = pq;
        Offsets_p(nContacts,:) = Offset_p;
        Offsets_q(nContacts,:) = Offset_q;
        Separations(nContacts,:) = Separation;
      end
%
    end
%
    Contact_List = Contact_List(1:nContacts,:);
    Contacting = Contacting(1:nContacts,:);
    Offsets_p = Offsets_p(1:nContacts,:);
    Offsets_q = Offsets_q(1:nContacts,:);
    Separations = Separations(1:nContacts,:);
  end
%
%------------------------------------------------------------------
% Now, create the contact quantities that will be needed to compute
% the contact and assembly stiffnesses
%
% number of contacts, which includes the pq and qp variants
  M = 2*size(Contact_List,1);

% list of pq and pq contacts
  V_pq = zeros(M,2);
  V_pq(1:2:M,:) = Contact_List(:,[1 2]);
  V_pq(2:2:M,:) = Contact_List(:,[2 1]);
%
% list of pq and qp components
  C_pq = zeros(M,2);
  C_pq(1:2:M,:) = Contact_List(:,[3 4]);
  C_pq(2:2:M,:) = Contact_List(:,[4 3]);
%
% contact frictional and stiffness properties. We allow for different
% properties for the disk-disk contacts
  mu_pq = mu * ones(M,1);
  k_pq = k * ones(M,1);
  alpha_pq = alpha * ones(M,1);
%
% the models for the contact force; and contact moments;
  FM_model_pq = 1*ones(M,1);
%
% find the numbers of contacts for each particle
  List = [Contact_List(:,1); Contact_List(:,2)];
  Coord_Number = hist(List,1:N);
%
% remove List
  List = [];
%
% branch vectors:  from center of p to center of q
  rx = zeros(M,3);
  rx(1:2:M,:) =  rxcntr(Contacting,:);
  rx(2:2:M,:) = -rxcntr(Contacting,:);
%
% offsets for all 2M contacts
  Offsets_p_all = zeros(M,3);
  Offsets_q_all = zeros(M,3);
  Offsets_p_all(1:2:M,:) = Offsets_p;
  Offsets_q_all(1:2:M,:) = Offsets_q;
  Offsets_p_all(2:2:M,:) = Offsets_q;
  Offsets_q_all(2:2:M,:) = Offsets_p;
%
% remove
  Offsets_p = []; Offsets_q = [];
%
% unit normal vectors. First, the vector from the center of the sphere
% component of p to the center of the sphere component of q
  R_pq = rx + Offsets_q_all - Offsets_p_all;
%
% normalize to find the unit normal vectors
  n_pq = R_pq ./ (sqrt(sum(R_pq.^2,2)) * [1 1 1]);
%
% contact vectors
  r_pq_p = Offsets_p_all ...
           + ((s_rad(C_pq(:,1)).*Size(V_pq(:,1)) ...
               ./ (  s_rad(C_pq(:,1)).*Size(V_pq(:,1)) ...
                   + s_rad(C_pq(:,2)).*Size(V_pq(:,2)))) * [1 1 1]) ...
             .* R_pq;
  r_pq_q = r_pq_p - rx;
%
% remove
  Offsets_p_all = []; Offsets_q_all = [];
% remove
  R_pq = []; rx = [];
%
% the vector normal contact forces of the contacts, stored in the rows of
% an Mx3 matrix "Fn_pq". Note that it acts in the direction opposite 
% the contact normal.
%
% if Initial_Assembly
%   with an initial assembly, compute the normal forces using the overlaps
%   at the contacts
    Overlap_pq = zeros(M,1);
    Overlap_pq(1:2:M) = -Separations;
    Overlap_pq(2:2:M) = -Separations;
%
%   normal forces
%   Fn_pq = -((Overlap_pq.*k_pq)*[1 1 1]) .* n_pq;
%
%   tangential forces, assumed zero
%   Ft_pq = zeros(M,3);
% else
%   Overlap_pq = [];
% end
%
% remove
  Separations = [];
%
% the contact forces at the contacts
% f_pq = Fn_pq + Ft_pq;
%
% no contact moments at the contacts
% m_pq = zeros(size(f_pq));
%
% curvatures and curvature tensors of particles p and q at their contact pq
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
% radii of curvatures of particles p and q and contact pq
  rho_pq = CurveFactor ...
           * [s_rad(C_pq(:,1)).*Size(V_pq(:,1)), ...
              s_rad(C_pq(:,2)).*Size(V_pq(:,2))];
%
% we store the curvature tensors in a M x 2 x 3 x 3 array (a 4D array), K_pq.
% Note that M is twice the number of contacts.  See KPD Eqs. 27 and 32
  K_pq = zeros(M,2,3,3);
%
% the pseud0-inverse of [k_p + K_q], a M x 3 x 3 array.  See KPD Eq. 27.
  K_inv = zeros(M,3,3);
%
% Each pq contact has a K_p curvature tensor and a K_q curvature tensor.
% Each of these two tensors is 3 x 3.  Refer to Kuhn & Bagi (2004),
% J. Engrg. Mech., vol. 130, no 7, 826-835.
%
  if idim==2
    p = 1; q = 2;
%
    for pq = 1:M
%     we will use this matrix in our calculations:
      Matrix = [ n_pq(pq,2)*n_pq(pq,2), -n_pq(pq,1)*n_pq(pq,2), 0; ...
                -n_pq(pq,1)*n_pq(pq,2),  n_pq(pq,1)*n_pq(pq,1), 0; ...
                                     0,                      0, 0];
%  
%     the 3x3 K_p tensor for contact pq
      K_pq(pq,p,:,:) = -(1/rho_pq(pq,p)) * Matrix;
%  
%     the 3x3 K_q tensor for contact pq
      K_pq(pq,q,:,:) = -(1/rho_pq(pq,q)) * Matrix;
%  
%     also compute the 3x3 pseudo-inverse of [K_p + K_q], which involves
%     the quotient (rho_p * rho_q) / (rho_p + rho_q), as in eq. 41 of
%     Kuhn and Bagi (2004).
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
%     now, compute the 3x3 pseudo-inverse of [K_p + K_q]
      K_inv(pq,:,:) = - Quotient * Matrix;
    end
  elseif idim==3
    Kp_ap = zeros(3);
    Kq_aq = zeros(3);
%
    Darboux1_p = zeros(M,3);
    Darboux2_p = zeros(M,3);
    Darboux3_p = zeros(M,3);
%
    Darboux1_q = zeros(M,3);
    Darboux2_q = zeros(M,3);
    Darboux3_q = zeros(M,3);
%
    Darboux3_p = n_pq;
%
    avect_o = [0 0 1];
%
    Darboux1_p(:,1) =   avect_o(2).*Darboux3_p(:,3) ...
                      - avect_o(3).*Darboux3_p(:,2);
    Darboux1_p(:,2) =   avect_o(3).*Darboux3_p(:,1) ...
                      - avect_o(1).*Darboux3_p(:,3);
    Darboux1_p(:,3) =   avect_o(1).*Darboux3_p(:,2) ...
                      - avect_o(2).*Darboux3_p(:,1);
%
    Darboux1_p = Darboux1_p.*((sqrt(sum((Darboux1_p.^2)'))').^(-1) *ones(1,3));
%
    Darboux2_p(:,1) =   Darboux3_p(:,2).*Darboux1_p(:,3) ...
                      - Darboux3_p(:,3).*Darboux1_p(:,2);
    Darboux2_p(:,2) =   Darboux3_p(:,3).*Darboux1_p(:,1) ...
                      - Darboux3_p(:,1).*Darboux1_p(:,3);
    Darboux2_p(:,3) =   Darboux3_p(:,1).*Darboux1_p(:,2) ...
                      - Darboux3_p(:,2).*Darboux1_p(:,1);
%
%   Darboux unit vectors for the 2nd (j) particle
%
    Darboux3_q = -n_pq;
%
    Darboux1_q(:,1) =   avect_o(2).*Darboux3_q(:,3) ...
                      - avect_o(3).*Darboux3_q(:,2);
    Darboux1_q(:,2) =   avect_o(3).*Darboux3_q(:,1) ...
                      - avect_o(1).*Darboux3_q(:,3);
    Darboux1_q(:,3) =   avect_o(1).*Darboux3_q(:,2) ...
                      - avect_o(2).*Darboux3_q(:,1);
%
    Darboux1_q = Darboux1_q.*((sqrt(sum((Darboux1_q.^2)'))').^(-1) *ones(1,3));
%
    Darboux2_q(:,1) =   Darboux3_q(:,2).*Darboux1_q(:,3) ...
                      - Darboux3_q(:,3).*Darboux1_q(:,2);
    Darboux2_q(:,2) =   Darboux3_q(:,3).*Darboux1_q(:,1) ...
                      - Darboux3_q(:,1).*Darboux1_q(:,3);
    Darboux2_q(:,3) =   Darboux3_q(:,1).*Darboux1_q(:,2) ...
                      - Darboux3_q(:,2).*Darboux1_q(:,1);
%
    for pq = 1:M
%
      p = V_pq(pq,1);
      q = V_pq(pq,2);
%
      Kp_ap(1,1) = -1/rho_pq(pq,1);
      Kp_ap(2,2) = Kp_ap(1,1);
%
      Kq_aq(1,1) = -1/rho_pq(pq,2);
      Kq_aq(2,2) = Kq_aq(1,1);
%
      T_ap_e = [Darboux1_p(pq,:);
                Darboux2_p(pq,:);
                Darboux3_p(pq,:)];
%
      T_aq_e = [Darboux1_q(pq,:);
                Darboux2_q(pq,:);
                Darboux3_q(pq,:)];
%
      T_aq_ap = T_aq_e * T_ap_e';
%
      Kq_ap = T_aq_ap' * Kq_aq * T_aq_ap;
%
      Kp_Kq_ap = Kp_ap + Kq_ap;
%
%     Kp_Kq_e = T_ap_e' * (Kp_Kq_ap) * T_ap_e;
%
      Kp_Kq_ap_inv = [inv(Kp_Kq_ap(1:2,1:2)) [0;0]; [0 0 0]];
      Kp_Kq_e_inv = T_ap_e' * Kp_Kq_ap_inv * T_ap_e;
%
      Kp_e = T_ap_e' * Kp_ap * T_ap_e;
      Kq_e = T_aq_e' * Kq_aq * T_aq_e;
%
%     the 3x3 K_p tensor for contact pq
      K_pq(pq,1,:,:) = Kp_e;
%  
%     the 3x3 K_q tensor for contact pq
      K_pq(pq,2,:,:) = Kq_e;
%  
%     now, 3x3 pseudo-inverse of [K_p + K_q]
      K_inv(pq,:,:) = Kp_Kq_e_inv;
    end
  end
%
