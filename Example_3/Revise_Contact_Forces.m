function [ f_pq, m_pq, b_p, w_p, Warning] = ...
         Revise_Contact_Forces( ...
           N, M, f_pq, m_pq, df_pq, dm_pq, n_pq, ...
           r_pq_p, Overlap_pq, k_pq, V_pq, V_pq_old, ...
           C_pq, C_pq_old, Map_p_Rel_to_All, Map_p_Rel_to_All_old, ...
           FM_model_pq, mu_pq, Time)
%
% this function does the following:
%   1) adds the increments of contact force to compute the new contact forces.
%      For the first time step, the contact forces are computed from the
%      contact overlaps, with the tangential forces equal to zero
%   2) compute the net force and moment on each particle
%
%---------------- INPUT --------------------------
% N =          scalar number of particles
% M =          scalar number of contacts
% f_pq =       Mx3 array of contact forces
% m_pq =       Mx3 array of contact moments
% df_pq =      Mx3 array of contact force increments
% dm_pq =      Mx3 array of contact moment increments
% n_pq =       Mx3 array of contac normals, outward from p
% r_pq_p =     Mx3 array of contact vectors from the center of p to contact pq
% Overlap_pq = Mx1 array of contact overlaps
% k_pq =       Mx1 array of contact normal contact stiffnesses
% V_pq =       Mx2 array of contact topology: p = V_pq(:,1), q = V_pq(:,2)
% V_pq_old     Mx2 array of contact topology
% ---- the following inputs are included to prevent the new contact forces
%      from exceeding the friction limit
% FM_model_pq =M x 1 array.  The contact model for the force of each contact.
%              Only one model is coded in function F_M_pq.m: 
%                = 1, the linear-frictional model of Section 2.5, with zero
%                     contact moments
% mu_pq =      M x 1 array of the friction coefficients of the contacts
%
%-------------- OUTPUT --------------------------
% f_pq = Mx3 array of contact forces
% m_pq = Mx3 array of contact moments
% b_p = Mx3 array of net particle forces
% w_p = Mx3 array of net particle moments
%
  Warning = 0;
%
  if isempty(f_pq)
%   the original contact forces have not yet be computed.  Compute them now
%
%   check whether an error is encountered
    if Time>1
      Warning = 1;
    end
%
%   normal forces
    Fn_pq = -((Overlap_pq.*k_pq)*[1 1 1]) .* n_pq;
%
%   tangential forces, assumed zero
    Ft_pq = zeros(M,3);
%
%   the contact forces at the contacts
    f_pq = Fn_pq + Ft_pq;
%
%   no contact moments at the contacts
    m_pq = zeros(size(f_pq));
%
  else
%   We add increments of contact force to the previous contact forces
%
    f_pq_new = zeros(size(V_pq,1),3);
    m_pq_new = zeros(size(V_pq,1),3);
%
%   the current and previous sets of contacts might be different.  Find the
%   intersection of the two sets
    [C, IA, IB] = intersect([Map_p_Rel_to_All_old(V_pq_old(:,1))', ...
                            Map_p_Rel_to_All_old(V_pq_old(:,2))', C_pq_old], ...
                            [Map_p_Rel_to_All(V_pq(:,1))', ...
                             Map_p_Rel_to_All(V_pq(:,2))', C_pq], 'rows');
%
%   add the increments of contact forces to the shared contacts
    f_pq_new(IB,:) = f_pq(IA,:) + df_pq(IA,:);
    m_pq_new(IB,:) = m_pq(IA,:) + dm_pq(IA,:);
%
%   we might have some new contacts.  Find them, compute the normal force
%   from the contact overlap, and add these forces
    [C,IA] = setdiff([Map_p_Rel_to_All(V_pq(:,1))', ...
                      Map_p_Rel_to_All(V_pq(:,2))', C_pq], ...
                     [Map_p_Rel_to_All_old(V_pq_old(:,1))', ...
                      Map_p_Rel_to_All_old(V_pq_old(:,2))', C_pq_old], ...
                     'rows');
    if ~isempty(IA)
      f_pq_new(IA,:) = -((Overlap_pq(IA).*k_pq(IA))*[1 1 1]) .* n_pq(IA,:);
      m_pq_new(IA,:) = zeros(length(IA),3);
    end
%
%   with the new contact forces, we must check whether contacts have 
%   disengaged and whether the frictional limit has not been reached
    if all(FM_model_pq) == 1
%     the standard linear-frictional contact model
%
%     normal contact force
      fn_new = sum(f_pq_new .* n_pq, 2);
%
%     the dx from the previous time step may have disengaged some of
%     the contacts.  Find these contacts, and set the contact force to
%     zero
      Not_touching = find(fn_new > 0);
      f_pq_new(Not_touching,:) = zeros(length(Not_touching),3);
%
%     check whether the friction limit is exceeded.  The normal and
%     tangential forces.  fn_new should be < 0
      fn_new = sum(f_pq_new .* n_pq, 2);
      ft_new = f_pq_new - (fn_new * [1 1 1]) .* n_pq;
%
%     these are the contacts that have exceeded the friction limit
      Past_limit = find(sqrt(sum(ft_new.^2,2)) > -mu_pq.*fn_new);
%
%     scale back the tangent force, so that it is at the friction limit
      ft_new(Past_limit,:) = ...
        ft_new(Past_limit,:) ...
        .* ...
        ((-mu_pq(Past_limit) .* fn_new(Past_limit) ...
         ./ sqrt(sum(ft_new(Past_limit,:).^2,2)))*[1 1 1]);
%
      f_pq_new(Past_limit,:) = ...
          (fn_new(Past_limit)*[1 1 1]).*n_pq(Past_limit,:) ...
        + ft_new(Past_limit,:);
    else
      disp(cstrcat('The values of FM_model_pq not accomodated in function', ...
                   ' Revise_Contact_Forces'))
      ERROR_FM_model_pq
    end
%
    f_pq = f_pq_new;
    m_pq = m_pq_new;
  end
%
% now that the contact forces have been adjusted, compute the net force
% and moment on the particles
  b_p = zeros(N,3);
  w_p = zeros(N,3);
%
% cross product or contact vector and contact force.  A contributor to the
% net moment on particles (KPD Eq. 6))
  r_cross_f = Cross(r_pq_p, f_pq);
%
% compute the net force on particles
  for pq = 1:size(V_pq,1)
%   for each contact, add its force and moment to the net particle force
%   and moment
    b_p(V_pq(pq,1),:) = b_p(V_pq(pq,1),:) + f_pq(pq,:);
    w_p(V_pq(pq,1),:) = w_p(V_pq(pq,1),:) + m_pq(pq,:) + r_cross_f(pq,:);
  end
%
% update the topology
% V_pq_old = V_pq;
% C_pq_old = C_pq;
