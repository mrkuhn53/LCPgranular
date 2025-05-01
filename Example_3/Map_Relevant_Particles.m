function [ N_All, p_All, N_Rel, p_Rel, N_NonRel, p_NonRel, ...
            Map_p_Rel_to_All, Map_p_All_to_Rel, ...
            Map_p_NonRel_to_All, Map_p_All_to_NonRel, ...
           M_All, pq_All, M_Rel, pq_Rel, M_NonRel, pq_NonRel, ...
            Map_pq_Rel_to_All, Map_pq_NonRel_to_All, ...
           x, V_pq, C_pq, r_pq_p, r_pq_q, n_pq, Overlap_pq, ...
           rho_pq, K_pq, K_inv, mu_pq, k_pq, alpha_pq, FM_model_pq] = ...
         Map_Relevant_Particles( ...
           N, M, Min_Coord_No, V_pq, C_pq, Coord_Number, ...
           x, r_pq_p, r_pq_q, n_pq, Overlap_pq, Map_p_Rel_to_All, ...
           rho_pq, K_pq, K_inv, mu_pq, k_pq, alpha_pq, FM_model_pq)
%
% some particles can be excluded from the analysis.  "Min_Coord_No" gives
% the minimum coordination number for a particle to be included.  We will
% ignore all particles with a Coord_Number < Min_Coord_No
%
% this function give mapping and reverse-mapping for the included
% particles & contacts, and for the excluded ones.
%
% "Map_p_Rel_to_All" maps the 1:N_Rel particles to the 1:N particles;
% p_Rel is a list of the relevant particles
  [Map_p_Rel_to_All] = find(Coord_Number >= Min_Coord_No);
%
% list of all particles
  N_All = N;
  p_All = [1:N_All]';
%
% number and list of relevant particles
  N_Rel = length(Map_p_Rel_to_All);
  p_Rel = [1:N_Rel]';
%
% reverse mapping
  Map_p_All_to_Rel = zeros(N,1);
  Map_p_All_to_Rel(Map_p_Rel_to_All) = p_Rel;
%
% map of all particles to the relevant particles
  [p_NonRel, Map_p_NonRel_to_All] = setdiff(p_All, Map_p_Rel_to_All);
%
% number of non-relevant particles
  N_NonRel = length(p_NonRel);
  p_NonRel = [1:N_NonRel]';
%
% reverse mapping
  Map_p_All_to_NonRel = zeros(N,1);
  Map_p_All_to_NonRel(Map_p_NonRel_to_All) = p_NonRel;
% 
% "Map_pq_Rel_to_All" maps the 1:M_Rel contacts to the 1:M contacts;
% pq_Rel is a list of the relevant contacts
  [Map_pq_Rel_to_All] = find(  Coord_Number(V_pq(:,1)) >= Min_Coord_No ...
                             & Coord_Number(V_pq(:,2)) >= Min_Coord_No);
%
% list of all contacts
  M_All = M;
  pq_All = [1:M_All]';
%
% number and list of relevant contacts
  M_Rel = length(Map_pq_Rel_to_All);
  pq_Rel = [1:M_Rel]';
%
% map of all particles to the relevant particles
  [pq_NonRel, Map_pq_NonRel_to_All] = setdiff(pq_All, Map_pq_Rel_to_All);
%
% number of non-relevant contacts
  M_NonRel = length(pq_NonRel);
%
% return the culled contact arrays
  x = x(Map_p_Rel_to_All,:);
  V_pq = [Map_p_All_to_Rel(V_pq(Map_pq_Rel_to_All,1)), ...
          Map_p_All_to_Rel(V_pq(Map_pq_Rel_to_All,2))];
%
  C_pq = C_pq(Map_pq_Rel_to_All,:);
%
  r_pq_p = r_pq_p(Map_pq_Rel_to_All,:);
  r_pq_q = r_pq_q(Map_pq_Rel_to_All,:);
  n_pq = n_pq(Map_pq_Rel_to_All,:);
  Overlap_pq = Overlap_pq(Map_pq_Rel_to_All,:);
  rho_pq = rho_pq(Map_pq_Rel_to_All,:);
  K_pq = K_pq(Map_pq_Rel_to_All,:,:,:);
  K_inv = K_inv(Map_pq_Rel_to_All,:,:);
  mu_pq = mu_pq(Map_pq_Rel_to_All,:);
  k_pq = k_pq(Map_pq_Rel_to_All,:);
  alpha_pq = alpha_pq(Map_pq_Rel_to_All,:);
  FM_model_pq = FM_model_pq(Map_pq_Rel_to_All,:);
%
