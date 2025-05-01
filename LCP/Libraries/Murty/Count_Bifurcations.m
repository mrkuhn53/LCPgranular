function [ Bif_I, Bif_II, Bif_III, Bif_II_Dim, Bif_III_Dim, ...
           Bif_I_Grp, Bif_II_Grp, Bif_III_Grp, Bif_I_X, ...
           Bif_II_Corners, Bif_III_Corners, ...
           Bif_I_list, Bif_II_III_list, Bif_all_list, Index_Type] = ...
         Count_Bifurcations( ...
           M, q, EPS1, ...
           X, Xr, DegenSoln, IsolatedSoln, NonIsolated, NonIsolatedGrp, ...
           NonIsolatedDim, ColCompetent, DegenMatrix, n_DegenMatrix);
%
% this function counts the number of Type I, II, and III bifurcations
%
% we eliminate solutions from this list
  Scratch = 1:size(M,1);
%
  Bif_I = 0;    Bif_I_Grp = [];   Bif_I_X = [];
  Bif_II = 0;   Bif_II_Grp = [];  Bif_II_Dim = []; Bif_II_Corners = [];
  Bif_III = 0;  Bif_III_Grp = []; Bif_III_Dim = []; Bif_III_Corners = [];
  Bif_I_list = []; Bif_II_III_list = []; Bif_all_list = [];
  Index_Type = [];
%
  if size(M,1) > 0
%
    Index_Type = zeros(size(X));
    for k = 1:size(X,2)
      x = X(:,k);
      w = M*x + q;
      Index_Type(find(x < EPS1 & w > EPS1),k) = 1;
      Index_Type(find(x > EPS1 & w < EPS1),k) = 2;
      Index_Type(find(x < EPS1 & w < EPS1),k) = 3; % degen. soln.
    end
%
%   begin with the non-isolated solutions of Types II and III bifurcations
    NonIsolated_list = find(~IsolatedSoln);
    NonIsolatedGrp_unique = unique(NonIsolatedGrp(NonIsolated_list));
    Remove = [];
    Remove2 = [];
%
    find_IsolatedSoln = find(IsolatedSoln);
%
    for this_NonIsolatedGrp = NonIsolatedGrp_unique
      this_NonIsolatedGrp_list = find(NonIsolatedGrp == this_NonIsolatedGrp);
      [Y, unique_Index, J] = ...
        unique(Index_Type(:,this_NonIsolatedGrp_list)', 'rows');
      for j = unique_Index'
        Comparator = X(:, this_NonIsolatedGrp_list(j));
        p = find(sum(abs(X - Comparator*ones(1,size(X,2))), 1) < EPS1);
        t = intersect(p, find_IsolatedSoln);
        Remove = [Remove, t];
%
        Bif_II_III_list = [Bif_II_III_list, this_NonIsolatedGrp_list(j)];
      end
%
      Corners = 0;
      for j = this_NonIsolatedGrp_list(unique_Index)
        I2 = find(Index_Type(:,j) == 2);
        Corners = Corners ...
                  + (rank(M(I2,I2), 1e-11)==length(I2) ...
                     && sum(abs(X(I2,j) + M(I2,I2) \ q(I2))) < 1e-11);
      end
%
      if Corners ... % length(setdiff(this_NonIsolatedGrp_list, Remove)) ...
          >= NonIsolatedDim(this_NonIsolatedGrp_list(1))+1
        Bif_II = Bif_II + 1;
        Bif_II_Grp = [Bif_II_Grp, this_NonIsolatedGrp];
        Bif_II_Dim = [Bif_II_Dim, ...
                      NonIsolatedDim(this_NonIsolatedGrp_list(1))];
        Bif_II_Corners = [Bif_II_Corners, Corners];
      elseif (Corners ... % (length(setdiff(this_NonIsolatedGrp_list, Remove)) ...
               < NonIsolatedDim(this_NonIsolatedGrp_list(1))+1) ...
             && length(setdiff(this_NonIsolatedGrp_list, Remove)) > 0
        Bif_III = Bif_III + 1;
        Bif_III_Grp = [Bif_III_Grp, this_NonIsolatedGrp];
        Bif_III_Dim = [Bif_III_Dim, ...
                       NonIsolatedDim(this_NonIsolatedGrp_list(1))];
        Bif_III_Corners = [Bif_III_Corners, Corners];
      end
%
      Remove2 = [Remove2, this_NonIsolatedGrp_list];
    end
%
    Remaining = setdiff(1:size(X,2), union(Remove, Remove2));
    X = X(:, Remaining);
    Xr = Xr(Remaining);
    Index_Type_ = Index_Type(:, Remaining);
    DegenSoln = DegenSoln(Remaining);
    IsolatedSoln = IsolatedSoln(Remaining);
    NonIsolated = NonIsolated(:, Remaining);
    NonIsolatedGrp = NonIsolatedGrp(Remaining);
    NonIsolatedDim = NonIsolatedDim(Remaining);
%
    
    [Y, unique_Index, J] = unique(Index_Type_', 'rows');
    Bif_I = Bif_I + length(unique_Index);
    Bif_I_Grp = Xr(unique_Index);
    Bif_I_X = X(:,unique_Index);
    Bif_I_list = Remaining(unique_Index);
%
    Bif_all_list = [Bif_I_list, Bif_II_III_list];
%
  end
