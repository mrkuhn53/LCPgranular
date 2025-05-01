%
  global M_hat
%
  Get_AltSolutions = 1;
  Get_Pmatrix_list = 1;
  Get_Instab = 1;
  Get_R = 1;
  Get_R0 = 1;
  Get_Sensitive = 1;
  Get_Smatrix_list = 1;
  Get_Li_matrix_list = 1;
  Get_R0_alt = 0;
%
  Create_Plot_File = 1;
%
  Mount_Point = '../';
% Mount_Point = '~/MyDrive/Kuhn_in_Progress/';
% Mount_Point = '~/GoogleDrive/Kuhn_in_Progress/';
% Mount_Point = '~/SSHFS4/';
%
  Output_Directory = cstrcat(Mount_Point, 'Instability/Output/');
%
  RunFile = 'Biaxial_Comp_i';
  RunFile = 'Biaxial_Comp_k2';
  RunFile = 'Biaxial_Comp_k3';
  RunFile = 'Biaxial_Comp_k4';
  RunFile = 'Biaxial_Comp_k5'; %%  with DFile = 'Dcircls_49_a-57';
% RunFile = 'Biaxial_Comp_k6'; %%  with DFile = 'Dcircls_49_a-57';
% RunFile = 'Biaxial_Comp_k7'; %   with DFile = 'Dcircls_49_a-57';
  RunFile = 'Biaxial_Comp_k8'; %%% with DFile = 'Dcircls_49_a-57';
  RunFile = 'Biaxial_Comp_k10'; %%% with DFile = 'Dcircls_49_a-57';
% RunFile = 'Biaxial_Comp_k12'; %%% with DFile = 'Dcircls_49_a-57';
% RunFile = 'Biaxial_Comp_k14'; %%% with DFile = 'Dcircls_49_a-57';
% RunFile = 'Biaxial_Comp_k15'; %%% with DFile = 'Dcircls_49_a-57';
  RunFile = 'Biaxial_Comp_k16'; %%% with DFile = 'Dcircls_49_a-57';
% RunFile = 'Biaxial_Comp_k18'; %%% with DFile = 'Dcircls_49_a-57';
% RunFile = 'Biaxial_Comp_k20'; %%% with DFile = 'Dcircls_49_a-57';
%
  DFile = 'Dcircls_64_a-1';
  DFile = 'Dcircls_49_a-57';
% DFile = 'Dcircls_256_a-100';
%
% output file path
  File_Name = cstrcat('Output_', RunFile, '_', DFile);
  File_Path = cstrcat(Output_Directory, File_Name, '.mat');
%
  load('-binary', File_Path, 'Headings', 'Output_Cell_Array');
%
  if Get_R0_alt
    File_Name_alt = cstrcat('Output_alt_', RunFile, '_', DFile);
    File_Path_alt = cstrcat(Output_Directory, File_Name_alt, '.mat');
%
    R0_alt_FilePath = cstrcat(Output_Directory, File_Name, '_R0_alt.mat');
    is_R0_alt_File = isfile(R0_alt_FilePath);
%
    load('-binary', File_Path_alt, 'Headings_alt', 'Output_Cell_Array_alt');
%
    if is_R0_alt_File
      load('-binary', R0_alt_FilePath, ...
           'NeutralVec_alt', 'NeutralGrp_alt', 'NeutralRank_alt');
    end
  end
%
  idim = 2;
%
  Headings;
%
  Initial_Time = 5;
  End_Time = size(Output_Cell_Array,1);
%
  Range = [Initial_Time:End_Time];
  Length = length(Range);
%
  Time_ = zeros(Length,1);
%
  Strain_11 = zeros(Length,1);
  Strain_22 = zeros(Length,1);
  Stress_11 = zeros(Length,1);
  Stress_22 = zeros(Length,1);
  q = zeros(Length,1);
  p = zeros(Length,1);
%
  M_lambda_all = zeros(Length,1);
  RCond_BD_all = zeros(Length,1);
%
  M_contacts = zeros(Length,1);
  M_sliding = zeros(Length,1);
%
  M_isPD = zeros(Length,1);
  M_isPmatrix = zeros(Length,1);
  M_isLi_matrix = zeros(Length,1);
  n_soln = zeros(Length,1);
  Method_used_all = zeros(Length,1);
%
  max_lambda = zeros(Length,1);
  Path_Stability = zeros(Length,1);
  W_2 = zeros(Length,1);
  B_2 = zeros(Length,1);
  E_2 = zeros(Length,1);
%
  M_isZmatrix = zeros(Length,1);
  M_isRule5 = zeros(Length,1);
  M_isSymmetric = zeros(Length,1);
%
% Get_Pmatrix_list = 1;
% Pmatrix_FilePath = cstrcat(Output_Directory, File_Name, '_Pmatrix.mat');
% if Get_Pmatrix_list; is_Pmatrix_File = isfile(Pmatrix_FilePath); end
%
  addpath('~/LCP/Libraries/Murty');
%
  AltSolutions_FilePath = cstrcat(Output_Directory, File_Name, '_LCP.mat');
  if Get_AltSolutions
    is_AltSolutions_File = isfile(AltSolutions_FilePath);
%
    if ~is_AltSolutions_File
%     the mat-file of results does not yet exist.  Initialize the vectors
%     and cell array that will hold this information
%
      M_all = cell(1,Length);
      q_all = cell(1,Length);
      Lambda_all = cell(1,Length);
      Lambda_r_all = cell(1,Length);
      DegenSoln_all = cell(1,Length);
      IsolatedSoln_all = cell(1,Length);
      NonIsolated_all = cell(1,Length);
      NonIsolatedGrp_all = cell(1,Length);
      ColCompetent_all = cell(1,Length);
      DegenMatrix_all = cell(1,Length);
      n_DegenMatrix_all = zeros(1,Length);
    end
  end
%
  M_alt_all = cell(1,Length);
%
  for i = 1:Length
    j = Range(i);
%
    Time_(i) = Output_Cell_Array{j,1};
%
    strain_ = Output_Cell_Array{j,3};
    Strain_11(i) = strain_(1,1);
    Strain_22(i) = strain_(2,2);
%
    def_gradient = strain_ + eye(3);
%
    strain_prev = Output_Cell_Array{j-1,3};
    def_gradient_prev = strain_prev + eye(3);
    ddef = (def_gradient - def_gradient_prev) * inv(def_gradient);
%
    dV(i) = prod(diag(def_gradient)) - 1;
%
    stress_ = Output_Cell_Array{j,4};
    Stress_11(i) = stress_(1,1);
    Stress_22(i) = stress_(2,2);
%
    stress_prev = Output_Cell_Array{j-1,4};
    dstress = stress_ - stress_prev;
%
    q = Stress_11 - Stress_22;
    p = 0.5 * (Stress_11 + Stress_22);
%
    M_lambda_all(i) = Output_Cell_Array{j,8};
    RCond_BD_all(i) = Output_Cell_Array{j,13};
%
    if ~isempty(Output_Cell_Array{j,14})
      Method_used_all(i) = Output_Cell_Array{j,14};
    else
      Method_used_all(i) = 0;
    end
    is_Soln_all(i) = Output_Cell_Array{j,15};
%
%   an Index_Type=2 means that the contact is sliding
    M_sliding(i) = sum(Output_Cell_Array{j,12}==2);
%
%   number of contacts
    M_contacts(i) = Output_Cell_Array{j,5};
%
%   whether matrix M_hat is positive definite
    M = Output_Cell_Array{j,9};
    if ~isempty(M)
%     M_isPD(i) = isdefinite(0.5*(M+M'),100*eps);
      M_isPD(i) = isdefinite(0.5*(M+M'));
    else
      M_isPD(i) = 1;
    end
%
%   test whether M is symmetric
    M_isSymmetric(i) = issymmetric(M, 100*eps);
%
%   test whether M is a Z-matrix
    M_isZmatrix(i) = all(all((M - diag(diag(M)))<=100*eps));
%
%   test of Rule 5
    M_isRule5(i) = all(all(M>=0)) & all(diag(M))>100*eps;
%
%   number of solutions, and largest lambda values
    if isempty(Output_Cell_Array{j,11})
      n_soln(i) = -1;
      max_lambda(i) = 0;
    else
      n_soln(i) = size(Output_Cell_Array{j,16}, 2);
      max_lambda(i) = max(Output_Cell_Array{j,11});
    end
%
    xcell = Output_Cell_Array{j,20};
    Vol = prod(diag(xcell)(1:idim));
%
%
    Path_Stability(i) = Output_Cell_Array{j,19};
    W_2(i) = (1/Vol)*Output_Cell_Array{j,19};
    B_2(i) = sum(sum(ddef .* dstress));
    E_2(i) = B_2(i) - W_2(i);
%
%---the above information was extracted from the original simulaiton output
%   file.  Now, we re-solve the LCPs and do other analyses based on the
%   origional M and q data. The new results are then stored
%
%   solve the LCP again
    if Get_AltSolutions && ~is_AltSolutions_File
      if mod(i, 100)==0 % || (i>4700 && i<4800)
        disp(i)
      end
%
      M_hat = M;
      q_hat = Output_Cell_Array{j,10};
%
      if i==7944
        Lambda_all{i} = [];
        Lambda_r_all{i} = [];
        DegenSoln_all{i} = [];
        IsolatedSoln_all{i} = [];
        NonIsolated_all{i} = [];
        NonIsolatedGrp_all{i} = [];
        ColCompetent_all{i} = [];
        DegenMatrix_all{i} = 0;
        n_DegenMatrix_all(i) = 0;
      elseif ~isempty(M)
        [ Lambda, Lambda_r, DegenSoln, IsolatedSoln, ...
          NonIsolated, NonIsolatedGrp, ...
          ColCompetent, DegenMatrix, n_DegenMatrix] = ...
        LCP_enumerate4(M_hat, q_hat, 1e-9, 1000*eps); % 10000
%
%       M_lambda_all(i);
        M_all{i} = M;
        q_all{i} = q_hat;
        Lambda_all{i} = Lambda;
        Lambda_r_all{i} = Lambda_r;
        DegenSoln_all{i} = DegenSoln;
        IsolatedSoln_all{i} = IsolatedSoln;
        NonIsolated_all{i} = NonIsolated;
        NonIsolatedGrp_all{i} = NonIsolatedGrp;
        ColCompetent_all{i} = ColCompetent;
        DegenMatrix_all{i} = DegenMatrix;
        n_DegenMatrix_all(i) = n_DegenMatrix;
      else
%       M_lambda_all(i);
        M_all{i} = [];
        q_all{i} = [];
        Lambda_all{i} = [];
        Lambda_r_all{i} = [];
        DegenSoln_all{i} = [];
        IsolatedSoln_all{i} = [];
        NonIsolated_all{i} = [];
        NonIsolatedGrp_all{i} = [];
        ColCompetent_all{i} = [];
        DegenMatrix_all{i} = 0;
        n_DegenMatrix_all(i) = 0;
      end
    end
%
    if Get_R0_alt
      M_alt = Output_Cell_Array_alt{j,3};
      if i==7944
      elseif ~isempty(M_alt)
        M_alt_all{i} = M_alt;
      else
        M_alt_all{i} = [];
      end
    end
%
  end
%
%
  Pmatrix_FilePath = cstrcat(Output_Directory, File_Name, '_Pmatrix.mat');
  if Get_Pmatrix_list; is_Pmatrix_File = isfile(Pmatrix_FilePath); end
  if Get_Pmatrix_list && ~is_Pmatrix_File
%   Boolean list of whether M is a P-matrix 7944
    for i = 1:Length
      j = Range(i);
      if mod(i,100)==0
        disp(i)
      end
      if i==7944
        M_isPmatrix(i) = 1;
      else
        M = Output_Cell_Array{j,9};
        if Get_Pmatrix_list && ~is_Pmatrix_File
          M_isPmatrix(i) = isempty(M) || isPmatrix(M);
        end
      end
    end
    save('-binary', Pmatrix_FilePath, 'M_isPmatrix');
  elseif Get_Pmatrix_list
    load('-binary', Pmatrix_FilePath, 'M_isPmatrix');
  end
%
%
  Li_matrix_FilePath = cstrcat(Output_Directory, File_Name, '_Li_matrix.mat');
  if Get_Li_matrix_list; is_Li_matrix_File = isfile(Li_matrix_FilePath); end
  if Get_Li_matrix_list && ~is_Li_matrix_File
%   Boolean list of whether M is a Li-matrix:
%   Cui-Xia Li & Shi-Liang Wu, "A note on the unique solution of linear
%   complementarity problem." Cogent Mathematics, 3:1, 1271268.
    for i = 1:Length
      j = Range(i);
      if mod(i,100)==0
        disp(i)
      end
%
      M = Output_Cell_Array{j,9};
      if Get_Li_matrix_list && ~is_Li_matrix_File
        if isempty(M)
          M_isLi_matrix(i) = 1;
        else
          M_isLi_matrix(i) = ...
            max(svd(abs(eye(size(M)) - M))) ...
            < min(svd(eye(size(M)) + M));
        end
      end
    end
    save('-binary', Li_matrix_FilePath, 'M_isLi_matrix');
  elseif Get_Li_matrix_list
    load('-binary', Li_matrix_FilePath, 'M_isLi_matrix');
  end
%
  if Get_AltSolutions && ~is_AltSolutions_File
    save('-binary', AltSolutions_FilePath, ...
         'M_all', 'q_all', 'Lambda_all', 'Lambda_r_all', ...
         'DegenSoln_all', 'IsolatedSoln_all', 'NonIsolated_all', ...
         'NonIsolatedGrp_all', 'ColCompetent_all', ...
         'DegenMatrix_all', 'n_DegenMatrix_all');
  elseif Get_AltSolutions
    load('-binary', AltSolutions_FilePath, ...
         'M_all', 'q_all', 'Lambda_all', 'Lambda_r_all', ...
         'DegenSoln_all', 'IsolatedSoln_all', 'NonIsolated_all', ...
         'NonIsolatedGrp_all', 'ColCompetent_all', ...
         'DegenMatrix_all', 'n_DegenMatrix_all');
%
    M_isDegen = zeros(Length, 1);
    nIsolated_Solns = zeros(Length, 1);
    nNonIsolated_Solns = zeros(Length, 1);
    Incongruous = zeros(Length, 1);
    Path_Instability = zeros(Length, 1);
    nDegen = zeros(Length, 1);
%
    for i = 1:Length
      M_isDegen(i) = n_DegenMatrix_all(i)>0;
      nDegen(i) = size(n_DegenMatrix_all(i), 2);
%
      M_hat = M_all{i};
      q_hat = q_all{i};
      Lambda = Lambda_all{i};
%
      Limit1 = 1e-9;
      Index_Type = zeros(size(Lambda));
      for k = 1:size(Lambda,2)
        lambda = Lambda(:,k);
        w = M_hat*lambda + q_hat;
        Index_Type(find(lambda < Limit1 & w > Limit1),k) = 1;
        Index_Type(find(lambda > Limit1 & w < Limit1),k) = 2;
        Index_Type(find(lambda < Limit1 & w < Limit1),k) = 3; % degen. soln.
      end
%
%     first, find all of the non-isolated solutions.  These should occur
%     in groups intercepts, each with the same NonIsolatedGrp number
      IsolatedSoln = IsolatedSoln_all{i};
      NonIsolatedGrp = NonIsolatedGrp_all{i};
%
%     a scratch list of the solutions
      Scratch = 1:size(Lambda, 2);
%
%     find duplicates of Index_Type and eliminate any IsolatedGrp solutions
%     among duplicates
      [Unique, Unique_I, Unique_J] = unique(Index_Type', "rows");
%
      for kGrp = setdiff(unique(NonIsolatedGrp), 0)
        N = dec2bin(kGrp-1, length(q_hat));
        Ones = N=='1';
        J = find(Ones);
        Ar = eye(length(q_hat));
        Ar(:,J) = -M_hat(:,J);
        Max_S = max(abs(svd(M_hat)));
        rankAr = rank(Ar, 10*1e-9*Max_S);
        dimNonIsol = length(q_hat) - rankAr;
%
        if dimNonIsol~=0 && sum(NonIsolatedGrp==kGrp) > dimNonIsol
%         there are a sufficient number of intercepts for this solution
%         to qualify as a non-isolated solution
          nNonIsolated_Solns(i) = length(setdiff(unique(NonIsolatedGrp),0));
%
%         remove these solutions from the solution list, including any
%         duplicate solutions that are masquerading as isolated solutions
          for l = find(NonIsolatedGrp==kGrp)
            for m = setdiff(Scratch, l)
              if all(Index_Type(:,m)==Index_Type(:,l))
                Scratch = setdiff(Scratch, m);
              end
            end
            Scratch = setdiff(Scratch, l);
          end
        elseif dimNonIsol~=0 && sum(NonIsolatedGrp==kGrp) == dimNonIsol ...
               && size(Lambda,2)>2
          l = find(NonIsolatedGrp>0);
          rmScratch = [];
          for m = 1:length(Scratch)
            if Index_Type(find(Index_Type(:,m)~=Index_Type(:,l)),m)==3
              rmScratch = [rmScratch, m];
            end
          end
          nNonIsolated_Solns(i) = 1;
          rmScratch = [rmScratch, l];
          Scratch = setdiff(Scratch, rmScratch);
        elseif dimNonIsol~=0
%         this solution is a bogus non-isolated solution.  Remove it.
          Scratch = setdiff(Scratch, find(NonIsolatedGrp==kGrp));
        end
      end
%
      if ismember(i, ...
          [5320 5368, 5380 5402, 5512, 5612, 5684, 5694, 5704, 5714, 5724])
        nNonIsolated_Solns(i) = 1;
        Scratch = [];
      end
%
      X = Index_Type(:,Scratch);
      X = unique(X', "rows")';
      nIsolated_Solns(i) = size(X, 2);
%
%     no solutions
      Incongruous(i) = (nIsolated_Solns(i) + nNonIsolated_Solns(i)) == 0 ...
                       & length(q_hat) ~= 0;
%
%     path instability
      Path_Instability(i) = (nIsolated_Solns(i) + nNonIsolated_Solns(i)) > 1;
    end
  end
%
  Smatrix_FilePath = cstrcat(Output_Directory, File_Name, '_Smatrix.mat');
  if Get_Smatrix_list; is_Smatrix_File = isfile(Smatrix_FilePath); end
  if Get_Smatrix_list && Get_Pmatrix_list && Get_AltSolutions ...
                      && ~is_Smatrix_File
%   read some functions that will be used in the optimization function "sqp"
%   sqp_functions_S;
    M_isSmatrix = zeros(1,Length);
%
    for i = 1:Length
      if mod(i, 100)==0
        disp(i)
      end
%
      M_hat = M_all{i};
%
%     if M_isPmatrix || isempty(M_hat);
      if isempty(M_hat);
        M_isSmatrix(i) = 1;
      elseif size(M_hat, 1)==1
        M_isSmatrix(i) = M_hat > 0;
      else
        Zeros = zeros(size(M_hat,1),1);
        Ones = ones(size(Zeros));
%
        warning off
        try
        [X, OBJ, INFO, ITER, NF, LAMBDA] = ...
          sqp(Ones, @Obj_S, [], @Con_S, Zeros, []);
%         sqp(rand(size(Ones)), @Obj_S, [], @Con_S, Zeros, []);
        catch
        [X, OBJ, INFO, ITER, NF, LAMBDA] = ...
          sqp(rand(size(Ones)), @Obj_S, [], @Con_S, Zeros, []);
        end_try_catch
        warning on
%
        M_isSmatrix(i) = -OBJ > 10*eps;
      end
    end
%
    save('-binary', Smatrix_FilePath, 'M_isSmatrix');
  elseif Get_Smatrix_list && is_Smatrix_File
    load('-binary', Smatrix_FilePath, 'M_isSmatrix');
  end
%
  Sensitive_FilePath = cstrcat(Output_Directory, File_Name, '_Sensitive.mat');
  if Get_Sensitive; is_Sensitive_File = isfile(Sensitive_FilePath); end
  if Get_Sensitive && Get_AltSolutions && ~is_Sensitive_File
    Sensitive_all = cell(1,Length);
    for i = 1:Length
      if mod(i, 101)==0
        disp(i)
      end
%
      [Sensitive_all{i}] = ...
        LCP_Sensitive(M_all{i}, q_all{i}, Lambda_all{i}, 1000*eps);
    end
    save('-binary', Sensitive_FilePath, 'Sensitive_all');
  elseif Get_Sensitive
    load('-binary', Sensitive_FilePath, 'Sensitive_all');
%
%   whether any of the solutions are sensitive
    Sensitive_any = zeros(Length,1);
    for i = 1:Length
      Sensistive_list = Sensitive_all{i};
      Sensitive_any(i) = any(Sensistive_list);
    end
  end
%
%
  R0_FilePath = cstrcat(Output_Directory, File_Name, '_R0.mat');
  if Get_R0; is_R0_File = isfile(R0_FilePath); end
  if Get_R0 && Get_AltSolutions && ~is_R0_File
    NeutralVec = cell(1,Length);
    NeutralGrp = cell(1,Length);
    NeutralRank = cell(1,Length);
%
    for i = 1:Length
      if mod(i, 100)==0
        disp(i)
      end
%
      if i == 7944
        NeutralVec{i} = [];
        NeutralGrp{i} = [];
        NeutralRank{i} = [];
      else
        [NeutralVec{i}, NeutralGrp{i}, NeutralRank{i}] = ...
           LCP_Neutral2(M_all{i}, 1e-9, 1000*eps);
      end
    end
    save('-binary', R0_FilePath, 'NeutralVec', 'NeutralGrp', 'NeutralRank');
  elseif Get_R0
    load('-binary', R0_FilePath, 'NeutralVec', 'NeutralGrp', 'NeutralRank');
%
%   whether there are any modes of neutral equilibrium
    R0_modes = zeros(Length,1);
    for i = 1:Length
      NeutralGrp_list = NeutralGrp{i};
      R0_modes(i) = length(NeutralGrp_list);
    end
  end
%
%
% Get_R0_alt = 0;
  if Get_R0_alt && ~is_R0_alt_File
    NeutralVec_alt = cell(1,Length);
    NeutralGrp_alt = cell(1,Length);
    NeutralRank_alt = cell(1,Length);
%
    for i = 1:Length
      if mod(i, 100)==0
        disp(i)
      end
%
      if i == 7944
        NeutralVec_alt{i} = [];
        NeutralGrp_alt{i} = [];
        NeutralRank_alt{i} = [];
      else
        [NeutralVec_alt{i}, NeutralGrp_alt{i}, NeutralRank_alt{i}] = ...
           LCP_Neutral2(M_alt_all{i}, 1e-9, 1000*eps);
      end
    end
    save('-binary', R0_alt_FilePath, ...
         'NeutralVec_alt', 'NeutralGrp_alt', 'NeutralRank_alt');
  elseif Get_R0_alt
%
%   whether there are any modes of neutral equilibrium
    R0_alt_modes = zeros(Length,1);
    for i = 1:Length
      Neutral_alt_Grp_list = NeutralGrp_alt{i};
      R0_alt_modes(i) = length(Neutral_alt_Grp_list);
    end
  end
%
  R_FilePath = cstrcat(Output_Directory, File_Name, '_R.mat');
  if Get_R; is_R_File = isfile(R_FilePath); end
  if Get_R && ~is_R_File
    M_isRmatrix = zeros(Length,1);
%
    for i = 1:Length
      j = Range(i);
      if mod(i,100)==0
        disp(i)
      end
%
      M = Output_Cell_Array{j,9};
%
      if i ~= 7944
        if ~isempty(M) 
          [ Lambda_, Lambda_r_, DegenSoln_, IsolatedSoln_, ...
            NonIsolated_, NonIsolatedGrp_, ...
            ColCompetent_, DegenMatrix_, n_DegenMatrix_] = ...
          LCP_enumerate4(M, 1e-5*ones(size(M,1),1), 1e-9, 1000*eps); % 10000
%
          if isempty(Lambda_) || all(all(Lambda_ < 1000*eps))
            M_isRmatrix(i) = 1;
          end
        else
          M_isRmatrix(i) = 1;
        end
      else
        M_isRmatrix(i) = 1;
      end
    end
    save('-binary', R_FilePath, 'M_isRmatrix');
  elseif Get_R
    load('-binary', R_FilePath, 'M_isRmatrix');
  end
%
  Instab_FilePath = cstrcat(Output_Directory, File_Name, '_Instab.mat');
  if Get_Instab; is_Instab_File = isfile(Instab_FilePath); end
  if Get_Instab && is_Instab_File
    load('-binary', Instab_FilePath, 'Headings_Instab', 'Output_Instab_Array');
%
    nLambda_Instab = zeros(Length,1);
%
    for i = 1:Length   %size(Output_Instab_Array,3)
      j = Range(i);
      if ~isempty(Output_Instab_Array{j,3})
        nLambda_Instab(i) = length(Output_Instab_Array{j,3});
      end
    end
  end
%
% correct the various matrix classes for round-off error
  if Get_AltSolutions && Get_Pmatrix_list
    M_isPmatrix(find(M_isPmatrix & M_isDegen)) = 0;
    M_isPD(find(M_isPD & ~M_isPmatrix)) = 0;
  end
%
% read the B-file from the DEM simulation
  if 1
    BfileDirectory = './';
    B_file_name = 'BBiaxial_Comp_k16_oval';
%
    [ Strain_11_DEM, q_DEM, p_DEM] = ...
    DEM_Read_Bfile(BfileDirectory, B_file_name);
  end
%
% create a plotting file
  if Create_Plot_File
    Path_Plot_Gen_File = cstrcat(Output_Directory, File_Name, '_Plot_Gen.txt'); 
    File_Plot_Gen_No = fopen(Path_Plot_Gen_File, 'w');
    fprintf(File_Plot_Gen_No, ...
            cstrcat('#  1) iOutput\n', ...
                    '#  2) Time\n', ...
                    '#  3) -Eps_11\n', ...
                    '#  4) q/p\n', ...
                    '#  5) dv\n', ...
                    '#  6) M_tacts\n', ...
                    '#  7) M_active\n', ...
                    '#  8) M_sliding\n', ...
                    '#  9) M_isSymmetric\n', ...
                    '# 10) M_isPD\n', ...
                    '# 11) M_isPmatrix\n', ...
                    '# 12) M_isZmatrix\n', ...
                    '# 13) M_isDegen\n', ...
                    '# 14) nIsolated_Solns\n', ...
                    '# 15) nNonIsolated_Solns\n', ...
                    '# 16) Path_Instability\n', ...
                    '# 17) Sensitive_any\n', ...
                    '# 18) R0_modes\n', ...
                    '# 19) Incongruous\n', ...
                    '# 20) nLambda_Instab\n', ...
                    '# 21) M_isRule5\n', ...
                    '# 22) M_isRmatrix\n'));
%
    M_active = M_lambda_all;
%
    for i = 1:Length
      j = Range(i);
      fprintf(File_Plot_Gen_No, ...
              cstrcat( ...
              '%5i  %5.1f  %10.5e  %10.5e  %10.5e ', ...
              '%3i  %2i  %2i  %1i  %1i  ', ...
              '%1i  %1i  %1i  %3i  %2i  ', ...
              '%2i  %1i  %2i  %1i  %4i  %1i %1i\n'), ...
              j, ...                     %  1
              Time_(i), ...              %  2
              -Strain_11(i), ...         %  3
              q(i)/p(i), ...             %  4
              dV(i), ...                 %  5
              0.5*M_contacts(i), ...     %  6
              M_active(i), ...           %  7
              M_sliding(i), ...          %  8
              M_isSymmetric(i), ...      %  9
              M_isPD(i), ...             % 10
              M_isPmatrix(i), ...        % 11
              M_isZmatrix(i), ...        % 12
              M_isDegen(i), ...          % 13
              nIsolated_Solns(i), ...    % 14
              nNonIsolated_Solns(i), ... % 15
              Path_Instability(i), ...   % 16
              Sensitive_any(i), ...      % 17
              R0_modes(i), ...           % 18
              Incongruous(i), ...        % 19
              nLambda_Instab(i), ...     % 20
              M_isRule5(i), ...          % 21
              M_isRmatrix(i));           % 22
    end
    fclose(File_Plot_Gen_No);
  else
    disp('NO PLOTTING FILE CREATED.');
  end
%
% I_2 values of the non-isolated solutions:
Q = ...
 [-1.7855e-09  -1.7855e-09;
  -1.9218e-09  -1.9218e-09;
  -1.9207e-09  -1.9207e-09;
  -1.9135e-09  -1.9135e-09;
  -1.9061e-09  -1.9060e-09;
  -1.4433e-09  -1.4433e-09;
  -1.4438e-09  -1.4438e-09;
  -1.4426e-09  -1.4426e-09;
  -1.4408e-09  -1.4408e-09;
  -1.4376e-09  -1.4376e-09;
  -1.4359e-09  -1.4359e-09;
  -1.4343e-09  -1.4343e-09;
  -1.4327e-09  -1.4327e-09;
  -1.4311e-09  -1.4311e-09;
  -1.4280e-09  -1.4280e-09;
  -1.4264e-09  -1.4264e-09];
%
  if 0
    iList = [];
    List = find(nIsolated_Solns>1);
%   List = find(nNonIsolated_Solns>0);
    for i = List'
      j = Range(i);
%
      M_hat = M_all{i};
      q_hat = q_all{i};
      Lambda = Lambda_all{i};
%
      Limit1 = 1e-9;
      Index_Type = zeros(size(Lambda));
      for k = 1:size(Lambda,2)
        lambda = Lambda(:,k);
        w = M_hat*lambda + q_hat;
        Index_Type(find(lambda < Limit1 & w > Limit1),k) = 1;
        Index_Type(find(lambda > Limit1 & w < Limit1),k) = 2;
        Index_Type(find(lambda < Limit1 & w < Limit1),k) = 3; % degen. soln.
      end
%
      disp(i)
      Index_Type;
      iList = [iList, i];
    end
  end
%
% Output_Cell_Array{50,12}
% [M_sliding RCond_BD_all<1e-17]
% plot(-Strain_11, q./p, '-')
% plot(Time_, q./p, '-')
% plot(-Strain_11, dV);
%
