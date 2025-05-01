function [ Output_Cell_Array, Headings, Output_Instab_Array, Headings_Instab,...
           init_output, istepo, iconto, ...
           last_Save_Time, last_Write_Time, iOutput, File_Path, ...
           File_Path_Instab, ...
           File_Path_alt, Output_Cell_Array_alt, Headings_alt] = ...
         Output_alt( ...
           Output_Cell_Array, Headings, Output_Instab_Array, Headings_Instab, ...
           Time, istep, istepo, ...
           icont, iconto, idim, ipts, ipts2, ...
           xcell, xcell_init, V_pq, r_pq_p, f_pq, RCond_BD, u, du, ...
           M_lambda, M_hat, q_hat, lambda, last_Save_Time, last_Write_Time, ...
           Method_used, is_Soln, Lambda, DegenSoln, DegenMatrix, ...
           Path_Stability, if_Inst_Output, ...
           P_L, H_x_x, H_x_lambda, H_lambda_x, H_lambda_lambda, ...
           C_wo_dxcell, H_x, H_lambda, Q_x, Q_lambda, ...
           RunFile, DFile, Output_Directory, init_output, iOutput, ...
           File_Path, File_Path_Instab, ...
           if_Output_File_alt, File_Path_alt, Output_Cell_Array_alt, ...
           Headings_alt, M_alt, q_alt, ...
           Warning, Time_to_Stop);
%
% this function creates a cell array that stores the results for each
% time-step
%
% number of cells in the cell array that holds output values for each time-step
  ncells = 20;
  ncells_Instab = 5;
%
% is this the first output?
  if init_output
%   output file path
    File_Name = cstrcat('Output_', RunFile, '_', DFile);
    File_Path = cstrcat(Output_Directory, File_Name, '.mat');
%
%   initialize the cell array
    Output_Cell_Array = cell(1, ncells);
%
%   names of cell columns
    Headings = cell(1,ncells);
%
%   if the RunFile gives if_Inst_Output = 1, then we search for loading
%   directions with negative 2nd-order work
    if if_Inst_Output
%     output file path
      File_Name_Instab = cstrcat('Output_', RunFile, '_', DFile, '_Instab');
      File_Path_Instab = cstrcat(Output_Directory, File_Name_Instab, '.mat');
%
%     initialize the cell array
      Output_Instab_Array = cell(1, ncells_Instab);
%
%     names of cell columns
      Headings_Instab = cell(1, ncells_Instab);
%
      Headings_Instab(:) = ...
        {'iOutput', ...       % 1
         'Lambda_Instab', ... % 2
         'NegEigenGrp', ...   % 3
         'NegEigenValue', ... % 4
         'NegEigenVector'};   % 5
    end
%
%   if_Output_File_alt = true, then we will create an additional file
%   that contains the alternative M_hat and q_hat files for a control program
%   with no displacement constraints.  During post-processing, this can 
%   be used for identifying neutral equilibrium under full force-control
    if if_Output_File_alt
      ncells_alt = 4;
      File_Name_alt = cstrcat('Output_alt_', RunFile, '_', DFile);
      File_Path_alt = cstrcat(Output_Directory, File_Name_alt, '.mat');
      Output_Cell_Array_alt = cell(1, ncells_alt);
      Headings_alt = cell(1,ncells_alt);
      Headings_alt = {'Time', ...
                      'iOutput', ...
                      'M_alt', ...
                      'q_alt'};
    end
%
    init_output = 0;
  end
%
% "Save" is whether to add the results to the cell array. "Write" is whether
% to do a periodic uploading of results to a file.
  Save = 1;
  if istep <= length(ipts)
%   Save = (Time - 1) >= last_Save_Time + ipts2(istep) ...
%          | istep ~= istepo;
    Write = (Time - 1) >= last_Write_Time + ipts(istep) ...
           | istep ~= istepo;
  else
    Save = 1;
    Write = 1;
  end
%
% save the istep value so that it can be reported in output files elsewhere
  if istep ~= istepo
    istep_out = istepo;
    istepo = istep;     % note that this value is returned as output
  else
    istep_out = istep;
  end
%
% save the icont value so that it can be reported in output files elsewhere
  if ~all(icont == iconto')
    icont_out = iconto;
    iconto = icont';     % not that this value is returned as output
  else
    icont_out = icont;
  end
%
% the Time-step for which the output applies.  Note that this value is not
% in this function's output list
  Time = Time - 1;                        % (1)
%
% the control-step for which the output applies. Note that this value is not
% in this function's output list
  istep = istep_out;                      % (2)
%
% displacement gradient, relative to t=0
  def_gradient = xcell * inv(xcell_init); % (3)
%
% current volume of assembly
  V = prod(diag(xcell)(1:idim));
%
% current stress. 
  stress = (1/V) * r_pq_p' * f_pq;        % (4)
%
% twice the number of contacts
  M = size(V_pq,1);                       % (5)
%
% type of control for the control-step
  icont = icont_out;                      % (6)
%
% last warning within this time-step
  LastWarn = lastwarn();                  % (7)
  lastwarn('');
%
% tests of M_bar, q_bar, and lambda
  M_lambda;                               % (8)
  M_hat;                                  % (9)
  q_hat;                                  % (10)
  lambda;                                 % (11)
%
% the types of complementarity condition for each element if lambda
  Limit1 = 1e-9;
  w = M_hat*lambda + q_hat;
  Index_Type = 1000*ones(M_lambda,1);
  Index_Type(find(lambda > Limit1 & w > Limit1)) = 0; % not a solution
  Index_Type(find(lambda < Limit1 & w > Limit1)) = 1;
  Index_Type(find(lambda > Limit1 & w < Limit1)) = 2;
  Index_Type(find(lambda < Limit1 & w < Limit1)) = 3; % degenerate solution
%
% the types of each of the M_lambda lambda's
  Index_Type;                             % (12)
%
% reciprocal condition number of the matrix used in the Bott-Duffin inverse
  RCond_BD;                               % (13)
%
  Method_used;                            % (14)
  is_Soln;                                % (15)
  Lambda;                                 % (16)
  DegenSoln;                              % (17)
  DegenMatrix;                            % (18)
  Path_Stability;                         % (19)
  xcell;                                  % (20)
%
% headings of the cell array
  Headings(:) = ...
    {'Time', ...
      'istep', ...
      'strain', ...
      'stress', ...
      'M', ...
      'icont', ...
     'LastWarn', ...
     'M_lambda', ...
     'M_hat', ...
     'q_hat', ...
     'lambda', ...
     'Index_Type', ...
     'RCond_BD', ...
     'Method_used', ...
     'is_Soln', ...
     'Lambda', ...
     'DegenSoln', ...
     'DegenMatrix', ...
     'Path_Stability', ...
     'xcell'};
%
  if Save
    last_Save_Time = Time;
    iOutput = iOutput + 1;
%
    Output_Cell_Array(iOutput,:) = ...
      {Time, istep, def_gradient-eye(3), stress, M, icont, ...
       LastWarn, M_lambda, M_hat, q_hat, lambda, Index_Type, RCond_BD, ...
       Method_used, is_Soln, Lambda, DegenSoln, DegenMatrix, ...
       Path_Stability, xcell};
%
    if if_Inst_Output
      if ~isempty(Q_lambda) ...
         && iOutput > 0 && iOutput~=7948 && M_lambda<=21
%
        if 0
%         Moore-Penrose inverse near KPD Eq. 68
          C_MP_inv_wo_dxcell = C_wo_dxcell' * inv(C_wo_dxcell*C_wo_dxcell');
%
%         the matrix P_Lperp in KPD Eq. 59
          P_Lperp_wo_dxcell = C_MP_inv_wo_dxcell * C_wo_dxcell;
%
%         the matrix P_L in KPD Eq. 58
          P_L_wo_dxcell = eye(size(C_wo_dxcell,2)) - P_Lperp_wo_dxcell;
%
          N = size(H_x,1) / 6;
          Rows2D = [ 1:3:3*N, 2:3:3*N, (3*N+3):3:6*N];
          Columns2D = Rows2D;
%
          [ Lambda_Instab, NegEigenGrp, NegEigenValue, NegEigenVector] = ...
          LCP_Instability2(P_L_wo_dxcell, ...
                           H_x(Rows2D,Columns2D), H_lambda(Rows2D,:), ...
                           Q_x(:,Columns2D), Q_lambda, 100*eps);
        elseif 1
          [ Lambda_Instab, NegEigenGrp, NegEigenValue, NegEigenVector] = ...
          LCP_Instability2(P_L, ...
                           H_x_x, H_x_lambda, ...
                           H_lambda_x, H_lambda_lambda, 100*eps);
        end
      else
        Lambda_Instab = [];
        NegEigenGrp = [];
        NegEigenValue = [];
        NegEigenVector = [];
      end
%
      Output_Instab_Array(iOutput,:) = ...
        {iOutput, Lambda_Instab, NegEigenGrp, NegEigenValue, NegEigenVector};
    end
%
    if if_Output_File_alt
      Output_Cell_Array_alt(iOutput,:) = {Time, iOutput, M_alt, q_alt};
    end
  end
%
  if Write || Warning || Time_to_Stop
    last_Write_Time = Time;
%
%   save -binary File_Path Output_Cell_Array
    save('-binary', File_Path, 'Headings', 'Output_Cell_Array');
%
    if if_Inst_Output
      save('-binary', File_Path_Instab, ...
           'Headings_Instab', 'Output_Instab_Array');
    end
%
    if if_Output_File_alt
      save('-binary', File_Path_alt, 'Headings_alt', 'Output_Cell_Array_alt');
    end
%
%   to retrieve the data:
%   load('-binary', File_Path, 'Headings', 'Output_Cell_Array');
  end
