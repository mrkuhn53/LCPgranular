%
% |+|+|+|+|+|+|+|+|+|  Revision 3   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
function [ x, X, fval, ...
           DegenSoln, IsolatedSoln, NonIsolated, ...
           ColCompetent, DegenMatrix, n_DegenMatrix] = ...
         LCP3_Solve_Problem_alt( ...
           M_hat, q_hat, x0, Method, PHI, Hc, System, LibrariesPath)
%
% this function solves the LCP that is associated with Problem 2
%      M*x + q >= 0, x >= 0, x'*(M*x) = 0
%
% with Octave, we can choose among several solvers. With Matlab, we use the
% quadprog() function. The solvers include ones that solve the
% quadratic program of Problem 3, and return a single solution, even when
% there are multiple solutions (this type of solver is the only option
% with Matlab).  With Octave, the solver "LCP_enumerate.m" is available,
% which will find all solutions of the LCP, using an enumerative approach
% (this method does not scale well with the problem size, and computation
% time can become intractable when the size of M is greater than 15)
%
% -------------- INPUT ------------------------
%  System = 1, for Octave; 2 for Matlab
%  Method =    with Octave: 1, 2, 3, 4.  With Matlab: 1.  If Method==0,
%              the Method 1 will be used
%  LibrariesPath = character string.  Some of the methods use functions
%                  that can be located in a folder of libraries (with Octave,
%                  Methods 3 and 4).  If all of the functions are located
%                  in the current folder (pwd), then input './'
%  x0 =            initial guess of solution "x".  If empty "[]", the
%                  zeros are used
%
% -------------- OUTPUT ------------------------
% x =         solution
% X =         solutions.  Can find multiple solutions with Octave and  
%             Method = 4
% fval =      value of the objective function. Can be empty with some methods
% DegenSoln, IsolatedSoln, NonIsolated, ColCompetent, DegenMatrix =
%             These are output from the LCP_enumerate() function, which
%             returns information, such as whether the matrix and solutions
%             are degenerate, whether solutions are isolated, etc.
%             See the function LCP_enumerate() for descriptions
%
% this code works with both Octave and Matlab.  Create the following
% two constants. The two softwares have subtle differences.
  Octave = 1; Matlab = 2;
%
if 0
% cell array of library folders
  Folder = cell(100,1);
%
% List of Octave Methods that can be used to solve the LCP, and the folder
% locations that contain the library
  qpOctave = 1; Folder{qpOctave} = [];         % Octave's qp() solver
  Quadprog = 2; Folder{Quadprog} = [];         % Octave's quadprog() solver
  COMPASS =  3; Folder{COMPASS} =  'COMPASS';  % COMPASS solver
  Murty =    4; Folder{Murty} =    'Murty';    % LCP_enumerate() solver
  Yuval =    5; Folder{Yuval} =    'Yuval';    % Mathworks code by Yuval
  Almqvist = 6; Folder{Almqvist} = 'Almqvist'; % Mathworks code, Andreas Almqvist
  CompEcon = 7; Folder{CompEcon} = 'CompEcon/CEtools';
  sqpOctave =8; Folder{sqpOctave} = [];        % Octave's sqp() solver
%  
  if Method==0 || isempty(Method)
%   default method for Octave and Matlab
    Method = 1;
  end
end
%
% problem size
  M_lambda = size(M_hat,1);
%
% M_hat = M;
% q_hat = q;
  x = []; lambda = [];
  X = []; Lambda = [];
  fval = [];
  DegenSoln = []; IsolatedSoln = []; NonIsolated = [];
  ColCompetent = []; DegenMatrix = [];
  n_DegenMatrix = 0;
%
% initial guest of the solution, as required with the qp(), quadprog(),
% and Compass() methods.  No initial guess is used with Murty's enumerative
% method.
  if isempty(x0)
    lambda0 = zeros(M_lambda,1);
  else
    lambda0 = x0;
  end
%
  if M_lambda~=0
%   therefore, there are active contacts and there is a possibility that 
%   the system is not entirely elastic, and there might be some non-zero 
%   values among the lambdas.  On the other hand, when M_lambda==0, 
%   the system is entirely elastic
%
    if System==Octave
%     with Octave, we have a choice of methods for solving the LCP:
%     1) solving the quadratic programming (QP) problem using the built-in
%        Octave function qp(). Note that this function requires an initial
%        guess "lambda0".  The function returns a single solution, even
%        when multiple solutions exist.  The value of lambda0 cna be
%        manually changed to search for other contacts.
%     2) solving the QP using the quadprog() function is part of the Octave's
%        "optim" package.  The package must be installed on your computer, 
%        along with the dependency packages: struct, statistics, and io.
%        Note that this function requires an initial guess "lambda0".  
%        The function returns a single solution, even when multiple 
%        solutions exist.  The value of lambda0 cna be manually changed to
%        search for other contacts.
%     3) solving the LCP using the COMPASS package of Stefan Schmelzer.
%        This COMPASS solver is a free clone of the robust PATH solver,
%        and COMPASS is licensed under the GNU General Public License 
%        as published by the Free Software Foundation, version 3 or
%        later.  I have modified the code, by changing the script file
%        COMPASSmain.m into a function.  The COMPASS package is among the
%        mode robust solvers that I have used.  However, the function 
%        returns a single solution, even when multiple solutions exist.
%        The COMPASS package is available here:
%          https://github.com/haraldschilly/compass-solver-matlab
%        Documenation of the COMPASS package is available here:
%          Stefan Schmelzer, “COMPASS: A Free Solver for Mixed 
%          Complementarity Problems", 2012.
%     4) an enumerative method that the author wrote, based upon the
%        algorithm described by Katta G. Murty and Feng-Tien Yu.  The
%        method finds all solutions of the LCP, determines whether
%        the M-matrix is degnerate, determines whether the solutions
%        are degenerate, and determines wither the solutions are isolated.
%        Because the method is enumerative, the solution time of order
%        2^n, so it is not very efficient.
%        The Murty/Yu method is described here, on page 8:
%          Katta G. Murty and Feng-Tien Yu, Linear Complementarity, Linear 
%          and Nonlinear Programming. Internet Edition, 1997(?).
%
%     cell array of library folders
      Folder = cell(100,1);
%
%     Path to the directory that holds all of the libraries
%     LibrariesPath = '~/MyDrive/Kuhn_in_Progress/LCP/Libraries/';
%
%     List of Octave libraries that can be used to solve the LCP
      qpOctave = 1; Folder{qpOctave} = [];
      Quadprog = 2; Folder{Quadprog} = [];
      COMPASS =  3; Folder{COMPASS} =  'COMPASS';
      Murty =    4; Folder{Murty} =    'Murty';
      Yuval =    5; Folder{Yuval} =    'Yuval';    % Mathworks code by Yuval
      Almqvist = 6; Folder{Almqvist} = 'Almqvist'; % Andreas Almqvist
      CompEcon = 7; Folder{CompEcon} = 'CompEcon/CEtools';
      sqpOctave =8; Folder{sqpOctave} = [];
%
if 0
%     Select the method
%     Method = qpOctave; % Octave's built-in qp() function
%     Method = Quadprog; % the quadprog() function in package "optim"
%     Method = COMPASS;  % the COMPASS package of Stefan Schmelzer
%     Method = Murty;    % enumerative method described by Murty
%
%     add the library (Method) folder to Octave's path
      if ~isempty(Folder{Method})
        AddPath = cstrcat(LibrariesPath,Folder{Method});
%       if this library is not already included in the path(), then add it
        if isempty(strfind(path(),Folder{Method}))
          if isfolder(AddPath)
            addpath(AddPath,0);
          else
            disp( ...
              cstrcat('Ooops! Can not find the folder for this library.', ...
                      Folder(Method))),
            ERROR_1
          end
        end
      end
end
%
      switch(Method)
        case qpOctave
%         the qp() function is built into standard Octave
          [lambda, fval, info, LMs] = ...
          qp(lambda0, ...               % initial guess
             2*M_hat, q_hat, ...        % min 0.5*lambda*M_hat*lambda
             ...                        %     + q_hat'*lambda
             [], [], ...                % no equality constraints
             zeros(M_lambda,1), [], ... % lambda >= 0, no upper bound
             -q_hat, M_hat, [] ...    % M_hat*lambda >= -q_hat, no upper bound
            ); 
        case Quadprog
%         the quadprog() function is part of the "optim" package.  The package 
%         must be installed on your computer, along with the dependency 
%         packages: struct, statistics, and io.  The following command loads 
%         the optim package, so that the quadprog() program can be used
          pkg load optim
%
          [lambda, fval, exitflag, output, LMs] = ...
          quadprog(2*M_hat, q_hat, ...        % min 0.5*lambda*M_hat*lambda 
                   ...                        %     + q_hat'*lambda
                   -M_hat, q_hat, ...         % inequaltity constraints, but
                   ...                        % with upper bound
                   [], [], ...                % no equality constraints
                   zeros(M_lambda,1), [], ... % lambda >= 0, no upper bound`
                   lambda0);                  % initial guess
%
        case COMPASS
          Mbar = M_hat;
          dQbar = q_hat;
          x0 = lambda0;
          n = length(dQbar);
%
          max_major=15*n;        % restart the algorithm after max_majoriter
                                 % major iterations.
          max_minor=15*n;        % stop the pivot procedure after
                                 % max_minoriter minor iterations.
          max_restarts=50;       % stop COMPASS after min12 to max60 restarts
                                 % without a solution.
          target=1.e-8;          % quit when a point with merit function value
                                 % <= target is found
          display_depth = 0;     % display details on solving process:
                                 % 0 -solution only
                                 % 1 -major iteration info
                                 % 2 -major and minor iteration info
%
          [lambda] = COMPASS_func(Mbar, dQbar, x0', ... % NOTE: x0' (transpose)
                                  max_major, max_minor, max_restarts, target);
          lambda = lambda';
%
        case Murty
%         [Lambda, DegenSoln, IsolatedSoln, NonIsolated, ...
%           ColCompetent, DegenMatrix] = ...
%         LCP_enumerate(M_hat,q_hat);
%
          [ Lambda, Lambda_r, DegenSoln, IsolatedSoln, NonIsolated, ...
            NonIsolatedGrp, ColCompetent, DegenMatrix, n_DegenMatrix] = ...
          LCP_enumerate4(M_hat, q_hat, 1e-8, 1000*eps);
%
          if ~isempty(Lambda)
            lambda = Lambda(:,1);
          else
            lambda = [];
          end
%
        case Yuval
%         using the code of Yuval:
%           Yuval (2022). LCP / MCP solver (Newton-based) 
%          (https://www.mathworks.com/matlabcentral/fileexchange/
%           20952-lcp-mcp-solver-newton-based), 
%          MATLAB Central File Exchange. Retrieved October 15, 2022.
%
          Mbar = M_hat;
          dQbar = q_hat;
%
          l = 0;
          u = 1e50;
          display = 0;
%
%         lambda = LCP(Mbar, dQbar, l, u, x0, display);
          lambda = LCP(Mbar, dQbar);
%
        case Almqvist
%         the code of Andreas Almqvist, available here:
%           Andreas Almqvist (2022). A pivoting algorithm solving linear
%           complementarity problems 
%           (https://www.mathworks.com/matlabcentral/fileexchange/41485-a-
%              pivoting-algorithm-solving-linear-complementarity-problems), 
%           MATLAB Central File Exchange. Retrieved October 26, 2022.
%
%         the code does not seem to allow for an initial guess
%
%         note that the code allows a 3rd and 4th inputs for tolerance and
%         number of iterations
%
          Mbar = M_hat;
          dQbar = q_hat;
%
          [w,lambda,retcode] = LCPSolve(Mbar, dQbar);
%
        case CompEcon
%         Paul L. Fackler & Mario J. Miranda
%
%
          Mbar = M_hat;
          dQbar = q_hat;
          a = zeros(size(dQbar));
          b = 1e20*ones(size(dQbar));
%
%         [lambda,z] = lcpbaard(Mbar, dQbar, a, b, lambda0);
%         [lambda,z] = lcplemke(Mbar, dQbar, a, b, lambda0);
          [lambda,z] = lcpsolve(Mbar, dQbar, a, b, lambda0);
%
        case sqpOctave
%         built-in sqp() function in Octave for sequential quadratic programming
%         Note:  the "PHI" and "Hc" objects are cell arrays that point to
%         functions that contain "global M_hat q_hat" statement, so that 
%         M_hat and q_hat are not seen in the argument list 
%         statements
          warning('off', 'all');
          [XX, OBJ, INFO, ITER, NF, LAMBDA] ...
            = sqp(lambda0, PHI, [], Hc, 0, []);
%         warning('on', 'all');
%
          lambda = XX;
          fval = OBJ;
%
        otherwise
          'The selected library does not exist', ERROR_2
        endswitch
%
    elseif System==Matlab
      lambda0 = zeros(M_lambda,1);
      [lambda, fval, exitflag, output, LMs] = ...
        quadprog(2*M_hat, q_hat, ...        % min 0.5*lambda*M_hat*lambda 
                 ...                        %     + q_hat'*lambda
                 -M_hat, q_hat, ...         % inequaltity constraints, but
                 ...                        % with upper bound
                 [], [], ...                % no equality constraints
                 zeros(M_lambda,1), [], ... % lambda >= 0, no upper bound`
                 lambda0);                  % initial guess
    end
  else
%   M_lambda==0, so the system is necessarily elastic.  We must use the
%   following for of an empty matrix "[]" so that the later multiplications
%   are legitimate
    lambda = zeros(M_lambda,1);
  end
%
% output values
  x = lambda;
  if isempty(Lambda)
    X = x;
  else
    X = Lambda;
  end
