function [NeutralVec, NeutralGrp, NeutralRank] = ...
         LCP_Neutral(M, EPS)
%
% solve the linear complementarity problem (LCP) using an enumerative
% approach that searches all 2^n possibilities to find all solutions
% of the LCP:
%   find "x" such that  M*x >= 0; x >=0; x'*(M*x) = 0
%
% Input: M = square matrix, nxn
%        q = column vector, nx1
%
% Output: X =            matrix with columns x that are solutions of the LCP
%                        Note that each column is unique, as redundant columns
%                        are removed.
%         DegenSoln =    row vector of 0s and 1s. A 0 means that the
%                        corresponding solution in X is not a degenerate
%                        solution. A 1 means that the corresponding solution 
%                        is degenerate.
%         IsolatedSoln = row vector of 0s and 1. A 1 means that the
%                        corresponding solution in X is an isolated
%                        solution.  A 2 means that the soltions is not
%                        isolated.
%         ColCompetent = 1 or 0, whether matrix M is column competent or not
%         DegenMatrix =  if the matrix is non-degenerate, then the empty
%                        matrix [] is assigned.  If M is degenerate, then
%                        DegenMatrix is a matrix with n rows.  
%                        The columns of DegenMatrix indicate the index
%                        sets for which the principal minor is zero.
%                        The index set is given by the rows of the
%                        1's in a column.  A matrix M can be degenerate
%                        with multiple prinicpal minors that are zero.
%                        Hence, possible multiple columns in DegenMatrix.
%         NonIsolated =  matrix of directions in which non-isolated solutions
%                        (like a polygonal edge) emanate from the solutions 
%                        in matrix X.  NonIsolated has n rows, and one column
%                        for each solution.  A column is filled with zeros, 
%                        if the solution is isolated.  Otherwise, the column
%                        contains a unit vector of the direction of the null
%                        space of the complemenary matrix, giving the
%                        direction "z" of solutions near to "x": x + eps*z.
%                        Note that the code only handles complemenary 
%                        matrices with a rank of n-1, not less.
%
% Limitation: for degenerate matrices, we currently are limited to cases
% in which complementary matrices are of rank n-1.
%
% Method: we use the method described my Murty and Yu, page 8:
%   Katta G. Murty and Feng-Tien Yu, Linear Complementarity: Linear and
%   Nonlinear Programming, Internet edition, 1997.

%   Murty, Katta G. "On the number of solutions to the complementarity 
%   problem and spanning properties of complementary cones." Linear 
%   Algebra and its Applications 5.1 (1972): 65-108.
%
% Limitation: the implementation time is of order 2^n, so large problems 
% (say, n>15) are intractable.  If you only need a single solution of a
% multi-solution LCP, or if you are sure that the LCP has only a single
% solution (for example, if M is a P-matrix), then there are several other
% Matlab/Octave codes that will solve the problem more efficiently than
% this function.  I suggest trying qp() or quadprog() in Octave, or 
% quadprog in Matlab.
%
  NeutralVec = [];
  NeutralGrp = [];
  NeutralRank = [];
%
% number of active contacts, M_lambda
  n = size(M, 1);
%
  N = dec2bin(0:(2^n-1),n);
  Ones = N=='1';
%
% check whether M_lambda is zero, and M_hat and q_hat are empty
  if n==0
    disp('M_hat is empty.  Invalid LCP in LCP_Neutral.m')
  else
%
%   search through all 2^n possible solutions
    for r = 1:2^n
%     construct the complementary matrix "Ar"
      Ar = eye(n);
      Jc = find(Ones(r,:));
      Ar(:,Jc) = -M(:,Jc);
%
%     if mod(r,5000)==0; 
%       disp(sprintf('Combination %6i of %6i in LCP_Neutral.m', r, 2^n));
%     end
%
%     check whether the complementar matrix is singular
      RCOND = rcond(Ar);
%
      if RCOND < EPS
%       the matrix is singular and has a null space not equal to the zero vector
%
%       find the null space
        NullAr = null(Ar, EPS);
%
%       rank of the null space of Ar
        RankNullAr = size(NullAr, 2);
%
        if RankNullAr == 1
%         the matrix is non-Ro, if and only if the single null-vector is
%         all-nonpositive or all-nonnegative
%
          if all(NullAr  > -EPS) || all(NullAr  < EPS)
%           the single-vector null space satisifies the Mohan criterion
%           non-Ro
%
            if all(NullAr  < EPS)
              NullAr = -NullAr;
            end
%
%           the zero and non-negative "x" values
            x = zeros(n,1);
            x(Jc) = NullAr(Jc);
%
%           append the neutral equilibrium solution "x" to the other solutions
            NeutralVec = [NeutralVec x];
            NeutralGrp = [NeutralGrp, r];
            NeutralRank = [NeutralRank, RankNullAr];
          end
        elseif RankNullAr > 1
%         we must consider the case a null-space in which each spanning
%         vector might not be a solution of the LCP, but combinaitons of
%         these vectors can be a solution
%
          if any(all(NullAr  > -EPS) | all(NullAr  < EPS))
%           check the obvious.  Note that this test also catches null vectors
%           that like on the boundary of the non-negative orthant, which is a
%           situation that will not be caught solving a linear program
%
            iNullAr = find(all(NullAr  > -EPS) | all(NullAr  < EPS));
%
            for i = iNullAr
              if all(NullAr(:,i)  < EPS)
                NullAr(:,i) = -NullAr(:,i);
              end
%
%             the zero and non-negative "x" values
              x = zeros(n,1);
              x(Jc) = NullAr(Jc);
%
%             append the neutral equilibrium solution "x" to the other solutions
              NeutralVec = [NeutralVec x];
              NeutralGrp = [NeutralGrp, r];
              NeutralRank = [NeutralRank, RankNullAr];
            end
          else
%           we will determine the existance of an intersection of the null
%           space with the non-negative orthant that is not the zero (trivial)
%           solution
%
%           we seek a solution of the linear program:
%               minimize  c'*NullAr*z
%               s.t.      M*z >= 1
%                         z unrestricted
%           where x = M*z
%
%           we use the octave glpk function in the "optimization" package
%
            c = abs(rand(size(NullAr,1),1));
            b = ones(size(NullAr,1),1);
            lb = []; ub = [];
%
%           constraints are lower bounds
            ctype = ''; 
            for i=1:size(NullAr,1); ctype = cstrcat(ctype,'L'); end
%
%           z is a continuous variable
            vartype = ''; 
            for i=1:size(NullAr,2); vartype = cstrcat(vartype,'C'); end
%
            param.msglev = 1;
            param.itlim = 100;
            s = 1;             % minimization, not maximization
%
            [ zopt, fmin, errnum, extra] = ...
              glpk(c'*NullAr, NullAr, b, lb, ub, ctype, vartype, s, param);
%
            if all(~isna(zopt))
%             a solution was found. Although a solution cone of dimension 
%             greater than 1 (larger than a ray) might exist, the linear
%             program only gives a single solution.  Report this solution
%
              y = NullAr * zopt;
%
%             the zero and non-negative "x" values
              x = zeros(n,1);
              x(Jc) = y(Jc);
%
%             append the neutral equilibrium solution "x" to the other solutions
              NeutralVec = [NeutralVec x];
              NeutralGrp = [NeutralGrp, r];
              NeutralRank = [NeutralRank, RankNullAr];
%
            end
          end % ~any(all(NullAr  > -EPS) | all(NullAr  < EPS))
        end % RankNullAr == 0
      end % if RCOND < EPS
    end % for r = 1:2^n
  end % if ~(n==0)
