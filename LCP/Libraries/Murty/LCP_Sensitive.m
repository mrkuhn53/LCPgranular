function [Sensitive] = LCP_Sensitive(M, q, Lambda, EPS);
%
% determine whether the Lambda solutions of the LCP(q,M) are path-sensitive
%
% the following LCP must be solved beforehand, with the solutions stored
% in the matrix Lambda
%   find "lambda" such that  M*lambda >= 0; lambda >=0; lambda'*(M*lambda) = 0
%
% Input: 
%   M =      square matrix, nxn
%   q =      column vector, nx1
%   Lambda = matrix, nxm, where "m" is the number of solutions
%   EPS =    small parameter for determining whether matrices are
%            non-singular:  when rcond() > EPS
%
% Output: 
%   Sensitive = row vector, 1xm, whether the ith solution, Lambda(:,i),
%               is path-sensitive (=1) or not (=0)
%
  Sensitive = -ones(1,size(Lambda,2));
%
% number of active contacts, M_lambda
  n = size(M, 1);
%
% number of solutions
  nLambda = size(Lambda, 2);
%
  if nLambda > 0
    w = M*Lambda + q*ones(1,nLambda);
    I_Type = 1000*ones(n,nLambda);
%
    I_Type(find(Lambda < -EPS | w < -EPS)) =  0; % not a solution
    I_Type(find(Lambda < EPS & w > EPS)) =  1;
    I_Type(find(Lambda > EPS & w < EPS)) =  2;
    I_Type(find(Lambda < EPS & w < EPS)) =  3; % degenerate solution
%
    for iLambda = 1:nLambda
%     for each solution "lambda", check each of 3 criteria for regularity
%
      if any(Lambda(:,iLambda) < -EPS | w(:,iLambda) < -EPS)
%       not even a solution, just claim that it is regular, so that it is
%       ignored later
        Sensitive(iLambda) = 0;
      else
%       the indices of the three types of lambda values
        I_1 = find(I_Type(:,iLambda)==1)';
        I_2 = find(I_Type(:,iLambda)==2)';
        I_3 = find(I_Type(:,iLambda)==3)';
%
%       first test of regularity
        Regular = isempty(I_2) || rcond(M(I_2,I_2)) > EPS;
%
%       second test of regularity
        Regular = Regular ...
                  && ...
                  (isempty([I_2, I_3]) || rcond(M([I_2, I_3],[I_2, I_3])) > EPS);
%
%       third test of regularity
        if Regular && ~isempty(I_3)
          A = [ M(I_2,I_2)',                     M(I_3,I_2)';
                M(I_2,I_3)',                     M(I_3,I_3)'];
%
          c = [zeros(length(I_2),1);
                ones(length(I_3),1)];
%
%         arranged as
%           u = u([I_2, I_3],1)
%
          C = c;
%
          B = zeros(length([I_2, I_3]), 1);
%
%
          LB = -Inf(length([I_2, I_3]), 1);
          UB = [Inf(length([I_2]), 1); zeros(length([I_3]), 1)];
%
          CTYPE = '';
          for i = 1:length(I_2); CTYPE = cstrcat(CTYPE, 'S'); end
          for i = 1:length(I_3); CTYPE = cstrcat(CTYPE, 'L'); end
          %
          VARTYPE = '';
          for i = 1:length(I_2); VARTYPE = cstrcat(VARTYPE, 'C'); end
          for i = 1:length(I_3); VARTYPE = cstrcat(VARTYPE, 'C'); end
          %
          SENSE =  1; % objective function is minimized
%
          PARAM.msglev = 1;
          PARAM.outfrq = 0;
          PARAM.itlim = 100;
          PARAM.presol = 0;
%
%         C, A, B, LB, UB, CTYPE, VARTYPE
          [XOPT, FMIN, ERRNUM, EXTRA] = ...
            glpk(C, A, B, LB, UB, CTYPE, VARTYPE, SENSE, PARAM);
%         XOPT
%         FMIN
%         ERRNUM
%         EXTRA
%         EXTRA.status
%
          Regular = Regular && (EXTRA.status==5 & all(XOPT==0))
        end
%
        Sensitive(iLambda) = ~Regular;
      end
%
    end
  else
    Sensitive = [];
  end
