%
% bundle the functions that will be used in the Octave function sqp()
%
% the objective function and its gradient and Hessian
  PHI = cell(3,1);
  PHI{1} = @Obj;
  PHI{2} = @ObjGrad;
  PHI{3} = @ObjHess;
% 
% the inequality constraint function and its gradient
  Hc = cell(2,1);
  Hc{1} = @Con;
  Hc{2} = @ConGrad;
