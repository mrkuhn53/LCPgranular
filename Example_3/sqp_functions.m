1;  % a dummy statement, so that functions can be defined
%
% a library of functions that will be used in Octave's sqp() function
%
% objective function
  function obj = Obj(x);
    global M_hat q_hat;
    obj = x' * (M_hat*x + q_hat);
  endfunction
%
% gradient of objective function
  function objGrad = ObjGrad(x);
    global M_hat q_hat
    objGrad = (M_hat + M_hat')*x + q_hat;
  endfunction
%
% gradient of objective function
  function objHess = ObjHess(x);
    global M_hat q_hat
    objHess = M_hat + M_hat';
  endfunction
%
% constraint function
  function h = Con(x)
    global M_hat q_hat
    h = M_hat*x + q_hat;
  endfunction
%
% gradient of constraint function
  function conGrad = ConGrad(x);
    global M_hat q_hat
    conGrad = M_hat;
  endfunction
%
  PHI = cell(3,1);
  PHI{1} = @Obj;
  PHI{2} = @ObjGrad;
  PHI{3} = @ObjHess;
%
  Hc = cell(2,1);
  Hc{1} = @Con;
  Hc{2} = @ConGrad;
