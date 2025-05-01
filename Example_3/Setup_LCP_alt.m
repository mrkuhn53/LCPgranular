function [ M_hat, q_hat, H_xx_BD_inv, RCond_BD, ...
           M_alt, q_alt, H_xx_BD_inv_wo_constr] = ...
         Setup_LCP_alt( ...
           H_x_x, H_x_lambda, H_lambda_x, H_lambda_lambda, ...
           P_L, P_Lperp, C_MP_inv, dp_dstress, dc_w_dxcell, ...
           P_L_wo_constr, P_Lperp_wo_constr, C_MP_inv_wo_constr, dc_wo_constr);
%
%   reciprocal condition number
    RCond_BD = rcond(H_x_x*P_L + P_Lperp);
%
%   pre-Bott-Duffin inverse. Output variables are local to this function.
%   Use the pseudoinverse "pinv" with a tolerance as the 2nd argument
    Pre_H_xx_BD_inv = pinv(H_x_x*P_L + P_Lperp, 10*eps);
%
%   Bott-Duffin inverse
    H_xx_BD_inv = P_L * Pre_H_xx_BD_inv;
%
%   the M and q matrices of the LCP problem
    M_hat = H_lambda_lambda - H_lambda_x * H_xx_BD_inv * H_x_lambda;
%
    q_hat =   H_lambda_x * (C_MP_inv - H_xx_BD_inv * H_x_x * C_MP_inv) ...
              * dc_w_dxcell ...
            + H_lambda_x * H_xx_BD_inv * dp_dstress;
%
%   pre-Bott-Duffin inverse. Output variables are local to this function.
%   Use the pseudoinverse "pinv" with a tolerance as the 2nd argument
    Pre_H_xx_BD_inv_wo_constr = ...
      pinv(H_x_x*P_L_wo_constr + P_Lperp_wo_constr, 10*eps);
%
%   Bott-Duffin inverse
    H_xx_BD_inv_wo_constr = P_L_wo_constr * Pre_H_xx_BD_inv_wo_constr;
%
%   the M and q matrices of the LCP problem
    M_alt = H_lambda_lambda - H_lambda_x * H_xx_BD_inv_wo_constr * H_x_lambda;
%
    q_alt =   H_lambda_x ...
              * (C_MP_inv_wo_constr ...
                 - H_xx_BD_inv_wo_constr * H_x_x * C_MP_inv_wo_constr) ...
              * dc_wo_constr ...
            + H_lambda_x * H_xx_BD_inv_wo_constr * dp_dstress;
