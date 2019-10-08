% @file runge_kutta4.m to Simaan + ECG + LVAD
% @author Edvonaldo Horácio (edvonaldohs@ic.ufal.br)
% @brief
% @version 0.1
% @date 2019-07-07
% 
% @copyright Copyright (c) IC 2019
function y = runge_kutta4(x, ts, A, B, E, U, W)
  dx = A * x + B * E + U * W^2;
  k1 = ts * dx;
  aux = x + 0.5 * k1;
  
  dx = A * aux + B * E + U * W^2;
  k2 = ts * dx;
  aux = x + 0.5 * k2;
  
  dx = A * aux + B * E + U * W^2;
  k3 = ts * dx;
  aux = x + k3;
  
  dx = A * aux + B * E + U * W^2;
  k4 = ts * dx;
  
  y = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end

