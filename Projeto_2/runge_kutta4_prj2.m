% @file runge_kutta4_prj2.m
% @author Edvonaldo Horácio (edvonaldohs@ic.ufal.br)
% @brief
% @version 0.1
% @date 2019-07-07
% 
% @copyright Copyright (c) IC 2019
function y = runge_kutta4_prj2(x, ts, z_0, w)
  dx = ECG_current(x, z_0, w);
  k1 = ts * dx;
  aux = x + 0.5 * k1;
  
  dx = ECG_current(aux, z_0, w);
  k2 = ts * dx;
  aux = x + 0.5 * k2;
  
  dx = ECG_current(aux, z_0, w);
  k3 = ts * dx;
  aux = x + k3;
  
  dx = ECG_current(aux, z_0, w);
  k4 = ts * dx;
  
  y = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end
