% @file ECG_current.m
% @author Edvonaldo Horácio (edvonaldohs@ic.ufal.br)
% @brief
% @version 0.1
% @date 2019-07-07
% 
% @copyright Copyright (c) IC 2019
function y = ECG_current(x, z_0, w)
  alpha = 1 - sqrt(x(1)^2 + x(2)^2);
  a = [1.2, -5, 30, -7.5, 0.75];
  b = [0.25, 0.1, 0.1, 0.1, 0.4];
  th = [-1/3, -1/12, 0, 1/12, 1/2] * pi;
  theta = atan2(x(2), x(1));
  z_diff = x(3) - z_0;
  
  x_ECG = alpha * x(1) - w * x(2);
  y_ECG = alpha * x(2) + w * x(1);
  z_ECG = 0;
  
  for i = 1:5
    delta_theta = theta - th(i);
    z_ECG = z_ECG - a(i) * delta_theta * exp(-(delta_theta^2)/(2 * b(i)^2));
  end
  
  z_ECG = z_ECG - z_diff;
  
  y = [x_ECG; y_ECG; z_ECG];
end 