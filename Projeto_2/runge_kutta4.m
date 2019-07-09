function y = runge_kutta4(x, ts, A, B, E)
  dx = A * x + B * E;
  k1 = ts * dx;
  aux = x + 0.5 * k1;
  
  dx = A * aux + B * E;
  k2 = ts * dx;
  aux = x + 0.5 * k2;
  
  dx = A * aux + B * E;
  k3 = ts * dx;
  aux = x + k3;
  
  dx = A * aux + B * E;
  k4 = ts * dx;
  
  y = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end

