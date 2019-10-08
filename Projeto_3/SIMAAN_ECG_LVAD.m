% @file SIMAAN_ECG_LVAD.m
% @author Edvonaldo Horácio (edvonaldohs@ic.ufal.br)
% @brief
% @version 0.1
% @date 2019-07-22
% 
% @copyright Copyright (c) IC 2019

clc;
clear all;

% PARÂMETROS INICIAIS
t_begin = 0;
t_step = 0.001;
t_end = 60;
t_step_simaan = 0.0001;
t_ECG = t_begin:t_step:t_end;
t_simaan = t_begin:t_step_simaan:t_end;
N = length(t_ECG);
quot = t_step / t_step_simaan;

% CTES. DO ECG SINTÉTICO E FILTRO MAMEMI PARA A DETECÇÃO DA AMPLITUDE EM R
amplitude = 0.15; 
HR_begin = 70;
sigma = 2;
delta = 2;
beta = 15;

% ELASTÂNCIA
var_aux = 0;
t_elementar = 0;
e_el = 0;

% CTES. DO MODELO SIMAAN
% Resistências do Modelo Simaan (mmHg.s/ml)
R_s = 1.0000;% Systemic Vascular Resistance (SVR) - R_s = R_s
R_m = 0.0050;% Mitral Valve Resistance - R_m = R_m
R_a = 0.0010;% Aortic Valve - R_a = R_a
R_c = 0.0398;% Characteristic Resistence - R_c = R_c
R_i = 0.0677;% Resistence of input of the Canula
R_o = 0.0677;% Resistence of output of the Canula

% Capacitâncias do Modelo Simaan (ml/mmHg)
C_ae = 4.4000;% Left Ventricular Compliance - C_ae = C_r
C_s = 1.3300;% Left Atrial Compliance - C_s = C_s
C_ao = 0.0800;% Aortic  Compliance - C_ao = C_a

% Indutância (mmHg.s²/ml)
L = 0.0005;% Inertance of blood in Aorta - L = L_s
L_i = 0.0127;% Inertancia de entrada da Canula
L_o = 0.0127;% Inertancia de saida da Canula
V_o = 10;

% Elastância
E_min = 0.06;% Elastância mínima
E_max = 2;% Elastância máxima

% OUTRAS CTES. DO LVAD
beta_zero = 0.1707;% Parametro de diferenca de pressao beta zero
beta_one = 0.02177;% Parametro de diferenca de pressao beta um
beta_two = -9.9025e-5;% Parametro de diferenca de pressao beta dois
W_min = 12e3 * 2 * pi/60;% Velocidade minima de rotacao do LVAD
W_max = 18e3 * 2 * pi/60;% Velocidade maxima de rotacao do LVAD
alpha = -3.5;% Parametro alfa da resistencia de succao
x1 = 1;% Parametro x1 barra da resistencia de succao

% VALORES INICIAIS DA ELASTÂNCIA, MODELO SIMAAN, VEL. DO LVAD E ECG
% SINTÉTICO
E_curr = E_min;% Elastância inicial
y = [140; 90; 0; 90; 5; 0];% Vetor com valores iniciais do Simaan: [V_ve, P_ao, Q_a, P_s, P_ae, Q_p]
ECG = [-1; 0; 0];
W_curr = W_min; %Valor inicial do Simaan


h = zeros(1,N);
max_index = zeros(1,N);
min_index = zeros(1,N);
a = zeros(1,N);
n = zeros(1,N+beta);
g = zeros(1,N);
r = zeros(1,N);
v = zeros(1,N);
w = zeros(1,N);
waveR_detect = zeros(1,N);

for index = 1:N-1
    HR = HR_begin + 2 * randn;
    RR(index) = 60 / HR;
    omega = 2 * pi / RR(index);
    z_0 = amplitude * sin(2*pi*t_ECG(index)*(60/(12+randn)));
    ECG(:,index+1) = runge_kutta4_prj2(ECG(:,index), t_step, z_0, omega);
    
    % DETECÇÃO DAS R-WAVES
    if ECG(3,index+1) > max_index(index)
        max_index(index+1) = max_index(index) + sigma * delta;
    else
        max_index(index+1) = max_index(index) - delta;
    end
    
    if ECG(3,index+1) < min_index(index)
        min_index(index+1) = min_index(index) - sigma * delta;
    else
        min_index(index+1) = min_index(index) + delta;
    end
    
    h(index+1) = ECG(3,index+1) - ((max_index(index+1) + min_index(index+1)) / 2);
    a(index+1) = (max_index(index+1) + min_index(index+1));
    
    if a(index+1) <= abs(h(index+1))
        n(index+1) = sign(h(index+1) * (abs(h(index+1)) - a(index+1)));
    else
        n(index+1) = 0;
    end
    
    if index > beta
        if n(index) > 0 && n(index) > n(index-beta) && n(index) > n(index+beta)
          g(index) = n(index) - max(n(index-beta),n(index+beta));
        elseif n(index) < 0 && n(index) < n(index-beta) && n(index) < n(index+beta)
          g(index) = n(index) + min(n(index-beta),n(index+beta));
        else
          g(index) = 0;
        end
            
        if g(index) > g(index-1) && g(index) > g(index+1)
          r(index) = g(index);
        else
          r(index) = 0;
        end
            
        if g(index) < g(index-1) && g(index) < g(index+1)
          v(index) = g(index);
        else
          v(index) = 0;
        end
    end
        
    if r(index) > 0
      w(index) = r(index);
    end
        
    if v(index) < 0
      w(index) = -v(index);
    end
        
    if w(index) == 1
      var_aux = 1;
    end
        
    if var_aux == 1 && ECG(3,index) >= 0.025
      var_aux = 2;
    end
        
    if var_aux == 2 && ECG(3,index) <= 0.025
      var_aux = 0;
      waveR_detect(index) = 1;
      e_el = 1;
    end
        % Elastância
    for j = 1:quot
       P_ve((index-1)*quot+j) = E_curr((index-1)*quot+j) * (y(1,(index-1)*quot+j) - V_o);
       P_ae = y(5,(index-1)*quot+j);
       P_ao = y(2,(index-1)*quot+j);
            
       if P_ae > P_ve((index-1)*quot+j)
         Dm = 1;
       else
         Dm = 0;
       end
       
       if P_ve((index-1)*quot+j) > P_ao
        Da = 1;
       else
        Da = 0;
       end
            
       if P_ve((index-1)*quot+j) > x1
         R_k = 0;
       else
         R_k = alpha * (P_ve((index - 1)*quot+j) - x1);
       end
        
       % GERAÇÃO DA ELASTÂNCIA A PARTIR DAS R-WAVE
       if e_el == 1
         if waveR_detect(index) == 1
           t_elementar = 1;
         end
         t_max = 0.2 + 0.15 * RR(index);
         t_n = t_elementar * t_step_simaan / t_max;
         tn_exp1 = (t_n/0.7) ^ 1.9;
         tn_exp2 = (t_n/1.17) ^ 21.9;
         E_n = 1.55 * (tn_exp1 / (1 + tn_exp1)) * (1 / (1 + tn_exp2));
         E_curr((index-1)*quot+j+1) = (E_max - E_min) * E_n + E_min;
         t_elementar = t_elementar + 1;
       else
          E_curr((index-1)*quot+j+1) = E_curr((index-1)*quot+j);
       end  
       W_curr((index-1)*quot+j+1) = W_min + ((W_max - W_min)/t_end) * t_simaan((index-1)*quot+j+1);
       c1 = Da/R_a;
       c2 = Dm/R_m;
       c3 = (c1 + c2);
       c4 = 1/C_ao;
       c5 = c1 * c4;
       c6 = 1/L;
       c7 = R_c/L;
       c8 = 1/C_s;
       c9 = c8/R_s;
       c10 = c2/C_ae;
       c11 = 1/(R_s * C_ae);
       c12 = (1/C_ae) * (1/R_s + c2);
       c13 = (c1 + c2) * V_o;
       c14 = c5 * V_o;
       c15 = c10 * V_o;
       % Outros parâmetros da matriz para gerar o LVAD
       L_x = L_i + L_o + beta_one;
       R_x = R_i + R_o + R_k + beta_zero;
       c16 = 1/L_x;
       c17 = R_x/L_x;
       c18 = c16 * V_o;
       c19 = beta_two/L_x;
       c20 = 1;
       % Matrizes do modelo   
       A = [-c3*E_curr((index-1)*quot+j), c1, 0, 0, c2, -c20;
            c5*E_curr((index-1)*quot+j), -c5, -c4, 0, 0, c4;
            0, c6, -c7, -c6, 0, 0;
            0, 0, c8, -c9, c9, 0;
            c10*E_curr((index-1)*quot+j), 0, 0, c11, -c12, 0;
            c16*E_curr((index-1)*quot+j), -c16, 0, 0, 0, -c17];
            
       B = [c13; -c14; 0; 0; -c15; -c18];
       
       U = [0; 0; 0; 0; 0; -c19];
      % Variáveis de estado    
      if index ~= N
        y(:,(index-1)*quot+j+1) = runge_kutta4(y(:,(index-1)*quot+j), ...
          t_step_simaan, A, B, E_curr((index-1)*quot+j), U,  W_curr((index-1)*quot+j));
      end
   end
end

P_ve(index*quot+1) = P_ve(index)

figure
subplot(2,1,1)
plot(t_ECG,ECG(3,:),'green',t_ECG,waveR_detect.*ECG(3,:),'-')
legend('ECG Signal')
xlabel('time (s)')
ylabel('mV (Volt)')
title('Synthetic ECG Signal')

subplot(2,1,2)
plot(t_simaan,E_curr,'black',t_ECG,waveR_detect.*max(E_curr),'--')
legend('Elastance Signal')
xlabel('time (s)')
ylabel('Elastance (mmHg/mL)')
title('Elastance of the ECG Signal')

figure
subplot(3,1,1)
ylabel('Press (mmHg)')
xlabel('time (s)')
hold on
plot(t_simaan, P_ve, 'red')
plot(t_simaan, y(5,:), 'green')
plot(t_simaan, y(2,:), 'blue')
legend('Left Ventricule', 'Left Atrium', 'Aort')
title('Params of Hemodynamic Wave Simulation of a Normal Heart - Subplot 1: Press')

subplot(3,1,2)
plot(t_simaan, y(3,:),'red')
ylabel('Flux on Aort (mL/s)')
xlabel('time (s)')
title('Params of Hemodynamic Wave Simulation of a Normal Heart - Subplot 2: Flux')

subplot(3,1,3)
plot(t_simaan, y(1,:), 'blue')
ylabel('Volume on Left Ventricule (mL)')
xlabel('time (s)')
title('Param of Hemodynamic Wave Simulation of a Normal Heart - Subplot 3: Volum')


%sgtitle('Hemodynamic Wave Simulation of a Normal Heart')

figure
plot(t_simaan, y(6,:))
ylabel('Pump flow (ml/s)')
xlabel('time (s)')
title('Pump flow to $R_s$ = 1.0 and RPM = 15500', 'Interpreter', 'latex')