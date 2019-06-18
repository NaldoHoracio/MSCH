% @file Simaan_model.m
% @author Edvonaldo Hor�cio (edvonaldohs@ic.ufal.br)
% @brief
% @version 0.1
% @date 2019-06-17
% 
% @copyright Copyright (c) IC 2018

clc;
clear all;
%% PAR�METROS DO MODELO SIMMAN ET. AL (2009)
% Resist�ncias (mmHg.s/ml)
R_s = 1.0000;% Systemic Vascular Resistance (SVR) - R_s = R_s
R_m = 0.0050;% Mitral Valve Resistance - R_m = R_m
R_a = 0.0010;% Aortic Valve - R_a = R_a
R_c = 0.0398;% Characteristic Resistence - R_c = R_c

% Capacit�ncias (ml/mmHg)
C_ae = 4.4000;% Left Ventricular Compliance - C_ae = C_r
C_s = 1.3300;% Left Atrial Compliance - C_s = C_s
C_ao = 0.0800;% Aortic  Compliance - C_ao = C_a

% Indut�ncia (mmHg.s�/ml)
L = 0.0005;% Inertance of blood in Aorta - L = L_s

% Outras constantes
E_min = 0.06;% Elast�ncia m�nima
E_max = 2;% Elast�ncia m�xima
FC = 75;% Frequ�ncia card�aca - FC = HR
T = 60/FC;% Intervalo de tempo de um ciclo card�aco - T = t_c
T_max = 0.2 + 0.15*T;% Tempo m�ximo
%t_n = T/T_max;% Tempo normalizado p/ um ciclo normalizado
V_0 = 10;% Volume inicial da Aorta (mL)
% Diodos (ideais)
% D_a = 1;
% D_m = 1;

%% DOUBLE HILL - FUN��O ELAST�NCIA NORMALIZADA
E_n = @(t_n) 1.55 * (((t_n/0.7)^1.9)/(1 + (t_n/0.7)^1.9)) * (1/(1 + (t_n/1.17)^21.9));


% FUN��O ELAST�NCIA
E_t = @(t_begin)(E_max - E_min)*E_n(mod(t_begin, T)/T_max) + E_min;

%  PLOT DA FUN��O ELAST�NCIA
% t_elast = 0:0.01:1;
% DH_n = 1.55 .* (((t_elast./0.7).^1.9)./(1 + (t_elast./0.7).^1.9)) .* (1./(1 + (t_elast./1.17).^21.9));
% EL_t = (E_max - E_min) .* DH_n + E_min;
% plot(t_elast,EL_t);

%%
% PRESS�O NO VENTR�CULO ESQUERDO - P_ve
P_ve = @(t_begin,V_ve)E_t(t_begin) * (V_ve - V_0);

% ESPA�O DE ESTADOS - MATRIZES
% Matriz A
A = @(t_begin, D_m, D_a) [-(D_m/R_m + D_a/R_a)*E_t(t_begin), D_m/R_m, 0, D_a/R_a, 0;
    (D_m * E_t(t_begin))/(R_m * C_ae), -(1/C_ae)*(1/R_s + D_m/R_m), 0, 0, (1/(R_s*C_ae));
    0, 0, -R_c/L, 1/L, -(1/L);
    (D_a * E_t(t_begin))/(R_a*C_ao), 0, -(1/C_ao), -D_a/(R_a*C_ao), 0;
    0, (1/(R_s*C_s)), 1/C_s, 0, -(1/(R_s*C_s))];
% Vetor p
p = @(t_begin, D_m, D_a) [(D_m/R_m + D_a/R_a)*E_t(t_begin)*V_0;
    -(D_m*E_t(t_begin)*V_0)/(R_m*C_ae);
    0;
    -(D_a*E_t(t_begin)*V_0)/(R_a*C_ao);
    0];

%% VALORES DA SIMULA��O
t_begin = 0;
t_end = 5;
step = 10^-4;
lines_matrix = t_end/step;

% Iniciando os dados com valor zero
simaan.x.V_ve = zeros(lines_matrix, 1);
simaan.x.P_ae = zeros(lines_matrix, 1);
simaan.x.Q_a = zeros(lines_matrix, 1);
simaan.x.P_ao = zeros(lines_matrix, 1);
simaan.x.P_s = zeros(lines_matrix, 1);
simaan.P_ve = zeros(lines_matrix, 1);
simaan.t = t_begin:step:t_end;% Tempo da simula��o

% Condi��es iniciais - bound rates
simaan.x.V_ve(1) = 140;% Volume no Ventr�culo Esquerdo - LVV
simaan.x.P_ae(1) = 5;% Press�o no �trio Esquerdo
simaan.x.Q_a(1) = 0;% Fluxo na Aorta
simaan.x.P_ao(1) = 90;% Press�o na Aorta
simaan.x.P_s(1) = 90;% Press�o Sist�mica
simaan.P_ve(1) = P_ve(t_begin, simaan.x.V_ve(1));% Press�o no Ventr�culo Esquerdo
% Diodos (ideais) na condi��o inicial
D_a = 1;
D_m = 1;
% %%
% clear all
% 
% %% TESTE - RUNGE-KUTTA
% t_begin = 0;
% t_end = 3;
% step = 10^-4;
% 
% index = 1;
% k_begin = 0;
% for i =  k_begin:step:(t_end-step)
%     t_begin = plus(t_begin,step);
%     index = index + 1;
% end

%% EXECUTANTANDO A SIMULA��O - M�TODO DE RUNGE-KUTTA
k_begin = 0;
idx = 1;
for i = k_begin:step:(t_end-step)
    
    t_begin = plus(t_begin, step);
    
    % Valores dos diodos (bloqueia ou n�o bloqueia a passagem)
    if (simaan.x.P_ao(idx) <= simaan.P_ve(idx))
        D_a = 1; D_m = 0;
    elseif (simaan.P_ve(idx) <= simaan.x.P_ae(idx))
        D_a = 0; D_m = 1;
    else D_m = 0; D_a = 0;
    end
    
    x = [simaan.x.V_ve(idx); simaan.x.P_ae(idx); simaan.x.Q_a(idx); simaan.x.P_ao(idx); simaan.x.P_s(idx)];
    
    % Runge - Kutta de 4� Ordem
    % Valores iniciais
    dx = A(t_begin, D_m, D_a) * x +  p(t_begin, D_m, D_a);
    kx1 = step * dx;
    x1 = x + 0.5*kx1;% Valor m�dio de x1
    
    % Incremento
    dx = A(t_begin, D_m, D_a) * x1 +  p(t_begin, D_m, D_a);
    kx2 = step * dx;
    x1 = x + 0.5*kx2;% Valor m�dio de x2
    
    % Incremento
    dx = A(t_begin, D_m, D_a) * x1 +  p(t_begin, D_m, D_a);
    kx3 = step * dx;
    x1 = x + kx3;% Apenas x3
    
    % Incremento
    dx = A(t_begin, D_m, D_a) * x1 +  p(t_begin, D_m, D_a);
    kx4 = step * dx;% Valor de k4 de acordo com o passo de 10^-4
    
    % Valor final de x (por aproxima��o)
    xf = x + (kx1 + 2 * kx2 + 2 * kx3 + kx4)/6;
    
    % Armazenando os dados
    simaan.x.V_ve(idx+1) = xf(1);
    simaan.x.P_ae(idx+1) = xf(2);
    simaan.x.Q_a(idx+1) = xf(3);
    simaan.x.P_ao(idx+1) = xf(4);
    simaan.x.P_s(idx+1) = xf(5);
    
    simaan.P_ve(idx+1) = P_ve(t_begin, simaan.x.V_ve(idx+1));
    idx = idx+1;
end

%% RESULTADOS
xmin = 0;
xmax = 5;
figure(1)

subplot(3,1,1);
plot(simaan.t, simaan.x.P_ao, '-.r', simaan.t, simaan.P_ve, 'g', simaan.t, simaan.x.P_ae, '--b');
title('Simulated hemodynamic waveforms for a normal heart using Simaan');
legend('AoP', 'LVP', 'LAP');
ylabel('Pressure (mmHg)');
axis([xmin, xmax, -5, 150]);
grid on;

subplot(3,1,2);
plot(simaan.t, simaan.x.Q_a, 'r');
ylabel('Aortic Flow (ml/s)');
axis([xmin, xmax, -50, 700]);
grid on;

subplot(3,1,3);
plot(simaan.t, simaan.x.V_ve, 'r');
ylabel('LVV (ml)');
xlabel('time (s)');
axis([xmin, xmax, 50, 160]);
grid on;
