% PARTE D: BÚSQUEDA DE RAÍCES Y ANÁLISIS DE SENSIBILIDAD
clear; clc; close all;
format longG;

% 1. Datos Originales del Experimento
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, ...
     460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, ...
     1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];

Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, ...
     132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, ...
     135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

N = length(f);
Z_th = 150; % Umbral crítico de impedancia

%% D.1: RECONSTRUCCIÓN DEL SPLINE NATURAL
A = zeros(N-2, N-2); B = zeros(N-2, 1);
for row = 1:(N-2)
    idx = row + 1;
    h_left  = f(idx) - f(idx-1);
    h_right = f(idx+1) - f(idx);
    if row > 1, A(row, row-1) = h_left; end
    A(row, row) = 2 * (f(idx+1) - f(idx-1));
    if row < N-2, A(row, row+1) = h_right; end
    B(row) = (6 / h_right) * (Z(idx+1) - Z(idx)) + (6 / h_left) * (Z(idx-1) - Z(idx));
end
M = [0; A \ B; 0]; 

%% D.2: BÚSQUEDA DE RAÍCES MEDIANTE MÉTODOS NUMÉRICOS
tol = 1e-6;
max_iter = 100;

% --- RAÍZ 1: BAJA FRECUENCIA (Intervalo [100, 120] Hz) ---
% A) Método de Bisección
[r1_bis, iter1_bis, tab1_bis] = bisection_method(100, 120, tol, max_iter, f, Z, M, Z_th);
% B) Método de Newton-Raphson (Aproximación inicial x0 = 110)
[r1_new, iter1_new, tab1_new] = newton_method(110, tol, max_iter, f, Z, M, Z_th);

% --- RAÍZ 2: ALTA FRECUENCIA (Intervalo [2160, 2340] Hz) ---
% A) Método de Bisección
[r2_bis, iter2_bis, tab2_bis] = bisection_method(2160, 2340, tol, max_iter, f, Z, M, Z_th);
% B) Método de Newton-Raphson (Aproximación inicial x0 = 2200)
[r2_new, iter2_new, tab2_new] = newton_method(2200, tol, max_iter, f, Z, M, Z_th);

%% D.3: CÁLCULO DE SENSIBILIDAD EN LA RAÍZ DE ALTA FRECUENCIA
[~, dZ_df_raiz2] = eval_spline(r2_new, f, Z, M);
sensibilidad_f = 1 / dZ_df_raiz2; % df/d|Z|

%% MOSTRAR TABLAS DE CONVERGENCIA EN CONSOLA
fprintf('========================================================================\n');
fprintf('        COMPARATIVA DE CONVERGENCIA: RAÍZ 1 (BAJA FRECUENCIA)\n');
fprintf('========================================================================\n');
fprintf('MÉTODO DE BISECCIÓN:\n');
fprintf('Iter\t\t\tX_mid\t\t\t\tF(X_mid)\n');
for i = 1:size(tab1_bis, 1)
    fprintf('%2d\t\t\t%11.5f\t\t\t%11.5e\n', tab1_bis(i,1), tab1_bis(i,2), tab1_bis(i,3));
end
fprintf('\nMÉTODO DE NEWTON-RAPHSON:\n');
fprintf('Iter\t\t\tX_n\t\t\t\tF(X_n)\n');
for i = 1:size(tab1_new, 1)
    fprintf('%2d\t\t\t%11.5f\t\t\t%11.5e\n', tab1_new(i,1), tab1_new(i,2), tab1_new(i,3));
end

fprintf('\n========================================================================\n');
fprintf('        COMPARATIVA DE CONVERGENCIA: RAÍZ 2 (ALTA FRECUENCIA)\n');
fprintf('========================================================================\n');
fprintf('MÉTODO DE BISECCIÓN:\n');
fprintf('Iter\t\t\tX_mid\t\t\t\tF(X_mid)\n');
for i = 1:size(tab2_bis, 1)
    fprintf('%2d\t\t\t%11.5f\t\t\t%11.5e\n', tab2_bis(i,1), tab2_bis(i,2), tab2_bis(i,3));
end
fprintf('\nMÉTODO DE NEWTON-RAPHSON:\n');
fprintf('Iter\t\t\tX_n\t\t\t\tF(X_n)\n');
for i = 1:size(tab2_new, 1)
    fprintf('%2d\t\t\t%11.5f\t\t\t%11.5e\n', tab2_new(i,1), tab2_new(i,2), tab2_new(i,3));
end

%% REPORTE GENERAL DE RESULTADOS
fprintf('\n========================================================================\n');
fprintf('                RESULTADOS \n');
fprintf('========================================================================\n');
fprintf('BANDA DE OPERACIÓN SEGURA (Impedancia < 150 Ohmios):\n');
fprintf('Límite Inferior (Raíz 1): %.4f Hz\n', r1_new);
fprintf('Límite Superior (Raíz 2): %.4f Hz\n', r2_new);
fprintf('------------------------------------------------------------------------\n');
fprintf('ANÁLISIS DE SENSIBILIDAD EN LA RAÍZ 2 (f = %.2f Hz):\n', r2_new);
fprintf('Pendiente de la curva d|Z|/df          : %.6f Ohmios/Hz\n', dZ_df_raiz2);
fprintf('Sensibilidad inversa df/d|Z|           : %.4f Hz/Ohmio\n', sensibilidad_f);
fprintf('========================================================================\n');

%% GRÁFICA DE RESPALDO DE INTERSECCIÓN
f_fina = linspace(100, 2730, 2000);
Z_fina = zeros(size(f_fina));
for i=1:length(f_fina)
    Z_fina(i) = eval_spline(f_fina(i), f, Z, M);
end

fig = figure('Name', 'Busqueda Raices', 'Renderer', 'painters');
plot(f_fina, Z_fina, 'b-', 'LineWidth', 2); hold on;
plot(f, Z, 'ko', 'MarkerFaceColor', 'y');
plot([100, 2730], [Z_th, Z_th], 'r--', 'LineWidth', 1.5);
plot([r1_new, r2_new], [Z_th, Z_th], 'gX', 'MarkerSize', 12, 'LineWidth', 3);
grid on;
xlabel('Frecuencia f (Hz)'); ylabel('Impedancia |Z| (\Omega)');
title('Identificación de la Banda de Operación Segura (|Z| = 150 \Omega)');
legend('Spline Cúbico', 'Datos Medidos', 'Umbral Z_{th} = 150 \Omega', 'Raíces Halladas');
saveas(fig, 'Banda_Segura_Raices.png');


%% =========================================================================
%   FUNCIONES LOCALES AUXILIARES
% =========================================================================
function [val, der] = eval_spline(x, f, Z, M)
    N = length(f);
    k = find(f <= x, 1, 'last');
    if isempty(k), k = 1; end
    if k == N, k = N - 1; end
    h = f(k+1) - f(k);
    
    val = (M(k)/(6*h))*(f(k+1)-x)^3 + (M(k+1)/(6*h))*(x-f(k))^3 + ...
          (Z(k)/h - (M(k)*h)/6)*(f(k+1)-x) + ...
          (Z(k+1)/h - (M(k+1)*h)/6)*(x-f(k));
      
    der = -(M(k)/(2*h))*(f(k+1)-x)^2 + (M(k+1)/(2*h))*(x-f(k))^2 - ...
          (Z(k)/h - (M(k)*h)/6) + ...
          (Z(k+1)/h - (M(k+1)*h)/6);
end

function [root, iter, history] = bisection_method(a, b, tol, max_iter, f, Z, M, Z_th)
    history = [];
    for iter = 1:max_iter
        x_mid = (a + b) / 2;
        f_mid = eval_spline(x_mid, f, Z, M) - Z_th;
        history = [history; iter, x_mid, f_mid];
        
        if abs(f_mid) < tol || (b - a)/2 < tol
            root = x_mid; return;
        end
        
        f_a = eval_spline(a, f, Z, M) - Z_th;
        if sign(f_mid) == sign(f_a)
            a = x_mid;
        else
            b = x_mid;
        end
    end
    root = x_mid;
end

function [root, iter, history] = newton_method(x0, tol, max_iter, f, Z, M, Z_th)
    x = x0; history = [];
    for iter = 1:max_iter
        [f_val, d_val] = eval_spline(x, f, Z, M);
        f_val = f_val - Z_th;
        history = [history; iter, x, f_val];
        
        if abs(f_val) < tol
            root = x; return;
        end
        x = x - f_val / d_val;
    end
    root = x;
end