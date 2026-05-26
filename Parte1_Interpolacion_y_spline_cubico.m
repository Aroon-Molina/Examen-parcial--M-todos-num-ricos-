%% PARTE 1: INTERPOLACIÓN BÁSICA Y SPLINE CÚBICO
% Proyecto: Caracterización de Front-End Analógico y Biosensor Implantable
% Rol: Ingeniero Responsable de Procesamiento Numérico

clear; clc; close all;
format shortG;

%% 1. Base de Datos Experimentales (40 mediciones tabuladas)
f = [10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, ...
     45.0, 47.5, 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, ...
     77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5];

V = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, 1.216, 1.048, 0.866, 0.689, ...
     0.521, 0.364, 0.223, 0.103, 0.012, -0.041, -0.057, -0.034, 0.018, 0.096, 0.197, 0.318, 0.452, ...
     0.579, 0.700, 0.809, 0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004];

Z = [182.4, 178.9, 175.1, 171.0, 166.8, 162.7, 158.9, 155.4, 152.0, 149.0, 146.1, 145.2, 145.8, 147.3, ...
     149.9, 153.5, 158.0, 163.2, 168.9, 174.8, 180.5, 186.2, 191.5, 196.2, 200.1, 203.1, 205.2, ...
     206.3, 206.1, 204.7, 198.0, 194.4, 190.9, 187.8, 185.1, 183.0, 181.6, 180.8, 180.6, 180.9];

%% 2. Interpolación de Lagrange de Segundo Grado (Local)
% Entorno óptimo para 41.0 kHz: Nodos 37.5, 40.0 y 42.5 kHz
f_41_pts = [37.5, 40.0, 42.5];
V_41_pts = [1.048, 0.866, 0.689];
Z_41_pts = [145.2, 145.8, 147.3];

% Entorno óptimo para 73.0 kHz: Nodos 70.0, 72.5 y 75.0 kHz
f_73_pts = [70.0, 72.5, 75.0];
V_73_pts = [0.197, 0.318, 0.452];
Z_73_pts = [200.1, 203.1, 205.2];

% Función anónima para calcular el polinomio cuadrático de Lagrange
lagrange2 = @(x_p, y_p, x_t) ...
    y_p(1) * ((x_t - x_p(2))*(x_t - x_p(3))) / ((x_p(1) - x_p(2))*(x_p(1) - x_p(3))) + ...
    y_p(2) * ((x_t - x_p(1))*(x_t - x_p(3))) / ((x_p(2) - x_p(1))*(x_p(2) - x_p(3))) + ...
    y_p(3) * ((x_t - x_p(1))*(x_t - x_p(2))) / ((x_p(3) - x_p(1))*(x_p(3) - x_p(2)));

V_41_lag = lagrange2(f_41_pts, V_41_pts, 41.0);
Z_41_lag = lagrange2(f_41_pts, Z_41_pts, 41.0);
V_73_lag = lagrange2(f_73_pts, V_73_pts, 73.0);
Z_73_lag = lagrange2(f_73_pts, Z_73_pts, 73.0);

%% 3. Interpolación por Spline Cúbico Natural (Global)
V_41_spl = spline_natural_eval(f, V, 41.0);
Z_41_spl = spline_natural_eval(f, Z, 41.0);
V_73_spl = spline_natural_eval(f, V, 73.0);
Z_73_spl = spline_natural_eval(f, Z, 73.0);

%% 4. Despliegue de Resultados en Consola
fprintf('========================================================================\n');
fprintf('             REPORTE DE PROCESAMIENTO NUMÉRICO - PARTE 1                \n');
fprintf('========================================================================\n');
fprintf('Estimación en f = 41.0 kHz:\n');
fprintf('  -> Lagrange 2° Grado:  V = %7.4f V   |   |Z| = %7.2f Ohm\n', V_41_lag, Z_41_lag);
fprintf('  -> Spline Cúbico Nat:  V = %7.4f V   |   |Z| = %7.2f Ohm\n', V_41_spl, Z_41_spl);
fprintf('------------------------------------------------------------------------\n');
fprintf('Estimación en f = 73.0 kHz:\n');
fprintf('  -> Lagrange 2° Grado:  V = %7.4f V   |   |Z| = %7.2f Ohm\n', V_73_lag, Z_73_lag);
fprintf('  -> Spline Cúbico Nat:  V = %7.4f V   |   |Z| = %7.2f Ohm\n', V_73_spl, Z_73_spl);
fprintf('========================================================================\n');

%% 5. Graficación de Curvas de Respuesta Continua
f_mesh = linspace(min(f), max(f), 1000);
V_mesh = arrayfun(@(fi) spline_natural_eval(f, V, fi), f_mesh);
Z_mesh = arrayfun(@(fi) spline_natural_eval(f, Z, fi), f_mesh);

figure('Name', 'Respuesta del Front-End Analógico y Biosensor', 'Color', 'w');

subplot(2,1,1);
plot(f, V, 'r.', 'MarkerSize', 12); hold on;
plot(f_mesh, V_mesh, 'b-', 'LineWidth', 1.5);
plot([41.0, 73.0], [V_41_spl, V_73_spl], 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
title('Voltaje de Salida V(f) vs Frecuencia');
xlabel('Frecuencia (kHz)'); ylabel('Voltaje (V)'); grid on;
legend('Mediciones', 'Spline Cúbico', 'Estimaciones (41 y 73 kHz)', 'Location', 'best');

subplot(2,1,2);
plot(f, Z, 'm.', 'MarkerSize', 12); hold on;
plot(f_mesh, Z_mesh, 'g-', 'LineWidth', 1.5);
plot([41.0, 73.0], [Z_41_spl, Z_73_spl], 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
title('Impedancia Equivalente |Z(f)| vs Frecuencia');
xlabel('Frecuencia (kHz)'); ylabel('Impedancia (\Omega)'); grid on;
legend('Mediciones', 'Spline Cúbico', 'Estimaciones (41 y 73 kHz)', 'Location', 'best');

%% Función Interna: Constructor y Evaluador de Spline Cúbico Natural
function y_est = spline_natural_eval(x, y, x_target)
n = length(x);
h = diff(x);
A = zeros(n, n);
B = zeros(n, 1);

% Fronteras Naturales (Segunda derivada nula en los extremos)
A(1,1) = 1; 
A(n,n) = 1;

% Construcción de la matriz tridiagonal interna
for i = 2:n-1
    A(i, i-1) = h(i-1);
    A(i, i)   = 2 * (h(i-1) + h(i));
    A(i, i+1) = h(i);
    B(i) = 6 * ((y(i+1) - y(i))/h(i) - (y(i) - y(i-1))/h(i-1));
end
M = A \ B; % Resolución de los momentos vectoriales

% Localización del intervalo correspondiente
idx = find(x_target >= x, 1, 'last');
if isempty(idx), idx = 1; elseif idx == n, idx = n - 1; end

% Evaluación analítica mediante trazador
hi = h(idx);
term_A = (x(idx+1) - x_target)^3 / (6 * hi);
term_B = (x_target - x(idx))^3 / (6 * hi);
term_C = (y(idx) - (M(idx) * hi^2)/6) * (x(idx+1) - x_target) / hi;
term_D = (y(idx+1) - (M(idx+1) * hi^2)/6) * (x_target - x(idx)) / hi;

% --- LÍNEA CORREGIDA AQUÍ ---
y_est = term_A * M(idx) + term_B * M(idx+1) + term_C + term_D;
end