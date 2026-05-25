% SUBPARTE B2: INTERPOLACIÓN POR SPLINES CÚBICOS NATURALES
clear; clc; close all;
format longG;

% 1. Datos Originales del Experimento
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, ...
     460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, ...
     1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];

Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, ...
     132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, ...
     135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

N = length(f); % Cantidad de puntos (30)

%% B2.1: CONSTRUCCIÓN MANUAL DE LA MATRIZ DE TRACERES (image_01cbce.png)
A = zeros(N-2, N-2);
B = zeros(N-2, 1);

for row = 1:(N-2)
    idx = row + 1; % Índice real en los vectores de datos (puntos internos)
    
    h_left  = f(idx) - f(idx-1);
    h_right = f(idx+1) - f(idx);
    
    % Llenado de la matriz tridiagonal A (Coeficientes de f'')
    if row > 1
        A(row, row-1) = h_left;
    end
    
    A(row, row) = 2 * (f(idx+1) - f(idx-1));
    
    if row < N-2
        A(row, row+1) = h_right;
    end
    
    % Término independiente B (image_01cbce.png)
    B(row) = (6 / h_right) * (Z(idx+1) - Z(idx)) + (6 / h_left) * (Z(idx-1) - Z(idx));
end

% Resolver el sistema para los puntos internos
f_double_prime_internal = A \ B;

% Añadir las condiciones de Spline Natural: f''(x_1) = 0 y f''(x_n) = 0
M = [0; f_double_prime_internal; 0];

%% B2.2: EVALUACIÓN EN MALLA FINA MEDIANTE LA FÓRMULA DE LA DIAPOSITIVA (image_370655.png)
f_fina = linspace(100, 2730, 2000);
Z_spline = zeros(size(f_fina));

for j = 1:length(f_fina)
    x_val = f_fina(j);
    
    % Encontrar a qué segmento pertenece x_val
    if x_val >= f(end)
        k = N - 1;
    else
        k = find(f <= x_val, 1, 'last');
        if k == N, k = N - 1; end
    end
    
    % Implementación exacta de la fórmula de la diapositiva para P_i(x)
    h = f(k+1) - f(k);
    
    term1 = (M(k) / (6 * h)) * (f(k+1) - x_val)^3;
    term2 = (M(k+1) / (6 * h)) * (x_val - f(k))^3;
    term3 = (Z(k)/h - (M(k) * h)/6) * (f(k+1) - x_val);
    term4 = (Z(k+1)/h - (M(k+1) * h)/6) * (x_val - f(k));
    
    Z_spline(j) = term1 + term2 + term3 + term4;
end

%% CÁLCULO DE VALOR EN f = 1000 Hz
target_f = 1000;
k_target = find(f <= target_f, 1, 'last');
h_t = f(k_target+1) - f(k_target);

Z_1000_spline = (M(k_target) / (6 * h_t)) * (f(k_target+1) - target_f)^3 + ...
                (M(k_target+1) / (6 * h_t)) * (target_f - f(k_target))^3 + ...
                (Z(k_target)/h_t - (M(k_target) * h_t)/6) * (f(k_target+1) - target_f) + ...
                (Z(k_target+1)/h_t - (M(k_target+1) * h_t)/6) * (target_f - f(k_target));

%% COMPARACIÓN GRÁFICA CON EL POLINOMIO DE GRADO 5
[p5, ~, mu5] = polyfit(f, Z, 5);
Z_p5 = polyval(p5, f_fina, [], mu5);
Z_p5_1000 = polyval(p5, target_f, [], mu5);

figure('Name', 'Polinomio de grado 5 vs Spline Cúbico', 'NumberTitle', 'off');
plot(f, Z, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Datos Reales'); hold on;
plot(f_fina, Z_p5, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Polinomio Grado 5');
plot(f_fina, Z_spline, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Spline Cúbico Natural');
plot(target_f, Z_1000_spline, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Interp. Spline 1000 Hz');
grid on;
xlabel('Frecuencia f (Hz)', 'FontWeight', 'bold');
ylabel('Impedancia |Z| (\Omega)', 'FontWeight', 'bold');
title('Polinomio de grado 5 vs. Splines Cúbicos');
legend('Location', 'best');

% Mostrar resultados en consola
fprintf('========================================================================\n');
fprintf('         RESULTADOS COMPARATIVOS EN f = 1000 Hz\n');
fprintf('========================================================================\n');
fprintf('Impedancia estimada con Polinomio Grado 5 : %.4f ohmios\n', Z_p5_1000);
fprintf('Impedancia estimada con Spline Cúbico     : %.4f ohmios\n', Z_1000_spline);
fprintf('Diferencia absoluta entre modelos         : %.4f ohmios\n', abs(Z_p5_1000 - Z_1000_spline));
fprintf('========================================================================\n');