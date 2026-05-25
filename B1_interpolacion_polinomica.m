% =========================================================================
% SUBPARTE B1: INTERPOLACIÓN POLINÓMICA Y FENÓMENO DE RUNGE
% =========================================================================
clear; clc; close all;
format longG; % Permite visualizar la máxima precisión numérica

% 1. Vectores de datos originales (30 puntos)
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, ...
    460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, ...
    1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];

Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, ...
    132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, ...
    135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

% Malla fina de frecuencias para evaluar y observar los gráficos continuamente
f_fina = linspace(100, 2730, 2000);

%% =======================================================================
% B1.1 & B1.2: AJUSTE DE POLINOMIOS Y DEMOSTRACIÓN DE RUNGE
% =======================================================================

% Polinomio de Grado 29 (Método Matricial / Ajuste exacto)
% Nota: MATLAB advertirá que la matriz está mal condicionada (Vandermonde)
p29 = polyfit(f, Z, 29);
Z_p29 = polyval(p29, f_fina);

% Polinomios de menor grado (Escalonados / Aproximaciones)
p5  = polyfit(f, Z, 5);   Z_p5  = polyval(p5, f_fina);
p10 = polyfit(f, Z, 10);  Z_p10 = polyval(p10, f_fina);
p15 = polyfit(f, Z, 15);  Z_p15 = polyval(p15, f_fina);

% Graficación comparativa
figure('Name', 'Evidencia del Fenómeno de Runge', 'NumberTitle', 'off');
plot(f, Z, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Datos Reales'); hold on;
plot(f_fina, Z_p5, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Grado 5');
plot(f_fina, Z_p10, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Grado 10');
plot(f_fina, Z_p15, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Grado 15');
plot(f_fina, Z_p29, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Grado 29 (Runge)');

grid on; ylim([100, 200]); % Ajustamos límites en Y para apreciar el fenómeno sin perder la escala
xlabel('Frecuencia f (Hz)', 'FontWeight', 'bold');
ylabel('Impedancia |Z| (\Omega)', 'FontWeight', 'bold');
title('Comparativa Polinómica: Evidencia del Fenómeno de Runge');
legend('Location', 'best');

%% =======================================================================
% B1.3: INTERPOLACIÓN EN f = 1000 Hz Y VALIDACIÓN LEAVE-ONE-OUT (LOO)
% =======================================================================
% Seleccionamos el polinomio de Grado 5 por ser el más estable globalmente
grado_seleccionado = 5; 
f_interp = 1000;
Z_interp_1000 = polyval(polyfit(f, Z, grado_seleccionado), f_interp);

% Fijamos una semilla aleatoria para que tus resultados sean reproducibles
rng(42); 
num_puntos_loo = 5;
indices_azar = randperm(length(f), num_puntos_loo);

% Inicialización de vectores para almacenar resultados del LOO
f_loo_puntos = zeros(num_puntos_loo, 1);
Z_reales     = zeros(num_puntos_loo, 1);
Z_predichos  = zeros(num_puntos_loo, 1);
err_relativos = zeros(num_puntos_loo, 1);

fprintf('========================================================================\n');
fprintf('   VALIDACIÓN LEAVE-ONE-OUT (LOO) PARA POLINOMIO DE GRADO %d\n', grado_seleccionado);
fprintf('========================================================================\n');
fprintf('%-10s %-12s %-15s %-15s %-15s\n', 'Iteración', 'Freq (Hz)', 'Z Real (\Omega)', 'Z Pred (\Omega)', 'Error Rel.');

for i = 1:num_puntos_loo
    idx_omitido = indices_azar(i);

    % Remover el punto i-ésimo del set de entrenamiento
    f_train = f; f_train(idx_omitido) = [];
    Z_train = Z; Z_train(idx_omitido) = [];

    % Re-entrenar el polinomio sin el punto omitido
    p_loo = polyfit(f_train, Z_train, grado_seleccionado);

    % Evaluar en el punto omitido
    f_loo_puntos(i) = f(idx_omitido);
    Z_reales(i)     = Z(idx_omitido);
    Z_predichos(i)  = polyval(p_loo, f_loo_puntos(i));

    % Calcular error relativo absoluto
    err_relativos(i) = abs(Z_predichos(i) - Z_reales(i)) / Z_reales(i);

    fprintf('%-10d %-12.4f %-15.4f %-15.4f %-15.4f\n', ...
        i, f_loo_puntos(i), Z_reales(i), Z_predichos(i), err_relativos(i));
end

fprintf('------------------------------------------------------------------------\n');
fprintf('Valor interpolado final en f = 1000 Hz: %.4f \Omega\n', Z_interp_1000);
fprintf('Error relativo promedio estimado (LOO): %.4f\n', mean(err_relativos));
fprintf('========================================================================\n');