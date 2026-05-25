%Definición de los vectores de las mediciones
%Frecuencia
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, ...
    460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, ...
    1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];

%Impedancia
Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, ...
    132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, ...
    135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

figure('Name', 'Análisis Exploratorio', 'NumberTitle', 'off');
plot(f, Z, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 5);
xlabel('Frecuencia f (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Magnitud de Impedancia |Z| (\Omega)', 'FontSize', 11, 'FontWeight', 'bold');
title('Impedancia vs Frecuencia', 'FontSize', 12);
grid on;

%Búsqueda del mínimo en la muestra
[min_Z, idx_min] = min(Z);
f_min = f(idx_min);

%resultados con 4 decimales
fprintf('====================================================\n');
fprintf('RESULTADOS DEL ANÁLISIS EXPLORATORIO\n');
fprintf('====================================================\n');
fprintf('Frecuencia del mínimo estimado: %10.4f Hz\n', f_min);
fprintf('Impedancia mínima observada:    %10.4f ohm\n', min_Z);
fprintf('====================================================\n');