% PARTE C: DERIVACIÓN NUMÉRICA MEDIANTE SPLINE ANALÍTICO 
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

%% C.1: RECONSTRUCCIÓN DEL SPLINE NATURAL (Cálculo de M = f'')
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

%% C.2: CÁLCULO DE LA PRIMERA DERIVADA EN LOS NODOS ORIGINALES
dZ_df_nodos = zeros(N, 1);
for i = 1:N
    k = i; 
    if i == N, k = N - 1; end 
    h = f(k+1) - f(k);
    x = f(i);
    dZ_df_nodos(i) = -(M(k) / (2 * h)) * (f(k+1) - x)^2 + ...
                      (M(k+1) / (2 * h)) * (x - f(k))^2 - ...
                      (Z(k)/h - (M(k) * h)/6) + ...
                      (Z(k+1)/h - (M(k+1) * h)/6);
end

%% C.3: EVALUACIÓN VECTORIZADA POR SEGMENTOS (¡Súper Rápido!)
f_fina = linspace(100, 2730, 100000); 
dZ_df_fina = zeros(size(f_fina));

for k = 1:(N-1)
    % Filtramos los puntos que pertenecen al segmento k actual
    if k == N-1
        idx_seg = (f_fina >= f(k) & f_fina <= f(k+1));
    else
        idx_seg = (f_fina >= f(k) & f_fina < f(k+1));
    end
    
    x_val = f_fina(idx_seg);
    h = f(k+1) - f(k);
    
    % Operación matemática vectorizada sobre el segmento
    dZ_df_fina(idx_seg) = -(M(k) / (2 * h)) * (f(k+1) - x_val).^2 + ...
                           (M(k+1) / (2 * h)) * (x_val - f(k)).^2 - ...
                           (Z(k)/h - (M(k) * h)/6) + ...
                           (Z(k+1)/h - (M(k+1) * h)/6);
end

% Localizar el cruce por cero (mínimo de impedancia)
idx_min = find(dZ_df_fina(1:end-1) < 0 & dZ_df_fina(2:end) > 0, 1);
f_minimo = f_fina(idx_min);

% Calcular impedancia exacta en el mínimo hallado
k_m = find(f <= f_minimo, 1, 'last'); if k_m == N, k_m = N-1; end
h_m = f(k_m+1) - f(k_m);
Z_minimo = (M(k_m)/(6*h_m))*(f(k_m+1)-f_minimo)^3 + (M(k_m+1)/(6*h_m))*(f_minimo-f(k_m))^3 + ...
           (Z(k_m)/h_m - (M(k_m)*h_m)/6)*(f(k_m+1)-f_minimo) + (Z(k_m+1)/h_m - (M(k_m+1)*h_m)/6)*(f_minimo-f(k_m));

% Segunda derivada analítica en el mínimo
d2Z_df2_minimo = (M(k_m) / h_m) * (f(k_m+1) - f_minimo) + (M(k_m+1) / h_m) * (f_minimo - f(k_m));

%% C.4: GRÁFICA DE LA PRIMERA DERIVADA
figure('Name', 'Primera Derivada de la Impedancia', 'NumberTitle', 'off');
plot(f_fina, dZ_df_fina, 'g-', 'LineWidth', 2, 'DisplayName', 'd|Z|/df (Spline Analítico)'); hold on;
plot(f, dZ_df_nodos, 'ko', 'MarkerFaceColor', 'r', 'DisplayName', 'Derivada en Nodos');
plot(f_minimo, 0, 'bx', 'MarkerSize', 12, 'LineWidth', 2.5, 'DisplayName', sprintf('Mínimo en %.2f Hz', f_minimo));
plot([100, 2730], [0, 0], 'k--', 'LineWidth', 0.8); 
grid on;
xlabel('Frecuencia f (Hz)', 'FontWeight', 'bold');
ylabel('d|Z|/df (\Omega/Hz)', 'FontWeight', 'bold');
title('Derivada Primera Analítica de la Impedancia Bioeléctrica');
legend('Location', 'best');
drawnow; % Fuerza a MATLAB a renderizar los gráficos en pantalla de inmediato

%% MOSTRAR REPORTES EN CONSOLA
fprintf('========================================================================\n');
fprintf('         DERIVACIÓN NUMÉRICA \n');
fprintf('========================================================================\n');
fprintf('Frecuencia exacta del mínimo (Z'' = 0)   : %.4f Hz\n', f_minimo);
fprintf('Impedancia mínima alcanzada en ese punto : %.4f ohmios\n', Z_minimo);
fprintf('Segunda derivada d^2|Z|/df^2 en el mínimo: %.8f ohmios/Hz^2\n', d2Z_df2_minimo);
if d2Z_df2_minimo > 0
    fprintf('Signo de la segunda derivada             : POSITIVO (+)\n');
    fprintf('Conclusión de estabilidad                : Mínimo Local Estable (Cóncavo hacia arriba)\n');
else
    fprintf('Signo de la segunda derivada             : NEGATIVO o CERO (-)\n');
    fprintf('Conclusión de estabilidad                : Inestable\n');
end
fprintf('========================================================================\n');