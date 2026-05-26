%% PARTE 2: DERIVACIÓN NUMÉRICA (VERSIÓN NATIVA COMPATIBLE)
clear; clc; close all;
format longG;

% Datos Originales del Ensayo (image_449763.png)
f = [10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, ...
    60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5];

V = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, 1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041, ...
    -0.057, -0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809, 0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004];

h = 2.5; 
idx = @(val) find(abs(f - val) < 1e-4);

%% 1 y 2. Diferencias Finitas
% Frecuencia: 40.0 kHz
i40 = idx(40.0);
dV_40_c2 = (V(i40+1) - V(i40-1)) / (2*h);
dV_40_c4 = (-V(i40+2) + 8*V(i40+1) - 8*V(i40-1) + V(i40-2)) / (12*h);

% Frecuencia: 70.0 kHz
i70 = idx(70.0);
dV_70_c2 = (V(i70+1) - V(i70-1)) / (2*h);
dV_70_c4 = (-V(i70+2) + 8*V(i70+1) - 8*V(i70-1) + V(i70-2)) / (12*h);

% Frecuencia: 100.0 kHz
i100 = idx(100.0);
dV_100_c2 = (V(i100+1) - V(i100-1)) / (2*h);
dV_100_c4 = (-V(i100+2) + 8*V(i100+1) - 8*V(i100-1) + V(i100-2)) / (12*h);

% Extremo inferior: 10.0 kHz (Fórmula progresiva de orden 2)
i10 = idx(10.0);
dV_10_p2 = (-3*V(i10) + 4*V(i10+1) - V(i10+2)) / (2*h);

%% 4. Derivada Analítica del Spline Cúbico (Estructura Nativa pp)
pp = spline(f, V); % Genera el spline cúbico de forma nativa
[breaks, coefs, l, k, d] = unmkpp(pp);

% Derivando los coeficientes de los polinomios cúbicos: [3*a, 2*b, c]
coefs_deriv = [3*coefs(:,1), 2*coefs(:,2), coefs(:,3)];
pp_deriv = mkpp(breaks, coefs_deriv);

% Evaluación de la derivada en los puntos de interés
dV_spline = fnval(pp_deriv, [10.0, 40.0, 70.0, 100.0]);

%% Impresión de Resultados en Consola
fprintf('========================================================================\n');
fprintf('                 REPORTE DE DERIVACIÓN NUMÉRICA              \n');
fprintf('========================================================================\n');
fprintf('f (kHz)   Dif. Centrada O2   Dif. Centrada O4   Prog. Orden 2   Spline Natural\n');
fprintf('------------------------------------------------------------------------\n');
fprintf(' 10.0           N/A                N/A          %11.6f   %14.6f\n', dV_10_p2, dV_spline(1));
fprintf(' 40.0       %11.6f        %11.6f          N/A          %14.6f\n', dV_40_c2, dV_40_c4, dV_spline(2));
fprintf(' 70.0       %11.6f        %11.6f          N/A          %14.6f\n', dV_70_c2, dV_70_c4, dV_spline(3));
fprintf('100.0       %11.6f        %11.6f          N/A          %14.6f\n', dV_100_c2, dV_100_c4, dV_spline(4));
fprintf('========================================================================\n');