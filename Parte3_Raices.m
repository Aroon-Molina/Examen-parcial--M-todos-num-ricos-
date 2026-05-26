%% PARTE 3: RAÍCES POR CAMBIO DE SIGNO Y BISECCIÓN
clear; clc; close all;
format longG;

% 1. Datos Originales del Ensayo (image_449763.png)
f = [10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, ...
    60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5];

V = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, 1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041, ...
    -0.057, -0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809, 0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004];

tol = 1e-6; % Tolerancia requerida para los métodos numéricos

%% 1, 2 y 3. Bisección Básica mediante Interpolación Lineal (Tabla Directa)
func_lin = @(x) interp1(f, V, x, 'linear');

r1_lineal = ejecutar_biseccion(func_lin, 55.0, 57.5, tol);
r2_lineal = ejecutar_biseccion(func_lin, 62.5, 65.0, tol);

%% 4. Refinamiento utilizando la curva continua del Spline Cúbico
func_spline = @(x) interp1(f, V, x, 'spline');

r1_spline = ejecutar_biseccion(func_spline, 55.0, 57.5, tol);
r2_spline = ejecutar_biseccion(func_spline, 62.5, 65.0, tol);

%% Despliegue de Resultados Oficiales en Consola
fprintf('========================================================================\n');
fprintf('             REPORTE DE LOCALIZACIÓN DE RAÍCES              \n');
fprintf('========================================================================\n');
fprintf('Cruce por Cero       Bisección Lineal       Spline Cúbico       Diferencia (kHz)\n');
fprintf('------------------------------------------------------------------------\n');
fprintf('Primer Cruce (r1)      %11.6f         %11.6f          %11.6f\n', r1_lineal, r1_spline, abs(r1_lineal - r1_spline));
fprintf('Segundo Cruce (r2)     %11.6f         %11.6f          %11.6f\n', r2_lineal, r2_spline, abs(r2_lineal - r2_spline));
fprintf('========================================================================\n');

%% Función Local: Algoritmo Estándar del Método de Bisección
function root = ejecutar_biseccion(func, a, b, tol)
if func(a) * func(b) > 0
    error('No existe cambio de signo en el intervalo seleccionado.');
end
while (b - a) / 2 > tol
    c = (a + b) / 2;
    if func(c) == 0
        break;
    end
    if func(a) * func(c) < 0
        b = c;
    else
        a = c;
    end
end
root = (a + b) / 2;
end