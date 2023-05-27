%PLANTEO CASO 2: CARRITO CON PENDULO INVERTIDO
% TP 2 
% Item 3 : carro con pendulo invertido
% Alumno: Salvatierra Herman Raymundo
%Limpieza de Pantalla
clear all, close all, clc
%Variables
m = 0.1 ; F = 0.1; l = 1.6 ; g = 9.8 ; M = 1.5;
%Despues de hacer calculos matematicos puedo plantear mi modelo con las
%matrices de estado:
% Definir las matrices de estado conocidas
A = [0  1  0  0 ; 0  -(F/M)  -((m*g)/M)  0 ; 0 0 0 1; 0  (F/(M*l))  (((M+m)*g)/(M*l))  0]
B = [0 ; 1/M ; 0  ; (-1)/(M*l)]
C = [1   0   0   0]
D = [0]

%%Voy a plantear mi funcion de transferncia a partir de mis matrices de
%%estado:
[n,m] = ss2tf(A,B,C,D)
G = tf(n,m)
pole(G)

%Luego de conocer los polos de mi sistema puedo calcular el tamaño de paso
%(h) y el tiempo de simulacion minimo (T):
%Para el tamaño de paso usaré el menor polo en la fórmula y eso lo dividiré
%por 3 para obtener una mejor resolución:
h = ((log(0.95))/(-2.5582))/(3)

%Para calcular el tiempo de simlacion (T) utilizo el polo más grande en la
%formula y lo aproximo a un tiempo que sea razonable para mostrar el
%comportamiento de mi sistema:
T = (log(0.05))/(-0.0625)


%Quiero ver el lugar de raices de mi funcion de tranferencia:
%figure(1)
%rlocus(G)

%Ahora veré la respuesta al escalon que tiene mi sistema
%figure(2)
%step(G)
