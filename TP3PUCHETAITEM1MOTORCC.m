%TP 3 CASO DE ESTUDIO 1
%ALUMNO: SALVATIERRA HERMAN RAYMUNDO
%motor CC - con Observador - zona muerta
clear all;clc ;
%variables
Laa=0.56e-3; 
J=0.0019;
Ra=1.35;
B=0.0082;         
Ki=0.1;
Km=0.1;
delta_t=1e-5;   %origianl 1e-6
u=12;           %ESCALON
Ts=0.00002;      %Tiempo de muestreo 20 us
KMAX=1/0.00001; 

%Torque
Tl=1.5;%va a ir variando para pi/2 y -pi/2

%Matrices
Mat_A=[-Ra/Laa -Km/Laa 0;Ki/J -B/J 0 ;0 1 0 ];
Mat_B=[1/Laa;0;0];
Mat_C=[0 0 1]; %x1=ia   x2=wr   x3=theta
Mat_D=[0];
I=eye(4);

%paso al sistema al tiempo discreto
H=[0;0;0] ; d_tao=Ts/100 ; tao=0;
for hh=1:100
    dH=expm(Mat_A*tao)*Mat_B*d_tao;
    H=H+dH;
    tao=tao+d_tao;
end
Mat_B=H;
Mat_A=expm(Mat_A*Ts);

% Construcción de las matrices ampliadas
Mat_Aa=[Mat_A,zeros(3,1);-Mat_C*Mat_A, eye(1)];
Mat_Ba=[Mat_B;-Mat_C*Mat_B];
Mat_Cc=[Mat_C 0];
Mat_Ma=[Mat_Ba Mat_Aa*Mat_Ba Mat_Aa^2*Mat_Ba Mat_Aa^3*Mat_Ba Mat_Aa^4*Mat_Ba];%Matriz Controlabilidad
rango=rank(Mat_Ma);

%calculo DLQR

Qc=diag([50 1/10000 10 1]); R=2; %conforme se disminuye R se aeclera el sistema, mejora

%Contrucción del Hamiltoniano para el cálculo del controlador
H=inv([eye(4) Mat_Ba*inv(R)*Mat_Ba'; zeros(4) Mat_Aa'])*[Mat_Aa zeros(4);-Qc eye(4)];
[V,D]=eig(H);MX1X2=[];
for ii=1:8
    if abs(D(ii,ii))<1
        MX1X2=[MX1X2 V(:,ii)];
    end
end
MX1=MX1X2(1:4,:); MX2=MX1X2(5:8,:);
Pc=real(MX2*inv(MX1)); 
Ka=inv(R+Mat_Ba'*Pc*Mat_Ba)*Mat_Ba'*Pc*Mat_Aa;
K=Ka(1:3);KI=-Ka(4);
aut_controlador=abs(eig(Mat_Aa-Mat_Ba*Ka));

%Cálculo del Observador de estados
Mat_AO=Mat_A';
Mat_BO=Mat_C';
Mat_CO=Mat_B';
Mat_Qobs=[Mat_C;Mat_C*Mat_A;Mat_C*Mat_A^2 ; Mat_C*Mat_A^3];
rango_matriz_obs=rank(Mat_Qobs);

Qobs=diag([1 .10 .09]); Ro=0.01; %inicialmente R era 100 para el controlador 

%Construcción del Hamiltoniano para el cálculo del Observador
Ho=[Mat_AO+Mat_BO*inv(Ro)*Mat_BO'*inv(Mat_AO')*Qobs -Mat_BO*inv(Ro)*Mat_BO'*inv(Mat_AO'); -inv(Mat_AO')*Qobs inv(Mat_AO')];
[Vo,Do]=eig(Ho);MX1X2=[];
for ii=1:6
    if abs(Do(ii,ii))<1
        MX1X2=[MX1X2 Vo(:,ii)];
    end
end
MX1o=MX1X2(1:3,:); MX2o=MX1X2(4:6,:);
Po=real(MX2o*inv(MX1o)); 
Kobs=(inv(Ro+Mat_BO'*Po*Mat_BO)*Mat_BO'*Po*Mat_AO)';
p_observador=abs(eig(Mat_A-Kobs*Mat_C)); %Verifica polos de observabilidad

ki=1;

J_(1)=0; V_(1)=0; 
psi(1)=0;
x_hat=[0;0;0];%inicializo el observador

delta_t=1e-4;   
tiempo=10;
pasos=round(tiempo/delta_t);
t=0:delta_t:(tiempo-delta_t);

%Hago la referencia:

ref1=(pi/2)*square(2*pi*t/10); 

ref2=0;
%Torque 
Tll=(Tl/2)+(Tl/2)*square(2*pi*t/10);
%condiciones iniciales
x=[0;0;0;0];

ua(1)=12;
u=12;
y=Mat_C*x(1:3);
ve1(1)=0  ;   ve2(1)=0   ;   y_hat=0;
Jn=0  ;  V_NL=[x(1:3);0]'*Pc*[x(1:3);0]   ;   Jmin_NL=V_NL   ;   i=1;

for ki=2:KMAX+1
    estado=[x(1,ki-1);x(2,ki-1);x(3,ki-1);x(4,ki-1)];%guardo el estado
    integracion=x(4,ki-1)+(ref1(1,ki-1)-Mat_Cc*estado);
    
        Y_  =   Mat_C*estado(1:3);
     y_hat  =   Mat_C*x_hat;                  %para el observador
    
    u_actual=-K*estado(1:3)+integracion*KI;color='r'  ;  %sin observador
  % u_actual=-K*x_hat+integracion*KI;color='g';          %Con observador
      
      %Agrego la no linealidad de la zona muerta de +-1
     if abs(u_actual)<20.5
         
         u_actual=0;
     end
      ua = [ua u_actual];
      
    %Ecuaciones del motor:
    ia_p=(-Ra/Laa)*estado(1)-(Km/Laa)*estado(2)+(1/Laa)*u_actual;
    w_p=(Ki/J)*estado(1)-(B/J)*estado(2)-(1/J)*Tll(ki-1);
    theta_p=estado(2);
    
    xp_actual=[ia_p; w_p; theta_p];
    
    xsig=estado(1:3)+delta_t*xp_actual;
    x(1,ki)=xsig(1);
    x(2,ki)=xsig(2);
    x(3,ki)=xsig(3);
    x(4,ki)=integracion;
   
    %___OBSERVADOR___
    x_hat=Mat_A*x_hat+Mat_B*u_actual+Kobs*(Y_-y_hat);
    y_sal_O(i)=Mat_C*x_hat;
    y_sal(i)=Mat_Cc*estado;
    i = i+1;
end

u(i)=u_actual;t=(1:i)*delta_t;
figure(1);hold on;
%plot(t,x(3,:),color);title('Angulo Theta y referencia');hold on;plot(t,ref1);hold on;

subplot(3,2,1);plot(t,x(3,:),color);title('Angulo Theta y referencia');hold on;plot(t,ref1);hold on;
subplot(3,2,2);plot(t,x(2,:),color);title('\Omega_r');hold on
subplot(3,2,3); plot(t,x(1,:),color);title('Corriente Ia');hold on;
subplot(3,2,4); plot(t,Tll,color);title('Torque');hold on;
subplot(3,1,3);plot(t,ua,color);title('Accion de control');hold on;