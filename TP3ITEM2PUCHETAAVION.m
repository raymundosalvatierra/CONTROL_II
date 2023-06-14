%TP 3 ÍTEM 2
%Profesor: Pucheta
%Alumno: SALVATIERRA HERMAN RAYMUNDO
%CON OBSERVADOR y DLQR 

clear all;clc
%Constantes
w=9;
a=0.07;
b=5;
c=150; % [m/seg]
delta_t=1e-3;
u=0;
Ts = 10e-3; KMAX=70/Ts; % es el tiempo de la simulacion

%referencias
ref1 = -100; % el primer caso era con -100
ref2=0;

alpha_0=0;
fi_0=0;
fip_0=0;
h_0 = 500; % altura inicial, el primer caso era con 500


%Planteo las matrices
Mat_A=[-a a 0 0 ; 0 0 1 0; w^2 -(w^2) 0 0; c 0 0 0];
Mat_B=[0;0;b*w^2;0];
Mat_C=[0 0 0 1; 0 1 0 0]; %2 salidas por que mido la altura y el angulo phi];
MAT_D=0;
% [n,d] = ss2tf(Mat_A,Mat_B,Mat_C,MAT_D)
I=eye(4);
%paso al sistema discreto
H=[0;0;0;0];d_tao=Ts/100;tao=0;
for hh=1:100
 dH=expm(Mat_A*tao)*Mat_B*d_tao;
 H=H+dH;
 tao=tao+d_tao;
end
Mat_B=H;
Mat_A=expm(Mat_A*Ts);

% Construccion del sistema ampliado
Mat_Aa=[Mat_A,zeros(4,2);-Mat_C*Mat_A, eye(2)];
Mat_Ba=[Mat_B;-Mat_C*Mat_B];
Mat_Ma=[Mat_Ba Mat_Aa*Mat_Ba Mat_Aa^2*Mat_Ba Mat_Aa^3*Mat_Ba Mat_Aa^4*Mat_Ba];%Matriz Controlabilidad
rango=rank(Mat_Ma);

%calculo DLQR
Qc=diag([1000 1000000 1000 1/10000 1/600 2000]); R=2000000000; 

%Contruccion del Hamiltoniano para el calculo del controlador
  H=inv([eye(6) Mat_Ba*inv(R)*Mat_Ba'; zeros(6) Mat_Aa'])*[Mat_Aa zeros(6);-Qc eye(6)];
[V,D]=eig(H);MX1X2=[];
for ii=1:12
 if abs(D(ii,ii))<1
 MX1X2=[MX1X2 V(:,ii)];
 end
end
MX1=MX1X2(1:6,:); MX2=MX1X2(7:12,:);
Pc=real(MX2*pinv(MX1)); 
Ka=inv(R+Mat_Ba'*Pc*Mat_Ba)*Mat_Ba'*Pc*Mat_Aa;
K=Ka(1:4);KI=-Ka(5);
aut_controlador=abs(eig(Mat_Aa-Mat_Ba*Ka));

%Calculo del Observador de estados
Mat_AO=Mat_A';
Mat_BO=Mat_C';
Mat_CO=Mat_B';
Mat_Qobs=[Mat_C;Mat_C*Mat_A;Mat_C*Mat_A^2;Mat_C*Mat_A^3];
rango_matriz_obs=rank(Mat_Qobs);

Qobs=diag([0.1 1000000 1000 1/10]);Ro=diag([900000000 100000000]);

%Contruccion del Hamiltoniano para el calculo del Observador
Ho=inv([eye(4) Mat_BO*inv(Ro)*Mat_BO'; zeros(4) Mat_AO'])*[Mat_AO zeros(4);-Qobs eye(4)];

[Vo,Do]=eig(Ho);MX1X2=[];

for ii=1:8
 if abs(Do(ii,ii))<1
 MX1X2=[MX1X2 Vo(:,ii)];
 end
end
MX1o=MX1X2(1:4,:); MX2o=MX1X2(5:8,:);
Po=real(MX2o*inv(MX1o)); 
Kobs=(inv(Ro+Mat_BO'*Po*Mat_BO)*Mat_BO'*Po*Mat_AO)';
p_observador=abs(eig(Mat_A-Kobs*Mat_C)); %Verifica polos de observabilidad

ki=1;

%Verificacion con el MODELO LINEAL
t=0; x=[0;0;0;h_0];x_hat=[0;0;0;0];
ve1(1)=0;ve2(1)=0;
alpha(1)=x(1); phi(1)=x(2); phip(1)=x(3); h(1)=x(4);

V_L=[x;0;0]'*Pc*[x;0;0];
x_hat=[0;0;0;0];
Jl=0;Jmin_L=V_L; %e1=0 ; e2=0;

for ki=2:KMAX
 t=[t (ki-1)*Ts];
 
 Y_=Mat_C*x; %Se mide ACA
 
 e1=ref1-Y_(1);
 e2=ref2-Y_(2);
 
 V_L=[V_L [x;ve1(ki-1);ve2(ki-1)]'*Pc*[x;ve1(ki-1);ve2(ki-1)]];%funcion de Liapunov
 
 ve1(ki)=ve1(ki-1)+e1;
 ve2(ki)=ve2(ki-1)+e2;
 
 % u=-Ka*[x_hat;ve1(ki);ve2(ki)];color='g';%con observador
  u=-Ka*[x;ve1(ki);ve2(ki)];color='r';%con conrtolador sin observador

 % Agrego la no linealidad de la zona muerta de 0.1 
 % if abs(u)<0.7
 %  u=0;
 %   end
 
 ys=Mat_C*x; %Se mide ACA
 x=Mat_A*x+Mat_B*u;
 
 Jl=[Jl Jl(ki-1)+[x;ve1(ki);ve2(ki)]'*Qc*[x;ve1(ki);ve2(ki)]+u'*R*u];
 %para el observador
 y_hat=Mat_C*x_hat;
 x_hat=Mat_A*x_hat+Mat_B*u+Kobs*(Y_-y_hat);%Se actualiza aca¡
 
 alpha(ki)=x(1);
 phi(ki)=x(2);
 phip(ki)=x(3);
 h(ki)=x(4); %Valores Observados
%  alpha_(ki)=x(1);
%  phi_(ki)=x(2);
%  phip_(ki)=x(3);
%  h_(ki)=x(4);
 u_k(ki)=u;
end
Jl=[Jl Jl(ki-1)+[x;ve1(ki);ve2(ki)]'*Qc*[x;ve1(ki);ve2(ki)]+u'*R*u];
V_L=[V_L [x;ve1(ki-1);ve2(ki-1)]'*Pc*[x;ve1(ki-1);ve2(ki-1)]];
u=u_k;

figure(1);hold on;
subplot(3,2,1);plot(t,alpha,color);title('\alpha');hold on;
subplot(3,2,2);plot(t,phi,color);title('\phi');hold on;
subplot(3,2,3); plot(t,phip,color);title('\phi_p');hold on;
subplot(3,2,4);plot(t,h,color);title('h');hold on;
subplot(3,1,3);plot(t,u,color);title('Accion de control');hold on; 

figure(2);hold on;
plot(alpha,h,color);grid on;xlabel('Angulo \alpha');ylabel('altura');hold on;