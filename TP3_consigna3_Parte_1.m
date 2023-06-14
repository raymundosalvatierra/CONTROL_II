%CARRITO GRUA - CONTROLADOR CON OBSERVADOR
% Trabajo Practico 3 - Item 1 - 2023
% Alumno: SALVATIERRA HERMAN RAYMUNDO
%clc

%Variables
m=.1;
Fricc=0.1;
long=1.6;
g=9.8; 
M=1.5;

%TIEMPOS de muestreo, simulacion, de Euler e integracion
Ts=0.01;
KMAX=1500-1;
Veces_Euler=100;
h=1e-4;  %tamaño de paso
t_d=(0:KMAX)*Ts;
colorc='b';

%referencias para posicion y angulo
ref1=10;
ref2=0;

%condiciones iniciales
phi(1)=pi;color='-r';color_='r';

% Mat_Ac=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0];
% Mat_Bc=[0; 1/M; 0; -1/(long*M)];
% Mat_C=[1 0 0 0;0 0 1 0];
% I=eye(4);


% %Version linealizada en el equilibrio inestable
Mat_Ac=[0 1 0 0;0 -Fricc/M -m*g/M 0;0 0 0 1; 0 -Fricc/(long*M) -g*(m+M)/(long*M) 0];
Mat_Bc=[0;1/M;0;1/(long*M)];
Mat_C=[1 0 0 0; 0 0 1 0]; %2 salidas porque mido la altura y el angulo phi
I=eye(4);

H=[0;0;0;0];d_tao=Ts/100;tao=0;
for hh=1:100
 dH=expm(Mat_Ac*tao)*Mat_Bc*d_tao;
 H=H+dH;
 tao=tao+d_tao;
end

Mat_B=H;
Mat_A=expm(Mat_Ac*Ts);

Mat_B=H;
Mat_A=expm(Mat_Ac*Ts);
Mat_Aa=[Mat_A,zeros(4,2);-Mat_C*Mat_A, eye(2)];
Mat_Ba=[Mat_B;-Mat_C*Mat_B];
Mat_Ma=[Mat_Ba Mat_Aa*Mat_Ba Mat_Aa^2*Mat_Ba Mat_Aa^3*Mat_Ba Mat_Aa^4*Mat_Ba];%Matriz Controlabilidad
rango=rank(Mat_Ma);

%Calculo del DLQR
Qc=diag([1000 10 100 10 1e-1 40]);
R=1; %Ts=0.01

%Construcción del Hamiltoniano para el cálculo del controlador
H=inv([eye(6) Mat_Ba*inv(R)*Mat_Ba'; zeros(6) Mat_Aa'])*[Mat_Aa zeros(6);-Qc eye(6)];
[V,D]=eig(H);MX1X2=[];
for ii=1:12
 if abs(D(ii,ii))<1
 MX1X2=[MX1X2 V(:,ii)];
 end
end
MX1=MX1X2(1:6,:); MX2=MX1X2(7:12,:);

Pc=real(MX2*pinv(MX1)); % [K1,P,E]=dlqr(Mat_A,Mat_B,Q,R); %En Octave
Ka=inv(R+Mat_Ba'*Pc*Mat_Ba)*Mat_Ba'*Pc*Mat_Aa;
K=Ka(1:4); KI=-Ka(5);    %revisar signo de esto
aut_controlador=abs(eig(Mat_Aa-Mat_Ba*Ka));%verifico polos con parte real negativa
%Fin cálculo del controlador-

t=0;x=[0;0;phi(1);0];
delta(1)=x(1);delta_p(1)=x(2);phi(1)=x(3);phi_p(1)=x(4);
ve1(1)=0;ve2(1)=0;
delta_(1)=0;delta_p_(1)=0;phi_(1)=0;phi_p_(1)=0;u_k(1)=0;

V_L=[x;0;0]'*Pc*[x;0;0];Jl=0;Jmin_L=V_L;%e1=0;e2=0;
ki=1;

%Verificacion de la solucion con el modelo no lineal en tiempo continuo
T=t(ki);x=[0;0;phi(1);0];
delta=x(1);delta_p=x(2);phi=x(3);phi_p=x(4);phi_pp(1)=0;delta_pp(1)=0;
u=[];
y=Mat_C*x;ve1(1)=0;ve2(1)=0;phi_(1)=phi(1);

Jn=0;V_NL=[x;0;0]'*Pc*[x;0;0];Jmin_NL=V_NL;i=1;

for ki=2:KMAX
    
    V_NL=[V_NL [x;ve1(ki-1);ve2(ki-1)]'*Pc*[x;ve1(ki-1);ve2(ki-1)]]; %funcion de Lyapunov
    Y_=Mat_C*x; %se mide aca
    
   ve1(ki)=ve1(ki-1) +ref1-Y_(1);
   ve2(ki)=ve2(ki-1) +ref2-Y_(2);
   
   u1(ki)=-Ka*[x;ve1(ki);ve2(ki)]; %sin observador
   
%    if abs(u1(ki))<10
%        u1(ki)=0;
%    end
   
   if (ki==((KMAX+1)/2))   
        m=10*m
        ref1=0;
   end
   
   for kii=1:Veces_Euler
    
    u(i)=u1(ki);   
       
    delta_pp=(1/(M+m))*(u(i)-m*long*phi_pp*cos(phi(i))+ m*long*phi_p(i)^2*sin(phi(i))- Fricc*delta_p(i));
    phi_pp=(1/long)*(g*sin(phi(i)) - delta_pp*cos(phi(i)));
    delta_p(i+1) = delta_p(i) + h*delta_pp;
    delta(i+1) =  delta(i) + h*delta_p(i);
    phi_p(i+1) = phi_p(i) + h*phi_pp;
    phi(i+1) = phi(i) + h*phi_p(i);
    phi_(i)=x(3);
    i=i+1;
   end
   
   x=[delta(i-1);delta_p(i-1);phi(i-1);phi_p(i-1)]; %Aca esta x(k+1)
   Jn=[Jn Jn(ki-1)+[x;ve1(ki);ve2(ki)]'*Qc*[x;ve1(ki);ve2(ki)]+u1(ki)'*R*u1(ki)];%incremental de costos
end

u(i)=u1(ki);t=(1:i)*h;

figure(1);hold on;
subplot(3,2,1);plot(t,phi,color);grid on;title('phi_t');hold on;
subplot(3,2,2);plot(t,phi_p,color);title('phi_p');grid on;hold on;
subplot(3,2,3);plot(t,delta,color);grid on;title('delta_t');hold on;
subplot(3,2,4);plot(t,delta_p,color);title('delta_p');grid on;hold on;
subplot(3,1,3);plot(t,u,color);grid on;title('Accion de Control');xlabel('Tiempo en Seg');hold on;

figure(2);hold on;
subplot(2,2,1);plot(phi,phi_p,color);grid on;
xlabel('phi_t');hold on;
ylabel('phi_p');hold on;
subplot(2,2,2);plot(delta,delta_p,color);grid on;xlabel('Posicion carro');hold on;
xlabel('delta_t');hold on;
ylabel('delta_p');hold on;

subplot(2,2,3);
semilogy(t_d(1:end-1),Jn,color);grid on;xlabel('Modelo no lineal');
xlabel('Tiempo en Seg');ylabel('Funcionales en J y V');hold on;
semilogy(t_d(1:end-1),Jmin_NL*ones(size(Jn)),colorc);
semilogy(t_d(1:end-1),V_NL,colorc);
set(gca,'fontsize');

% figure(4)
% plot(t,ua,color);title('Accion de control');hold on;
% figure(5)
% plot(t,Tll,color);title('Torque');hold on;
%