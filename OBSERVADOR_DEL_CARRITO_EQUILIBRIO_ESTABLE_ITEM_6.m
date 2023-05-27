%TP 2  - ITEM 6 - OBSERVADOR - CARRITO - EQUILIBRIO ESTABLE
%ALUMNO: HERMAN RAYMUNDO SALVATIERRA
%clear all, close all, clc
%Variables
m =0.1 ;   Fricc=0.1  ;    long=1.6  ; g=9.8    ; M=1.5;
h=0.005; tiempo=(80/h);  delta_pp=0  ; phi_pp=0 ; t=0:h:tiempo*h;
phi_p=0:h:tiempo*h    ;   phi=0:h:tiempo*h      ; delta=0:h:tiempo*h;
delta_p=0:h:tiempo*h  ; u=linspace(0,0,tiempo+1);
delta_pp_o=0;phi_pp_o=0; 
%Condiciones iniciales
phi(1)=pi    ; color='b'; % sigo usando las matrices con algunos elementos cambiados de signo
phi_p(1)=0   ; delta_p(1)=0    ; u(1)=0    ; delta(1)=0   ; i=1  ; indice=0;
phi_o(1)=pi  ; color_o='g'     ;
phi_p_o(1)=0 ; delta_p_o(1)=0  ;  uo(1)=0  ; delta_o(1)=0;
ref = 2 ;
xop=[0; 0; pi;0]; % aca tambien uso las matrices modificadas
%Matrices para version linealzada del equilibrio inestable
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0;0 0 0 1; 0 -Fricc/(long*M) -g*(m+M)/(long*M) 0]; 
Mat_B=[0;1/M ; 0; 1/(long*M)]; 
Mat_C=[1 0 0 0];% la salida es la posicion 
Mat_D=[0];
%Matrices ampliadas
Mat_Aa=[Mat_A zeros(4,1); -Mat_C 0];
Mat_Ba=[Mat_B; 0];
Mat_Cc=[Mat_C 0];
Q=diag([2 1.0 10 01 0.1]);
R=0.1; 

%Construccion del Hamiltoniano para el calculo del controlador
Ha=[Mat_Aa -Mat_Ba*inv(R)*Mat_Ba'; -Q -Mat_Aa'];
[n,va]=size(Ha);

[V,D]=eig(Ha);MX1X2=[];
for ii=1:n
    if real(D(ii,ii))<0
        MX1X2=[MX1X2 V(:,ii)];
    end
end

MX1=MX1X2(1:n/2,:);
MX2=MX1X2(n/2+1:end,:);

P=real(MX2*inv(MX1));%Tomo la parte real por un tema de residuos
Ka=inv(R)*Mat_Ba'*P;
K=Ka(1:4); 
KI=-Ka(5);
eig(Mat_Aa-Mat_Ba*Ka)%Verifico polos con parte real negativa
% break
%Fin cálculo del controlador
J_(1)=0
V_(1)=0;
psi(1)=0;

%Calculo del observador

Mat_A_d = Mat_A';
Mat_B_d= Mat_C';
Mat_C_d=Mat_B';

Qo=diag([1000000   1    1    1/1000]); 
Ro=100;

%Construccion del Hamiltoniano para el calculo del observador
Ho=[Mat_A_d -Mat_B_d*inv(Ro)*Mat_B_d'; -Qo -Mat_A_d'];
[no,va]=size(Ho);

[V,D]=eig(Ho);MX1X2=[];
for ii=1:no
    if real(D(ii,ii))<0
        MX1X2=[MX1X2 V(:,ii)];
    end
end

MX1=MX1X2(1:no/2,:);
MX2=MX1X2(no/2+1:end,:);

Po=real(MX2*inv(MX1));%Tomo la parte real por un tema de residuos
Ko=(inv(Ro)*Mat_B_d'*Po)';
eig(Mat_A - Ko*Mat_C);

%Fin del observador

J_(1)=0; V_(1)=0;
psi(1)=0;
x_hat=[0;0;0;0];

while(i<(tiempo+1))
    estado=[delta(i);delta_p(i);phi(i);phi_p(i)];
    psi_p=ref-Mat_C*(estado - xop);
    psi(i+1)=psi(i) + psi_p*h;
    u(i) = -K*(x_hat - xop) + KI*psi(i+1);
    
    if(i == 20/h)
        m = 10*m;
        ref=0;
    end
    
    delta_pp=(1/(M+m))*(u(i)-m*long*phi_pp*cos(phi(i))+ m*long*phi_p(i)^2*sin(phi(i))- Fricc*delta_p(i));
    phi_pp=(1/long)*(g*sin(phi(i)) - delta_pp*cos(phi(i)));
    delta_p(i+1) = delta_p(i) + h*delta_pp;
    delta(i+1) =  delta(i) + h*delta_p(i);
    phi_p(i+1) = phi_p(i) + h*phi_pp;
    phi(i+1) = phi(i) + h*phi_p(i);
    
   %OBSERVADOR 
    y_sal_O(i) = Mat_C*x_hat;
    y_sal(i) = Mat_C*estado;
    x_hatp = Mat_A*x_hat + Mat_B*u(i) + Ko*(y_sal(i) - y_sal_O(i));
    x_hat = x_hat+ h*x_hatp;
    
    i= i+1;
end

figure(1); hold on;
%subplot(3,2,1);plot(t,phi_p,color);grid on; title('Velocidad angulo');hold on;
subplot(3,2,2);plot(t,phi,color);grid on; title('Angulo phi');hold on;
subplot(3,2,3);plot(t,delta,color);grid on; title('Delta (posicion carro)');hold on;

%subplot(3,2,4);plot(t,delta_p,color);grid on; title('Velocidad carro');hold on;
subplot(3,1,3);plot(t,u,color);grid on; title('Accion de control, u');xlabel('Tiempo en seg');hold on;

%figure(2);hold on;subplot(2,2,1);plot(phi,phi_p,color);grid on; xlabel('Angulo');
%ylabel('Velocidad angular'); hold on; subplot(2,2,2); plot(delta,delta_p,color);grid on;
%xlabel('Posicion carro');ylabel('Velocidad carro');hold on;