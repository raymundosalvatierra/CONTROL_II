figure(1);
%primero grafico al controlador inicial en color azul
plot(t,x(3,1:end),'b');title('Angulo Theta y referencia');grid on;hold on;
plot(t,ref);legend('Angulo Theta','ref');

%despues grafico al observador del item 2 en color rojo
plot(t,x_o(3,:),'r');title('Angulo Theta y referencia');hold on;
plot(t,ref);