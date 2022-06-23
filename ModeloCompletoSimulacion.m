% % % % Modelo de completo para una señal de entrada  % % % % 

% Para la realización de pruebas se simula la señal obtenida, para ello se
% crea como señal de entrada una senoidal con una amplitud de 20, una fase 
% de 5 rad/s y un desfase de 1 rad, inicialmente libre de ruido.

% % Parametros de la simulación:
t_sim=2;
h=0.001;    %fs=1kHz

t=0:h:t_sim;

% % Parametros de la entrada:
A=20;
w=5;
phi=1;
K=10;

y_valores=A*sin(w*t+phi)+K;

% % Parametros del ruido:
% A_ruido=0.01;
% ruido=-A_ruido+(2*A_ruido)*rand(1,length(t));

% y_valores=y_valores+ruido;


% % Configuración del centinela.
rango=120; 
condicion=0.0097;

% % Configuración del límite de tiempo y de omisión de datos.
tiempo_critico=2;
margen=0;

t=0:h:tiempo_critico;
% % Parametros del filtro paso bajo:
fc=10;

% % Modelo teorico de la frecuencia:
tic;

n=zeros(1,length(t));
d=zeros(1,length(t));
f=zeros(1,length(t));

resultados=zeros(1,length(t));

tic;
pos=2;
tolerancia=Inf();


y_valores=lowpass(y_valores,fc,1/h);

% Numerador
n1=(t(pos).^3).*y_valores(pos);
var=y_valores(1:pos);
n2=-3*(t(pos).^2).*trapz(t(1:pos),var);
var=var.*t(1:pos);
n3=+24*t(pos).*trapz(t(1:pos),var);
var=var.*t(1:pos);
n4=-30.*trapz(t(1:pos),var);
n(pos)=abs(n1+n2+n3+n4);

% Denominador
d1=1.5*(t(pos).^2).*trapz(t(1:pos),var);
var=var.*t(1:pos);
d2=-4*t(pos).*trapz(t(1:pos),var);
var=var.*t(1:pos);
d3=2.5*trapz(t(1:pos),var);
d(pos)=abs(d1+d2+d3);

% Frecuencia
f1=trapz(t(1:pos),n(1:pos));
f2=trapz(t(1:pos),d(1:pos));

f(pos)=sqrt(f1/f2);      
resultados(pos)=tolerancia;

pos=pos+1;

while( pos<=length(t)  &&  tolerancia>=condicion)
    % Numerador
    n1=(t(pos).^3).*y_valores(pos);
    var=y_valores(pos-1:pos);
    n2=(t(pos).^2).*(n2./(t(pos-1).^2)-3*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    n3=t(pos).*(n3./t(pos-1)+24*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    n4=n4-30.*trapz(t(pos-1:pos),var);
    n(pos)=abs(n1+n2+n3+n4);

    % Denominador
    d1=(t(pos).^2).*(d1/(t(pos-1).^2)+1.5*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    d2=t(pos).*(d2./t(pos-1)-4*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    d3=d3+2.5*trapz(t(pos-1:pos),var);
    d(pos)=abs(d1+d2+d3);
    
    % Frecuencia
    f1=f1+trapz(t(pos-1:pos),n(pos-1:pos));
    f2=f2+trapz(t(pos-1:pos),d(pos-1:pos)); 
    f(pos)=sqrt(f1/f2);
    
   
    % Seleccion de la frecuencia
    if (pos > rango + margen)
        s1=sum(f(pos-rango:pos));
        media=s1/rango;

        if (media ~= Inf())
            s2=sum((f(pos-rango:pos)-media).^2);
            tolerancia=sqrt(s2/rango)/abs(media);
        end
    end
    resultados(pos)=tolerancia;

    pos=pos+1;
end

freq=f(pos-1);
frecuencia=freq/(2*pi)

% Una vez obtenida la frecuencia se determina la amplitud, fase y el offset.        
t1=1;
t2=pos-1;

B1=[trapz(y_valores(t1:t2).*sin(freq*t(t1:t2))); trapz(y_valores(t1:t2).*cos(freq*t(t1:t2))); trapz(y_valores(t1:t2))];
A1=[trapz(sin(freq*t(t1:t2)).^2),trapz(cos(freq*t(t1:t2)).*sin(freq*t(t1:t2))), trapz(sin(freq*t(t1:t2))); trapz(cos(freq*t(t1:t2)).*sin(freq*t(t1:t2))),trapz(cos(freq*t(t1:t2)).^2), trapz(cos(freq*t(t1:t2)));  trapz(sin(freq*t(t1:t2))), trapz(cos(freq*t(t1:t2))), t2-t1];
x=A1\B1;

fase=atan(x(2)/x(1))
Amplitud=sqrt(x(1)^2+x(2)^2)
Const=x(3)

y_aprox=Amplitud*sin(freq*t+fase)+Const;
toc;

% Se representa la frecuencia.
figure;
plot(t,f/(2*pi)); hold on;
axis([0 2 0 8]);
xlabel('t(s)'); ylabel('frecuencia (Hz)')

% Se compara la señal estimada con la señal original.
figure;
plot(t,y_valores(1:length(t))); hold on;
plot(t,y_aprox);
