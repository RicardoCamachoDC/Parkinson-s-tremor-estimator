% Lectura de los datos de la señal.

filename='Tremor.xlsx';
datos=xlsread(filename);

% La columna 1 representa el numero de muestra, la columna 2 el valor en el
% eje X, la columna 3 en el eje Y y la columna 4 en el eje Z.

%Se define el marco de tiempo en el que se trabaja.
t_sim=1;
h=0.01; % Tiempo de muestreo.

t=0:h:t_sim-h;

%Se representa los valores que se obtiene del sensor:
figure;
plot(t,datos(:,2)); hold on;
plot(t,datos(:,3)); hold on;
plot(t,datos(:,4)); hold on;
xlabel('t(s)'); ylabel('Aceleración (g)');
legend('Eje X','Eje Y','Eje Z');

% % Configuración del centinela.
rango=20;
condicion=0.0708;

% % Configuración del límite de tiempo y el margen de omisión.
margen=10;
limite=1;

t=0:h:limite-h;

% % Configuración del filtro paso bajo.
fc=10;

% % Estimación de las frecuencias.
% Se inicia un contador para ver cuanto tiempo tarda en procesar los datos.
tic;

%Condiciones iniciales y asignación de los datos a cada eje.
pos=2;
tolerancia_x=Inf();
tolerancia_y=Inf();
tolerancia_z=Inf();

y_valores_x=datos(:,2);
y_valores_x=y_valores_x.';
y_valores_x=lowpass(y_valores_x,fc,1/h);

n_x=zeros(1,length(t));
d_x=zeros(1,length(t));
f_x=zeros(1,length(t));
resultados_x=zeros(1,length(t));

y_valores_y=datos(:,3);
y_valores_y=y_valores_y.';
y_valores_y=lowpass(y_valores_y,fc,1/h);

n_y=zeros(1,length(t));
d_y=zeros(1,length(t));
f_y=zeros(1,length(t));
resultados_y=zeros(1,length(t));

y_valores_z=datos(:,4);
y_valores_z=y_valores_z.';
y_valores_z=lowpass(y_valores_z,fc,1/h);

n_z=zeros(1,length(t));
d_z=zeros(1,length(t));
f_z=zeros(1,length(t));
resultados_z=zeros(1,length(t));

% Ejecucion inicial del programa.
% Numerador Eje x
n1_x=(t(pos).^3).*y_valores_x(pos);
var=y_valores_x(1:pos);
n2_x=-3*(t(pos).^2).*trapz(t(1:pos),var);
var=var.*t(1:pos);
n3_x=+24*t(pos).*trapz(t(1:pos),var);
var=var.*t(1:pos);
n4_x=-30.*trapz(t(1:pos),var);
n_x(pos)=abs(n1_x+n2_x+n3_x+n4_x);
% Denominador Eje x
d1_x=1.5*(t(pos).^2).*trapz(t(1:pos),var);
var=var.*t(1:pos);
d2_x=-4*t(pos).*trapz(t(1:pos),var);
var=var.*t(1:pos);
d3_x=2.5*trapz(t(1:pos),var);
d_x(pos)=abs(d1_x+d2_x+d3_x);
% Frecuencia Eje x
f1_x=trapz(t(1:pos),n_x(1:pos));
f2_x=trapz(t(1:pos),d_x(1:pos));
f_x(pos)=sqrt(f1_x/f2_x);      
resultados_x(pos)=tolerancia_x;

% Numerador Eje Y
n1_y=(t(pos).^3).*y_valores_y(pos);
var=y_valores_y(1:pos);
n2_y=-3*(t(pos).^2).*trapz(t(1:pos),var);
var=var.*t(1:pos);
n3_y=+24*t(pos).*trapz(t(1:pos),var);
var=var.*t(1:pos);
n4_y=-30.*trapz(t(1:pos),var);
n_y(pos)=abs(n1_y+n2_y+n3_y+n4_y);
% Denominador Eje Y
d1_y=1.5*(t(pos).^2).*trapz(t(1:pos),var);
var=var.*t(1:pos);
d2_y=-4*t(pos).*trapz(t(1:pos),var);
var=var.*t(1:pos);
d3_y=2.5*trapz(t(1:pos),var);
d_y(pos)=abs(d1_y+d2_y+d3_y);
% Frecuencia Eje Y
f1_y=trapz(t(1:pos),n_y(1:pos));
f2_y=trapz(t(1:pos),d_y(1:pos));
f_y(pos)=sqrt(f1_y/f2_y);      
resultados_y(pos)=tolerancia_y;

% Numerador Eje Z
n1_z=(t(pos).^3).*y_valores_z(pos);
var=y_valores_z(1:pos);
n2_z=-3*(t(pos).^2).*trapz(t(1:pos),var);
var=var.*t(1:pos);
n3_z=+24*t(pos).*trapz(t(1:pos),var);
var=var.*t(1:pos);
n4_z=-30.*trapz(t(1:pos),var);
n_z(pos)=abs(n1_z+n2_z+n3_z+n4_z);
% Denominador Eje Z
d1_z=1.5*(t(pos).^2).*trapz(t(1:pos),var);
var=var.*t(1:pos);
d2_z=-4*t(pos).*trapz(t(1:pos),var);
var=var.*t(1:pos);
d3_z=2.5*trapz(t(1:pos),var);
d_z(pos)=abs(d1_z+d2_z+d3_z);
% Frecuencia Eje Z
f1_z=trapz(t(1:pos),n_z(1:pos));
f2_z=trapz(t(1:pos),d_z(1:pos));
f_z(pos)=sqrt(f1_z/f2_z);      
resultados_z(pos)=tolerancia_z;


pos=pos+1;
% Ejecución de forma recursiva
while( pos<=length(t) &&  (tolerancia_x>=condicion && tolerancia_y>=condicion && tolerancia_z>=condicion))
    % Numerador Eje X
    n1_x=(t(pos).^3).*y_valores_x(pos);
    var=y_valores_x(pos-1:pos);
    n2_x=(t(pos).^2).*(n2_x./(t(pos-1).^2)-3*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    n3_x=t(pos).*(n3_x./t(pos-1)+24*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    n4_x=n4_x-30.*trapz(t(pos-1:pos),var);
    n_x(pos)=abs(n1_x+n2_x+n3_x+n4_x);
    % Denominador Eje X
    d1_x=(t(pos).^2).*(d1_x/(t(pos-1).^2)+1.5*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    d2_x=t(pos).*(d2_x./t(pos-1)-4*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    d3_x=d3_x+2.5*trapz(t(pos-1:pos),var);
    d_x(pos)=abs(d1_x+d2_x+d3_x);
    % Frecuencia Eje X
    f1_x=f1_x+trapz(t(pos-1:pos),n_x(pos-1:pos));
    f2_x=f2_x+trapz(t(pos-1:pos),d_x(pos-1:pos)); 
    f_x(pos)=sqrt(f1_x/f2_x);

    % Numerador Eje Y
    n1_y=(t(pos).^3).*y_valores_y(pos);
    var=y_valores_y(pos-1:pos);
    n2_y=(t(pos).^2).*(n2_y./(t(pos-1).^2)-3*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    n3_y=t(pos).*(n3_y./t(pos-1)+24*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    n4_y=n4_y-30.*trapz(t(pos-1:pos),var);
    n_y(pos)=abs(n1_y+n2_y+n3_y+n4_y);
    % Denominador Eje Y
    d1_y=(t(pos).^2).*(d1_y/(t(pos-1).^2)+1.5*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    d2_y=t(pos).*(d2_y./t(pos-1)-4*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    d3_y=d3_y+2.5*trapz(t(pos-1:pos),var);
    d_y(pos)=abs(d1_y+d2_y+d3_y);
    % Frecuencia Eje Y
    f1_y=f1_y+trapz(t(pos-1:pos),n_y(pos-1:pos));
    f2_y=f2_y+trapz(t(pos-1:pos),d_y(pos-1:pos)); 
    f_y(pos)=sqrt(f1_y/f2_y);
    
     % Numerador Eje Z
    n1_z=(t(pos).^3).*y_valores_z(pos);
    var=y_valores_z(pos-1:pos);
    n2_z=(t(pos).^2).*(n2_z./(t(pos-1).^2)-3*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    n3_z=t(pos).*(n3_z./t(pos-1)+24*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    n4_z=n4_z-30.*trapz(t(pos-1:pos),var);
    n_z(pos)=abs(n1_z+n2_z+n3_z+n4_z);
    % Denominador Eje Z
    d1_z=(t(pos).^2).*(d1_z/(t(pos-1).^2)+1.5*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    d2_z=t(pos).*(d2_z./t(pos-1)-4*trapz(t(pos-1:pos),var));
    var=var.*t(pos-1:pos);
    d3_z=d3_z+2.5*trapz(t(pos-1:pos),var);
    d_z(pos)=abs(d1_z+d2_z+d3_z);
    % Frecuencia Eje Z
    f1_z=f1_z+trapz(t(pos-1:pos),n_z(pos-1:pos));
    f2_z=f2_z+trapz(t(pos-1:pos),d_z(pos-1:pos)); 
    f_z(pos)=sqrt(f1_z/f2_z);
    
    % Seleccion de la frecuencia mediante el centinela
    %Eje X
    if (pos > rango + margen)
        s1=sum(f_x(pos-rango:pos));
        media=s1/rango;
        if (media ~= Inf())
            s2=sum((f_x(pos-rango:pos)-media).^2);
            tolerancia_x=sqrt(s2/rango)/abs(media);
        end
    end
    resultados_x(pos)=tolerancia_x;
    
    %Eje Y
    if (pos > rango + margen)
        s1=sum(f_y(pos-rango:pos));
        media=s1/rango;
        if (media ~= Inf())
            s2=sum((f_y(pos-rango:pos)-media).^2);
            tolerancia_y=sqrt(s2/rango)/abs(media);
        end
    end
    resultados_y(pos)=tolerancia_y;
    
    %Eje Z
    if (pos > rango + margen)
        s1=sum(f_z(pos-rango:pos));
        media=s1/rango;

        if (media ~= Inf())
            s2=sum((f_z(pos-rango:pos)-media).^2);
            tolerancia_z=sqrt(s2/rango)/abs(media);
        end
    end
    resultados_z(pos)=tolerancia_z;
    
    pos=pos+1;
end

%Asignación de la frecuencia para cada eje.
freq_x=f_x(pos-1);
frecuencia_x=freq_x/(2*pi)
freq_y=f_y(pos-1);
frecuencia_y=freq_y/(2*pi)
freq_z=f_z(pos-1);
frecuencia_z=freq_z/(2*pi)

% Si ningún eje cumple la condición del centinela, se realiza una media de
%tres frecuencia. En caso contrario, se selecciona la frecuencia que antes
%cumpla la condición del centinela.
if (pos>length(t))
    freq=(freq_x+freq_y+freq_z)/3;
    frecuencia=freq/(2*pi)
else
    real=[tolerancia_x<condicion tolerancia_y<condicion tolerancia_z<condicion];
    freq=max(real.*[freq_x freq_y freq_z]);
    frecuencia=freq/(2*pi)
end
% Se estima los parámetros de la amplitud, fase y offset y se recrea
% cada señal de forma independiente.
t1=1;
t2=pos-1;

% Eje X
B1=[trapz(y_valores_x(t1:t2).*sin(freq*t(t1:t2))); trapz(y_valores_x(t1:t2).*cos(freq*t(t1:t2))); trapz(y_valores_x(t1:t2))];
A1=[trapz(sin(freq*t(t1:t2)).^2),trapz(cos(freq*t(t1:t2)).*sin(freq*t(t1:t2))), trapz(sin(freq*t(t1:t2))); trapz(cos(freq*t(t1:t2)).*sin(freq*t(t1:t2))),trapz(cos(freq*t(t1:t2)).^2), trapz(cos(freq*t(t1:t2)));  trapz(sin(freq*t(t1:t2))), trapz(cos(freq*t(t1:t2))), t2-t1];
x_x=A1\B1;

fase_x=atan(x_x(2)/x_x(1))
Amplitud_x=sqrt(x_x(1)^2+x_x(2)^2)
Const_x=x_x(3)

%Eje Y
B1=[trapz(y_valores_y(t1:t2).*sin(freq*t(t1:t2))); trapz(y_valores_y(t1:t2).*cos(freq*t(t1:t2))); trapz(y_valores_y(t1:t2))];
A1=[trapz(sin(freq*t(t1:t2)).^2),trapz(cos(freq*t(t1:t2)).*sin(freq*t(t1:t2))), trapz(sin(freq*t(t1:t2))); trapz(cos(freq*t(t1:t2)).*sin(freq*t(t1:t2))),trapz(cos(freq*t(t1:t2)).^2), trapz(cos(freq*t(t1:t2)));  trapz(sin(freq*t(t1:t2))), trapz(cos(freq*t(t1:t2))), t2-t1];
x_y=A1\B1;

fase_y=atan(x_y(2)/x_y(1))
Amplitud_y=sqrt(x_y(1)^2+x_y(2)^2)
Const_y=x_y(3)


%EJE Z
B1=[trapz(y_valores_z(t1:t2).*sin(freq*t(t1:t2))); trapz(y_valores_z(t1:t2).*cos(freq*t(t1:t2))); trapz(y_valores_z(t1:t2))];
A1=[trapz(sin(freq*t(t1:t2)).^2),trapz(cos(freq*t(t1:t2)).*sin(freq*t(t1:t2))), trapz(sin(freq*t(t1:t2))); trapz(cos(freq*t(t1:t2)).*sin(freq*t(t1:t2))),trapz(cos(freq*t(t1:t2)).^2), trapz(cos(freq*t(t1:t2)));  trapz(sin(freq*t(t1:t2))), trapz(cos(freq*t(t1:t2))), t2-t1];
x_z=A1\B1;

fase_z=atan(x_z(2)/x_z(1))
Amplitud_z=sqrt(x_z(1)^2+x_z(2)^2)
Const_z=x_z(3)

t=0:h:t_sim-h;
y_aprox_x=Amplitud_x*sin(freq*t+fase_x)+Const_x;
y_aprox_y=Amplitud_y*sin(freq*t+fase_y)+Const_y;
y_aprox_z=Amplitud_z*sin(freq*t+fase_z)+Const_z;

%Se finaliza el contador para observar cuanto tiempo computacional ha
%tardado en estimar las frecuencias.
toc;


% Se representa los diferentes datos que se ha obtenido para la frecuencia:
figure; 
subplot(3,1,1); plot(t,f_x./(2*pi));
xlabel('t(s)'); ylabel('frecuencia_x (Hz)');
axis([0 t_sim 0 10]);
subplot(3,1,2); plot(t,f_y./(2*pi));
xlabel('t(s)'); ylabel('frecuencia_y (Hz)');
axis([0 t_sim 0 10]);
subplot(3,1,3); plot(t,f_z./(2*pi));
xlabel('t(s)'); ylabel('frecuencia_z (Hz)');
axis([0 t_sim 0 10]);

% Se representa los datos del estimador algebraico en comparación con los
% datos de entrada.
figure;
subplot(3,1,1); plot(t,datos(:,2)); hold on; plot(t,y_aprox_x);
xlabel('t(s)'); ylabel('Aceleración_x (g)'); 
subplot(3,1,2); plot(t,datos(:,3)); hold on; plot(t,y_aprox_y);
xlabel('t(s)'); ylabel('Aceleración_y (g)'); 
subplot(3,1,3); plot(t,datos(:,4)); hold on; plot(t,y_aprox_z);
xlabel('t(s)'); ylabel('Aceleración_z (g)'); 
legend('Real','Estimada');

% Se realiza la inversa para de la señal estimada sobre los datos reales.
figure;
subplot(3,1,1); plot(t,y_aprox_x-datos(:,2)'); 
xlabel('t(s)'); ylabel('Ruido_x (g)'); 
subplot(3,1,2); plot(t,y_aprox_y-datos(:,3)');
xlabel('t(s)'); ylabel('Ruido_y (g)'); 
subplot(3,1,3); plot(t,y_aprox_z-datos(:,4)');
xlabel('t(s)'); ylabel('Ruido_z (g)'); 