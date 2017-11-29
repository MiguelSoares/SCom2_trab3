%% 2.1

close all;
clc;
clear all;
m=600;  %numero simbolos a simular
n=100;   %numero de amostras a simular
KI=8;   %numero niveis do sinal em banda base
KQ=8;   %numero niveis do sinal em banda quadrante
M=KQ*KI;    %numero de simbolos da constelaçao
d=1;    %diferenca entre niveis
T=1;

Rc=4;   %relação entre freq da porta e a taxa de simbolos

ns=4;   %numeros de simbolos a visualizar
R=1e6;  %taxa de tramissao em bits/s
Rs=R/log2(M);   %taxa de transmissao em simbolos/s
fa=n*R; %freq de amostragem
fp=Rs*Rc;   %frequencia de portadora

%yd -> sinal banda base em fase
%yq -> sinal banda base em quadratura
%ydm -> portadora 1 modulada
%yqm -> portadora 2 modulada
%s -> sinal modulado

xd=mary(KI,m,n);
xq=mary(KQ,m,n);

nbits=log2(KI);
nsym=(0:(m*n-1))*nbits/n;

yd=xd*d/2;
yq=xq*d/2;

figure
stem(nsym,yd)
axis([0 9 -4 4])

phi1=sqrt(2/T)*cos(2*pi*Rc*nsym/nbits);
phi2=sqrt(2/T)*sin(2*pi*Rc*nsym/nbits);

ydm=yd.*phi1;

%plot(nsym,ydm)
%axis([0 9 -8 8])
%figure
yqm=yq.*phi2;

s=ydm+yqm;



%% 2.2 a)
close all

axisx=(-(KI-1):2:(KI-1))*d/2;
axisy=(-(KQ-1):2:(KQ-1))*d/2;

for i=axisx
    for k=axisy
        plot(i,k,'b o')
        hold on
    end
end
grid on;
axis(d*[min(-KI/2, -KQ/2) max(KI/2, KQ/2) min(-KI/2, -KQ/2) max(KI/2, KQ/2)])
hold on
p=0;

YD=zeros(1,m);
YQ=zeros(1,m);
YDM=zeros(1,m);
YQM=zeros(1,m);

for i=1:n:length(yd)
    p=p+1;
    plot(yd(i),yq(i),'*k')
    text(yd(i),yq(i),['  simb',int2str(p)]);
    YD(p)=yd(i);
    YQ(p)=yq(i);
    YDM(p)=ydm(i);
    YQM(p)=yqm(i);
    
end
title('Constelação 64-QAM')
xlabel('Componente em fase')
ylabel('Componente em quadratura')
figure

%% 2.2 b)
close all

As_energy_theoretical=(M-1)*(d^2)/6;
As_energy_pratical_constelation=(1/m)*sum(YD.^2+YQ.^2);
As_energy_pratical_carrier=(1/m)*sum(YDM.^2+YQM.^2);


%% 2.2 c)

%feito anteriormente

%% 2.2 d)

%comprova-se

%% 2.2 e)

%feito anteriormente

%% 2.2 f)

% [Pss,f]=pwelch(yd,1024,0,1024,fa,'onesided');
% figure
% plot(f,10*log10(Pss*1e3))
% hold on
windowL=16384;
windowT=hamming(windowL);
[Pxx,f]=pwelch(yd,windowT,windowL/2,windowL,n,'onesided');
plot(f/6,10*log10(Pxx/1e-3))
xlabel('f / R_{bit}');
ylabel('DEP bilateral (dB)')
title('DEP componente em fase banda base')
axis([0 1.2 0 60]), grid on
hold on


%% 2.2 f3)

plot(f/6, 10*log10((As_energy_theoretical.*((sin(pi.*f.*T)./(pi.*f.*T)).^2))./1e-3))
xlabel('f / R_{bit}');
ylabel('DEP bilateral (dB)')
title('DEP componente em fase banda base - Teorico')
axis([0 1.2 0 60])
grid on
hold on

%% 2.2 g)

windowL=16384;
windowT=hamming(windowL);
[Pxx,f]=pwelch(s,windowT,windowL/2,windowL,n,'onesided');
plot(f/6,10*log10(Pxx/1e-3))
xlabel('f / R_{bit}');
ylabel('DEP bilateral (dB)')
title('DEP componente em fase banda base')
axis([0 1.2 0 60]), grid on
hold on

%% 2.2 g) teorico


plot(f/6, 10*log10(((5.25.*((sin(pi.*(f-Rc).*T)./(pi.*(f-Rc).*T)).^2))./1e-3)+((5.25.*((sin(pi.*(f+Rc).*T)./(pi.*(f+Rc).*T)).^2))./1e-3)))
xlabel('f / R_{bit}');
ylabel('DEP bilateral (dB)')
title('DEP componente em fase banda base - Teorico')
axis([0 1.2 0 60])
grid on
hold on

%% h)

%Repetir para 1096 QAM comparar os valores dos globos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3)
 %% a)
x=inline('Es/(((erfcinv((sqrt(-Pes+1)-1)*-2))^2)*2)','Es','Pes');
d=5.5;
Es=(d^2)/2;
Pes=10e-4;
n0=x(Es,Pes);
n0=10*log10(n0/1e-3);

%% b)
map=[-1 -1 0 0;-1 1 0 1; 1 1 1 1; 1 -1 1 0];
map=map(:,1:2)*d/2;
M=4;
m=1000;

s=randint(1,m*log2(M));
s=reshape(s,log2(M),m)';

temp=zeros(size(s,1),1);
for i=1:log2(M)
    temp=s(:,i)*2^(log2(M)-i)+temp;
end

s_phase=zeros(size(temp,1),1);
for i=1:size(temp,1)
    s_phase(i)=map(temp(i)+1,1);
end

s_quadrature=zeros(size(temp,1),1);
for i=1:size(temp,1)
    s_quadrature(i)=map(temp(i)+1,2);
end

Pnoise=(10^(n0/10))/1000;   %Pruido (W/Hz)

noise_f=wgn(m,1,10*log10(Pnoise/2)); %ruido entre -d e d
noise_q=wgn(m,1,10*log10(Pnoise/2)); %ruido entre -d e d

s_phase_noise=s_phase+noise_f;
s_quadrature_noise=s_quadrature+noise_q;

figure();
plot(s_phase_noise,s_quadrature_noise,'r+');
hold on;
grid on;

plot(s_phase,s_quadrature,'go','LineWidth',2);
title('Constelation without noise(green) and with noise(red)');
xlabel('Phase component');
ylabel('Quadrature component');
hold off;

%% c)

Es_WithoutNoise=1/m*(sum(s_phase.^2+s_quadrature.^2));
Eb_WithoutNoise=Es_WithoutNoise/log2(M);
Pes_Pratical=(1-(1-0.5*erfc(sqrt(Es_WithoutNoise/(2*Pnoise))))^2);

clc;

fprintf('\n*************** 4-QAM System ****************\n');
fprintf('Distance between the constellation points: %g\n',d);
fprintf('Average energy per symbol: %g\n',Es_WithoutNoise);
fprintf('Average energy per bit, without noise: %g\n',Eb_WithoutNoise);
fprintf('Theoretical probability of symbol error: %g\n',Pes);
fprintf('Pratical probability of symbol error: %g\n',Pes_Pratical);

%% d)

map=[-1 -1 0 0;-1 1 0 1; 1 1 1 1; 1 -1 1 0];
map(:,1:2)=map(:,1:2)*d/2;
M=4;
m=3;

s=randint(1,m*log2(M));
s=reshape(s,log2(M),m)';

temp=zeros(size(s,1),1);
for i=1:log2(M)
    temp=s(:,i)*2^(log2(M)-i)+temp;
end

s_phase=zeros(size(temp,1),1);
for i=1:size(temp,1)
    s_phase(i)=map(temp(i)+1,1);
end

s_quadrature=zeros(size(temp,1),1);
for i=1:size(temp,1)
    s_quadrature(i)=map(temp(i)+1,2);
end

Pnoise=n0;   %Pruido (W/Hz)

noise_f=wgn(m,1,10*log10(Pnoise/2)); %ruido entre -d e d
noise_q=wgn(m,1,10*log10(Pnoise/2)); %ruido entre -d e d

s_phase_noise=s_phase+noise_f;
s_quadrature_noise=s_quadrature+noise_q;

% MAPPING OF THE CONSTELATION

figure();
plot(map(:,1),map(:,2),'go','LineWidth',2,'DisplayName','Constellation points');
% for i=1:length(map(:,1))
%    text(map(i,1),map(i,2),['  bits:',int2str(map(i,3))',int2str(map(i,4))']); 
% end
hold on

%PLOT OF THE SYMBOLS (OUTPUT OF QAM TRANSMITTER)
plot(s_phase,s_quadrature,'b+','MarkerSize',15,'LineWidth',2,'DisplayName','Output of QAM Transmitter');
for i=1:length(s_phase)
   text(s_phase(i),s_quadrature(i),['  BB:',int2str(i)]); 
end

%PLOT OF THE SYMBOLS (INPUT OF THE RECEIVER)
plot(s_phase_noise,s_quadrature_noise,'r+','MarkerSize',7,'DisplayName','Input of the Receveir');
for i=1:length(s_phase_noise)
   text(s_phase_noise(i),s_quadrature_noise(i),['          IR:',int2str(i)]); 
end
hold on;
grid on;

title('Constelation without noise(green) and with noise(red)');
xlabel('Phase component');
ylabel('Quadrature component');
axis([-6 6 -6 6])
hold on;

% RECEVEIR "MAXIMUM LIKELIHOOD" OUTPUT 
Ei=[map(:,1).^2+map(:,2).^2]*ones(1,m);
comp=map(:,1)*s_phase_noise'+map(:,2)*s_quadrature_noise'-Ei/2;
[x,idx]=max(comp);
phase_out=map(idx,1)';
quadrature_out=map(idx,2)';

plot(phase_out(1:3),quadrature_out(1:3),'mX','LineWidth',2,'DisplayName','Output of the Receveir','MarkerSize',10);
for i=1:length(phase_out(1:3))
   text(phase_out(i),quadrature_out(i),['OR:',int2str(i),'     '],'HorizontalAlignment','right'); 
end
legend('show')
hold on;

%% e)

%Esta certo

%% 3.3 a)
r=1; %escolhe o primeiro dos dois tipos de erro
d_calculator=inline('sqrt(Es*2)','Es');
Es=1;
d=d_calculator(Es);
Pes=[10e-3 0.5];
n0_calculator=inline('Es/(((erfcinv((sqrt(-Pes+1)-1)*-2))^2)*2)','Es','Pes');
n0=n0_calculator(Es,Pes(r));
m=10000;
Pesa=inline('erfc(sqrt(Es/(2*n0)))','Es','n0');
P_esa=Pesa(Es,n0);

err=zeros(1,5);
Pe_pratical=zeros(1,5);
for k=1:5
    map=[-1 -1 0 0;-1 1 0 1; 1 1 1 1; 1 -1 1 0];
    map(:,1:2)=map(:,1:2)*d/2;
    M=4;

    s=randint(1,m*log2(M));
    s=reshape(s,log2(M),m)';

    temp=zeros(size(s,1),1);
    for i=1:log2(M)
        temp=s(:,i)*2^(log2(M)-i)+temp;
    end

    s_phase=zeros(size(temp,1),1);
    for i=1:size(temp,1)
        s_phase(i)=map(temp(i)+1,1);
    end

    s_quadrature=zeros(size(temp,1),1);
    for i=1:size(temp,1)
        s_quadrature(i)=map(temp(i)+1,2);
    end

    Pnoise=n0;   %Pruido (W/Hz)

    noise_f=wgn(m,1,10*log10(Pnoise/2)); %ruido entre -d e d
    noise_q=wgn(m,1,10*log10(Pnoise/2)); %ruido entre -d e d

    s_phase_noise=s_phase+noise_f;
    s_quadrature_noise=s_quadrature+noise_q;


    % RECEVEIR "MAXIMUM LIKELIHOOD" OUTPUT 
    Ei=[map(:,1).^2+map(:,2).^2]*ones(1,m);
    comp=map(:,1)*s_phase_noise'+map(:,2)*s_quadrature_noise'-Ei/2;
    [x,idx]=max(comp);

    phase_out=map(idx,1)';
    quadrature_out=map(idx,2)';

    err(k) = sum((phase_out(1:m) ~= s_phase(1:m)')|(s_quadrature(1:m)' ~= quadrature_out(1:m)));
    Pe_pratical(k) = err(k)/length(s_phase);
    
end
fprintf('-------------------4-QAM SIMULATION INI ------------------ \n')
fprintf('Noise (dBm/Hz): %.5f \n',10*log10(n0/1e-3))
fprintf('Number of symbols used: %.0f \n',m)
fprintf('SER (average of 5): %.5f \n',sum(Pe_pratical)/5)
fprintf('Pes: %.5f \n',Pes(r))
fprintf('|SER-Pes|/Pes (percentage): %.5f \n',((abs((sum(Pe_pratical)/5)-Pes(r)))/Pes(r))*100)
fprintf('P_esa: %.5f \n',P_esa)
fprintf('Ratio Pes/Pesa: %.5f \n',Pes(r)/P_esa)

