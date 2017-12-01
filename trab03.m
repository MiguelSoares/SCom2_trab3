%% ------------------------------ Trabalho 3 ------------------------------
% Diogo Correia (dv.correia@ua.pt)

close all; clear; clc;
savePlots = 0; % '1' to save plots

% 2. QAM Signal ----------------------------------------------------------
% 2.1 --------------------------------------------------------------------

% INPUT PROGRAM VARIABLES 

M = 4;    % Number of levels per symbol
m = 3;    % Number of symbols to simulate 
n = 200;  % Number of samples per symbol 
d = 1;    % Difference between consecutive levels in the baseband modulating signals
Rc = 4;   % Ration between the carrier frequency and symbol rate
T = 1;    % T is the symbol duration

% QAM generation Testing (optional)
nbits = log2(M);
nsym = m*n*nbits/n;
dsym = nbits/n;

t = 0 : dsym : nsym - dsym;
figure('Name','2.1 - 16-QAM and Modulated');
set(gcf, 'Position', [100, 100, 800, 400]);

A = mary(M,m,n) .* d/2;
subplot(2,4,1);
stem(t,A,'Marker','none', 'Color', 'b'); grid on;
axis([0 nsym (-M+1)*d/2 (M-1)*d/2]);
title(['A(t) - ',num2str(m),' Symbols ',num2str(M^2),'-QAM'],'Color','b');
xlabel('t(bits)'); ylabel('A(t)');

B = mary(M,m,n) .* d/2;
subplot(2,4,5);
stem(t,B,'Marker','none','Color','r'); grid on;
axis([0 nsym (-M+1)*d/2 (M-1)*d/2]);
title(['B(t) - ',num2str(m),' Symbols ',num2str(M^2),'-QAM'],'Color','r');
xlabel('t(bits)'); ylabel('B(t)');

% Modulation
r = sqrt(2/T)*cos(2*pi*Rc*t/nbits);
i = sqrt(2/T)*sin(2*pi*Rc*t/nbits);

% Modulated Signals
subplot(2,4,2);
plot(t,A.*r,'b');
title(['C(t) - A(t) Modulated'],'Color','b');
xlabel('t(bits)'); ylabel('A(t) Modulated');
axis([0 nsym (-M-1)*d/2 (M+1)*d/2]);

subplot(2,4,6);
plot(t,B.*i,'r');
title(['D(t) - B(t) Modulated'],'Color','r');
xlabel('t(bits)'); ylabel('B(t) Modulated');
axis([0 nsym (-M-1)*d/2 (M+1)*d/2]);

subplot(2,4,[3,4,7,8]);
hold on; grid minor;
xlabel('A'); ylabel('B');
title('Constelation');

% Centering the axis on origin
try
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
catch
    warning('Your version of matlab migth not have the centering axis functions, so they will not run. They will be substituted with blue lines!');
    plot([-1.5 1.5],[0 0],'b');
    plot([0 0],[-1.5 1.5],'b');
end

% % Generating random symboles to populate the constelation
Q = [mary(M,1000,1) .*d/2 ; mary(M,1000,1) .*d/2];
plot(Q(1,:),Q(2,:),'*b');

for tt = 1 : n : m*n
    plot(A(tt),B(tt),'*r');
    text(A(tt),B(tt),['\leftarrow simb',num2str((tt-1)/n + 1)])
end

% Save plot
if savePlots == 1
    saveas(gcf,'plots/1-teste.png');
end

%% 2.2 --------------------------------------------------------------------
close all; clearvars -except savePlots; clc;

m = 1000;       % Number of symbols to simulate 
n = 200;        % Number of samples per symbol 
KI = 8;         % Number of levels of the in-phase baseband modulating signal 
KQ = 8;         % Number of levels of the quadrature baseband modulating signal 
d = 1;          % Difference between consecutive levels in the baseband modulating signals
Rc = 4;         % Ration between the carrier frequency and symbol rate
T = 1;          % T is the symbol duration
M = KI * KQ;    % Number of levels per symbol

ns = 3;         % Number of symbols to visualize
R = 1e6;        % Transmition rate in bits/s;
Rs = R/log(M);  % Transmition rate in symbols/s
fa = n*R;       % Sampling frequency
fp = Rs*Rc;     % Carrier frequency

% a) ----------------------------------------------------------------------
% Generating random symboles to populate the constelation
figure('Name','2.2 - 64-QAM and constellation');
hold on; grid on;

A = mary(KQ,m,n) * d/2;
B = mary(KI,m,n) * d/2;

scatter(A,B,15,'b','filled');

for i = 1 : 3
    N = i * n;
    scatter(A(N),B(N),30,'r','filled');
    text(A(N),B(N),['\leftarrow simb',num2str(i)])
    xlabel('A'); ylabel('B');
end

xlabel('A'); ylabel('B');
title('Constelation 64-QAM');

% Centering the axis on origin
try
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
catch
    warning('Your version of matlab migth not have the centering axis functions, so they will not run. They will be substituted with blue lines!');
    plot([-4 4],[0 0],'b');
    plot([0 0],[-4 4],'b');
end

% Save plot
if savePlots == 1
    saveas(gcf,'plots/2-constelation.png');
end

% b) ----------------------------------------------------------------------

nbits = log2(M);
nsym = m*n*nbits/n;
dsym = nbits/n;

t = 0 : dsym : nsym - dsym;

% Modulation
r = sqrt(2/T)*cos(2*pi*Rc*t/nbits);
i = sqrt(2/T)*sin(2*pi*Rc*t/nbits);

C = A.*r;
D = B.*i;
E = C + D;

% Averange symbol energy
AET = (M - 1) * (d^2) / 6;
AEPconstelation = 1/m * sum(A(1:1000).^2 + B(1:1000).^2);
AEPcarrier = 1/m * sum(C(1:1000).^2 + D(1:1000).^2);

disp(['Averange symbol energy (Theoretical) = ',num2str(AET)]);
disp(['Averange symbol energy (simulated constelation) = ', num2str(AEPconstelation)]);
disp(['Averange symbol energy (simulated carrier) = ', num2str(AEPcarrier)]);

% c) ----------------------------------------------------------------------

figure('Name','2.2 - 64-QAM 3 symbol baseband signal');

subplot(1,2,1);
stem(t(1:n*3),A(1:n*3),'Marker','none', 'Color', 'b'); grid on;
title(['A(t) - ',num2str(m),' Symbols ',num2str(M),'-QAM'],'Color','b');
xlabel('t(bits)'); ylabel('A(t)');
axis([0 3*nbits -KI*d/2 KI*d/2]);

subplot(1,2,2);
stem(t(1:n*3),B(1:n*3),'Marker','none', 'Color', 'r'); grid on;
title(['B(t) - 3 Symbols ',num2str(M),'-QAM'],'Color','r');
xlabel('t(bits)'); ylabel('B(t)');
axis([0 3*nbits -KQ*d/2 KQ*d/2]);

% Save plot
if savePlots == 1
    saveas(gcf,'plots/2-binary.png');
end

% d) ----------------------------------------------------------------------

points = [A((2:3).*n) ; B((2:3).*n)];

fprintf('\nVerifica-se comparando com a constelação:');
fprintf('\nConstelation = (A , B)\n');
disp('        sym1        sym2');
disp(['A :     ',num2str(points(1,:))]);
disp(['B :     ',num2str(points(2,:))]);

% e) ----------------------------------------------------------------------

% Modulation
r = sqrt(2/T)*cos(2*pi*Rc*t/nbits);
i = sqrt(2/T)*sin(2*pi*Rc*t/nbits);

Cs = A(1:3*n).*r(1:3*n);
Ds = B(1:3*n).*i(1:3*n);

figure('Name','Modulation: C(t), D(t) and E(t)');

subplot(3,1,1);
plot(t(1:3*n),Cs,'b');
title(['C(t)'],'Color','b');
xlabel('t(bits)'); ylabel('C(t)');
grid minor;

subplot(3,1,2);
plot(t(1:3*n),Ds,'r');
title('D(t)','Color','r');
xlabel('t(bits)'); ylabel('D(t)');
grid minor;

Em = Cs+Ds;
subplot(3,1,3);
plot(t(1:3*n),Em,'m');
title('E(t)','Color','m');
xlabel('t(bits)'); ylabel('E(t)');
grid minor;

% Save plot
if savePlots == 1
    saveas(gcf,'plots/2-modelação.png');
end

% f) ----------------------------------------------------------------------
% f.1) --------------------------------------------------------------------

figure('Name','DEP phase component baseband');

windowL = 16384;
[Pxx,f] = pwelch(A,hamming(windowL),windowL/2,windowL,n,'onesided');

plot(f/nbits, 10*log10(Pxx/1e-3),'b');
xlabel('f / R_{bit}'); ylabel('PSD (dB)');
title('PSD - phase component baseband');
axis([0 1.2 0 60]); grid on;

% f.2) --------------------------------------------------------------------
% <insert explanation here>

% f.3) --------------------------------------------------------------------
hold on;
plot(f/nbits, 10*log10((AET.*((sin(pi.*f.*T)./(pi.*f.*T)).^2))./1e-3),'r');
legend('Simulated PSD','Theorical PSD curve');

% Save plot
if savePlots == 1
    saveas(gcf,'plots/2-psd-baseband.png');
end

% f.4) --------------------------------------------------------------------
% <insert something here>

% g) ----------------------------------------------------------------------
figure('Name','DEP E component baseband');
grid on; hold on;

windowL=16384;
[Pxx,f] = pwelch(E,hamming(windowL),windowL/2,windowL,n,'onesided');

plot(f/nbits, 10*log10(Pxx/1e-3),'b');
xlabel('f / R_{bit}');
ylabel('unilateral PSD (dB)');
title('Signal E PSD (modulated carrier at the QAM transmitter output)');
axis([0 1.2 0 60]); 

% theorical PSD
plot(f/nbits, 10*log10(((5.25.*((sin(pi.*(f-Rc).*T)./(pi.*(f-Rc).*T)).^2))./1e-3)+((5.25.*((sin(pi.*(f+Rc).*T)./(pi.*(f+Rc).*T)).^2))./1e-3)),'r');
legend('Simulated PSD','Theorical PSD curve');

% Save plot
if savePlots == 1
    saveas(gcf,'plots/2-psd-modulatedcarrier.png');
end

% h) ----------------------------------------------------------------------
KI = 64;         % Number of levels of the in-phase baseband modulating signal 
KQ = 64;         % Number of levels of the quadrature baseband modulating signal 
M = KI * KQ;     % Number of levels per symbol

d = sqrt((6 * 10.5)/(M - 1));

A = mary(KQ,m,n) * d/2;
B = mary(KI,m,n) * d/2;

nbits = log2(M);
nsym = m*n*nbits/n;
dsym = nbits/n;

t = 0 : dsym : nsym - dsym;

% Modulation
r = sqrt(2/T)*cos(2*pi*Rc*t/nbits);
i = sqrt(2/T)*sin(2*pi*Rc*t/nbits);

C = A.*r;
D = B.*i;
E = C + D;

figure('Name','DEP E 4096-QAM');
grid on; hold on;

windowL=16384;
[Pxx,f] = pwelch(E,hamming(windowL),windowL/2,windowL,n,'onesided');

plot(f/nbits, 10*log10(Pxx/1e-3),'b');
xlabel('f / R_{bit}');
ylabel('unilateral PSD (dB)');
title('Signal E PSD 4096-QAM');
axis([0 1.2 0 60]);

% Save plot
if savePlots == 1
    saveas(gcf,'plots/2-psd-E-4096qam.png');
end

%% 3. System performance --------------------------------------------------
%% 3.1 --------------------------------------------------------------------
close all; clear; clc;

n = 1000;       % Number of symbols to simulate   
d = 2;          % Difference between consecutive levels in the baseband modulating signals
M = 4;

map = [-1 -1 0 0 ;
       -1 1 0 1 ; 
       1 1 1 1 ; 
       1 -1 1 0];

% 3.2 ---------------------------------------------------------------------
% a) ----------------------------------------------------------------------   
dsource = randi([1,4], 1, n);

qam = zeros(n,2);
binary = zeros(n,2);

% Coversão para binário e QAM
for i = 1:n
    qam(i,:) = map(dsource(i),1:2) * d/2;
    binary(i,:) = map(dsource(i),3:4);
end

Es = (d^2)/2;
dv = 5.5;

Pes = 10e-4;
n0 = Es/(2*(erfcinv((-sqrt(-Pes+1)+1)*2))^2);

disp(['Energy per symbol = ',num2str(Es)]);
disp(['Noise Power = ',num2str(n0)]);
disp(['Pes = ',num2str(Pes)]);

noise = wgn(n,2,10*log10(n0));
qamt = qam.*(dv/2);

r = qamt + noise;

% b) ----------------------------------------------------------------------
figure; hold on; grid on;
scatter(r(:,1),r(:,2),3,'r','filled');
scatter(qamt(:,1),qamt(:,2),25,'b','filled');
%axis([-1.5 1.5 -1.5 1.5]);
xlabel('A'); ylabel('B');
title('4-QAM Constelation');

% Centering the axis on origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% Naming the quadrant binarys
ylim = get(gca,'ylim');
xlim = get(gca,'xlim');
text(xlim(1),ylim(1),'[0,0]','Color','b','FontSize',14);
text(xlim(1),ylim(2),'[0,1]','Color','b','FontSize',14);
text(xlim(2),ylim(2),'[1,1]','Color','b','FontSize',14)
text(xlim(2),ylim(1),'[1,0]','Color','b','FontSize',14);

% c) ----------------------------------------------------------------------
Es = (M-1)*(d^2)/6;
n0 = mean(noise(1,:).^2+noise(2,:).^2);
Pes = 1 - (1-1/2*erfcinv(sqrt(Es/(2*n0))))^2;

disp(['Energy per symbol = ',num2str(Es)]);
disp(['Averange noise Power = ',num2str(n0)]);
disp(['Pes = ',num2str(Pes)]);

% d) ----------------------------------------------------------------------
j = zeros(1,3);

% Optimum QAM receiver

% e) ----------------------------------------------------------------------

% 3.3 ---------------------------------------------------------------------
