% ELEC 4700 MNA Building

R1 = 1;
Cap = 0.25;
R2 =  2;
L = 0.2;
R3 = 10;
alpha = 100;
R4 = 0.1;
RO = 1000;


G = [
    1.0000   -1.0000         0         0         0         0         0    1.0000    ;
   -1.0000    1.5000         0         0         0    1.0000         0         0    ;
         0         0    0.1000         0         0   -1.0000         0         0    ;
         0         0         0   10.0000  -10.0000         0    1.0000         0    ;
         0         0         0  -10.0000   10.0010         0         0         0    ;
         0    1.0000   -1.0000         0         0         0         0         0    ;
         0         0  -10.0000    1.0000         0         0         0         0    ;
    1.0000         0         0         0         0         0         0         0    ;
    ];
    
C = [
       Cap      -Cap         0         0         0         0         0         0    ;
      -Cap       Cap         0         0         0         0         0         0    ;
         0         0         0         0         0         0         0         0    ;
         0         0         0         0         0         0         0         0    ;
         0         0         0         0         0         0         0         0    ;
         0         0         0         0         0        -L         0         0    ;
         0         0         0         0         0         0         0         0    ;
         0         0         0         0         0         0         0         0    ;
    ];

F = [
    0   ;
    0   ;
    0   ;
    0   ;
    0   ;
    0   ;
    0   ;
    1   ;
    ];


%%%%% DC Sweep %%%%%
fig_dc = figure;
hold on;
V1 = [];
V5 = [];
for Vin=-10:0.1:10
    F(8) = Vin; % Set DC voltage
    % Solve for DC
    e = G\F;
    V1 = [V1 e(1)];
    V5 = [V5 e(5)];
end
title('DC Sweep');
xlabel('Vin (V)');
ylabel('Node Voltage (V)');
plot(-10:0.1:10, V1);
plot(-10:0.1:10, V5);
legend('V1', 'VO');

%%%%% AC Sweep %%%%%
fig_ac = figure;
hold on;
V5 = [];
F(8) = 1; % Set DC voltage
for w=10:100000
    e = (G+2*pi*w*1j*C)\F;
    V5 = [V5 20*log10(abs(e(5)/F(8)))];
end
title('AC Sweep');
xlabel('f (Hz)');
ylabel('Gain (dB)');
semilogx((10:100000)/180*pi, V5);
    
%%%%% Cap Sweep %%%%%
fig_cap = figure;
hold on;
V5 = [];
F(8) = 1; % Set DC voltage
std = 0.05;
c_dist = std.*randn(50000,1) + Cap;
w = pi;
for index=1:50000
    C(1,1) = c_dist(index);
    C(1,2) = -c_dist(index);
    C(2,1) = -c_dist(index);
    C(2,2) = c_dist(index);
    e = (G+2*pi*w*1j*C)\F;
    V5 = [V5 20*log10(abs(e(5)/F(8)))];
end
title('Capacitor Sweep');
xlabel('Gain (dB)');
histogram(V5);