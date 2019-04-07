%% ELEC 4700 Assignment 4 Question 3
% Andrew Branicki 100973961

R1 = 1;
Cap = 0.25;
R2 =  2;
L = 0.2;
R3 = 10;
alpha = 100;
R4 = 0.1;
RO = 1000;
Cn_1 = 0.00001;
Cn_2 = 0.0001;
Cn_3 = 0.01;

G = [
    1.0000   -1.0000         0         0         0         0         0    1.0000  ;
   -1.0000    1.5000         0         0         0    1.0000         0         0  ;
         0         0    0.1000         0         0   -1.0000         0         0  ;
         0         0         0   10.0000  -10.0000         0    1.0000         0  ;
         0         0         0  -10.0000   10.0010         0         0         0  ;
         0    1.0000   -1.0000         0         0         0         0         0  ;
         0         0  -10.0000    1.0000         0         0         0         0  ;
    1.0000         0         0         0         0         0         0         0  ;
    ];
    
C_1 = [
       Cap      -Cap         0         0         0         0         0         0  ;
      -Cap       Cap         0         0         0         0         0         0  ;
         0         0         Cn_1      0         0         0         0         0  ;
         0         0         0         0         0         0         0         0  ;
         0         0         0         0         0         0         0         0  ;
         0         0         0         0         0        -L         0         0  ;
         0         0         0         0         0         0         0         0  ;
         0         0         0         0         0         0         0         0  ;
    ];

C_2 = [
       Cap      -Cap         0         0         0         0         0         0  ;
      -Cap       Cap         0         0         0         0         0         0  ;
         0         0         Cn_2      0         0         0         0         0  ;
         0         0         0         0         0         0         0         0  ;
         0         0         0         0         0         0         0         0  ;
         0         0         0         0         0        -L         0         0  ;
         0         0         0         0         0         0         0         0  ;
         0         0         0         0         0         0         0         0  ;
    ];

C_3 = [
       Cap      -Cap         0         0         0         0         0         0  ;
      -Cap       Cap         0         0         0         0         0         0  ;
         0         0         Cn_3      0         0         0         0         0  ;
         0         0         0         0         0         0         0         0  ;
         0         0         0         0         0         0         0         0  ;
         0         0         0         0         0        -L         0         0  ;
         0         0         0         0         0         0         0         0  ;
         0         0         0         0         0         0         0         0  ;
    ];

F = [
    0   ;
    0   ;
    0   ;
    0   ;
    0   ;
    0   ;
    0   ;
    0   ;
    ];
disp('QUESTION 3 G and C matrices:')
disp(G)
disp(C_1)

time = 1;
timestep = 1000;
delta = time/timestep;

% -----------------------------------------------------------------------
% GAUSSIAN PULSE FUNCTION
% -----------------------------------------------------------------------
magnitude = 1;
std_dev = 0.03;
delay = 0.06;
% Set up the list as before
F_List = zeros(8,1,timestep);

current_storage = zeros(1,timestep);
% This time it's not a step function so we will need to populate input
% voltage a bit differently using the given sine function
for counter = 1:1:timestep
    F_List(8,1,counter) = exp(-((counter*delta - delay)/std_dev)^2);
    % Add noise to node 3
    F_List(3,1,counter) = -0.001*randn;
    % for histogram
    current_storage(counter) = F_List(3,1,counter);
end

V_List_1 = zeros(8,1,timestep);
V_List_2 = zeros(8,1,timestep);
V_List_3 = zeros(8,1,timestep);

% Set up a bunch of matrices containing the voltage value for each timestep
for counter = 2:1:timestep
    e_1 = C_1/delta + G;
    e_2 = C_2/delta + G;
    e_3 = C_3/delta + G;
    V_List_1(:,:,counter) = e_1\(C_1*V_List_1(:,:,counter-1)/delta +F_List(:,:,counter));
    V_List_2(:,:,counter) = e_2\(C_1*V_List_2(:,:,counter-1)/delta +F_List(:,:,counter));
    V_List_3(:,:,counter) = e_3\(C_1*V_List_3(:,:,counter-1)/delta +F_List(:,:,counter));
end

% Save our input and output voltages in new matrices for plotting
% Cn 1
V_Out_List_1(1,:) = V_List_1(5,1,:);
V_In_List_1(1,:) = V_List_1(1,1,:);

% Cn 2
V_Out_List_2(1,:) = V_List_2(5,1,:);
V_In_List_2(1,:) = V_List_2(1,1,:);
% Cn 3
V_Out_List_3(1,:) = V_List_3(5,1,:);
V_In_List_3(1,:) = V_List_3(1,1,:);


figure(12)
plot((1:timestep).*delta, V_Out_List_1(1,:))
hold on
plot((1:timestep).*delta, V_In_List_1(1,:))
title('Question 3 CN_1 = 0.00001 Gaussian Pulse Function with Added Noise Source')
xlabel('Time (s)')
ylabel('Voltage (v)')
grid on
legend('Output Voltage','Input Voltage')
hold off

figure(13)
histogram(current_storage)
title('Question 3 Noise Source Histogram')
grid on
xlabel('Current (A)')

% Fourier Transform
figure(14)
FF = abs(fftshift(fft(V_Out_List_1(1,:))));
plot(((1:length(FF))/timestep)-0.5,FF)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
xlim([-0.04 0.04])
title('Question 3 CN_1 = 0.00001 Gaussian Pulse Function (With Noise) Fourier Transform')

% CN 2 -----------------------------------------

figure(15)
plot((1:timestep).*delta, V_Out_List_2(1,:))
hold on
plot((1:timestep).*delta, V_In_List_2(1,:))
title('Question 3 CN_2 = 0.0001 Gaussian Pulse Function with Added Noise Source')
xlabel('Time (s)')
ylabel('Voltage (v)')
grid on
legend('Output Voltage','Input Voltage')
hold off

% Fourier Transform
figure(16)
FF = abs(fftshift(fft(V_Out_List_2(1,:))));
plot(((1:length(FF))/timestep)-0.5,FF)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
xlim([-0.04 0.04])
title('Question 3 CN_2 = 0.0001 Gaussian Pulse Function (With Noise) Fourier Transform')

% CN 3 -----------------------------------------

figure(17)
plot((1:timestep).*delta, V_Out_List_3(1,:))
hold on
plot((1:timestep).*delta, V_In_List_3(1,:))
title('Question 3 CN_3 = 0.01 Gaussian Pulse Function with Added Noise Source')
xlabel('Time (s)')
ylabel('Voltage (v)')
grid on
legend('Output Voltage','Input Voltage')
hold off

% Fourier Transform
figure(18)
FF = abs(fftshift(fft(V_Out_List_3(1,:))));
plot(((1:length(FF))/timestep)-0.5,FF)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
xlim([-0.04 0.04])
title('Question 3 CN_3 = 0.01 Gaussian Pulse Function (With Noise) Fourier Transform')

disp('QUESTION 3: We can see how increasing the value of Cn reduces the overall output of the circuit. With C3 = 0.01, the output is actually lower than the input.')
disp('QUESTION 3: As before, changing the timestep will have an effect on the accuracy of the simulation.')