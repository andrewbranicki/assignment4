%% ELEC 4700 Assignment 4 Question 2
% Andrew Branicki 100973961

disp('Question 2 a) We can see from the frequency response plot in the MNPA that this is a low pass filter circuit.')
disp('As the frequency increases, the inductor will resist the change in current (it will start to block current) and the capacitor will create a parallel resistance with the resistor.')
disp('Question 2 b) The expected frequency response would be a constant gain at low frequencies, and then a -3 dB cutoff at the specific cutoff frequency (fc = 1/(2piRC)) with a drop in gain after that.')

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
    0   ;
    ];

disp('QUESTION 2 G and C matrices:')
disp(G)
disp(C)

time = 1;
timestep = 1000;
delta = time/timestep;

% -----------------------------------------------------------------------
% STEP FUNCTION
% -----------------------------------------------------------------------
F_List = zeros(8,1,timestep);
% At 0.03s we need to set it to 1
F_List(8,1,30:timestep) = 1;

V_List = zeros(8,1,timestep);

% Set up a bunch of matrices containing the voltage value for each timestep
for counter = 2:1:timestep
    e = C/delta + G;
    V_List(:,:,counter) = e\(C*V_List(:,:,counter-1)/delta +F_List(:,:,counter));
end

% Save our input and output voltages in new matrices for plotting
V_Out_List_1(1,:) = V_List(5,1,:);
V_In_List_1(1,:) = V_List(1,1,:);

% Plot it up
figure(4)
plot((1:timestep).*delta, V_Out_List_1(1,:))
hold on
plot((1:timestep).*delta, V_In_List_1(1,:))
title('Question 2 Step Function')
xlabel('Time (s)')
ylabel('Voltage (v)')
grid on
legend('Output Voltage','Input Voltage')

% Fourier Transform
figure(5)
FF = abs(fftshift(fft(V_Out_List_1(1,:))));
plot(((1:length(FF))/timestep)-0.5,FF)
xlim([-0.03 0.03])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
title('Question 2 Step Function Fourier Transform')

% -----------------------------------------------------------------------
% SINE FUNCTION
% -----------------------------------------------------------------------
f1 = 1/0.03;
f2 = 1/0.1;
f3 = 1/0.5;
% Set up the list as before
F_List_1 = zeros(8,1,timestep);
F_List_2 = zeros(8,1,timestep);
F_List_3 = zeros(8,1,timestep);

% This time it's not a step function so we will need to populate input
% voltage a bit differently using the given sine function
for counter = 1:1:timestep
    F_List_1(8,1,counter) = sin(2*pi*f1*counter*delta);
    F_List_2(8,1,counter) = sin(2*pi*f2*counter*delta);
    F_List_3(8,1,counter) = sin(2*pi*f3*counter*delta);
end

V_List_1 = zeros(8,1,timestep);
V_List_2 = zeros(8,1,timestep);
V_List_3 = zeros(8,1,timestep);

% Set up a bunch of matrices containing the voltage value for each timestep
for counter = 2:1:timestep
    e = C/delta + G;
    V_List_1(:,:,counter) = e\(C*V_List_1(:,:,counter-1)/delta +F_List_1(:,:,counter));
    V_List_2(:,:,counter) = e\(C*V_List_2(:,:,counter-1)/delta +F_List_2(:,:,counter));
    V_List_3(:,:,counter) = e\(C*V_List_3(:,:,counter-1)/delta +F_List_3(:,:,counter));
end

% Save our input and output voltages in new matrices for plotting
% Freq 1
V_Out_List_1(1,:) = V_List_1(5,1,:);
V_In_List_1(1,:) = V_List_1(1,1,:);
% Freq 2
V_Out_List_2(1,:) = V_List_2(5,1,:);
V_In_List_2(1,:) = V_List_2(1,1,:);
% Freq 3
V_Out_List_3(1,:) = V_List_3(5,1,:);
V_In_List_3(1,:) = V_List_3(1,1,:);

% Plot it up
figure(6)
plot((1:timestep).*delta, V_Out_List_1(1,:))
hold on
plot((1:timestep).*delta, V_In_List_1(1,:))
title('Question 2 Sine Function f = 1/0.03')
xlabel('Time (s)')
ylabel('Voltage (v)')
grid on
legend('Output Voltage','Input Voltage')
hold off

% Fourier Transform
figure(7)
FF = abs(fftshift(fft(V_Out_List_1(1,:))));
plot(((1:length(FF))/timestep)-0.5,FF)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
xlim([-0.1 0.1])
title('Question 2 Sine Function f = 1/0.03 Fourier Transform')

figure(8)
plot((1:timestep).*delta, V_Out_List_2(1,:))
hold on
plot((1:timestep).*delta, V_In_List_2(1,:))
title('Question 2 Sine Function f = 1/0.1')
xlabel('Time (s)')
ylabel('Voltage (v)')
grid on
legend('Output Voltage','Input Voltage')
hold off

figure(9)
plot((1:timestep).*delta, V_Out_List_3(1,:))
hold on
plot((1:timestep).*delta, V_In_List_3(1,:))
title('Question 2 Sine Function f = 1/0.5')
xlabel('Time (s)')
ylabel('Voltage (v)')
grid on
legend('Output Voltage','Input Voltage')
hold off

disp('Question 2: SINE WAVE INPUT: We can see that at slower frequencies the gain of the circuit increases. We are able to get a larger voltage output.')

% -----------------------------------------------------------------------
% GAUSSIAN PULSE FUNCTION
% -----------------------------------------------------------------------
magnitude = 1;
std_dev = 0.03;
delay = 0.06;
% Set up the list as before
F_List = zeros(8,1,timestep);

% This time it's not a step function so we will need to populate input
% voltage a bit differently using the given sine function
for counter = 1:1:timestep
    F_List(8,1,counter) = exp(-((counter*delta - delay)/std_dev)^2);
end

V_List = zeros(8,1,timestep);

% Set up a bunch of matrices containing the voltage value for each timestep
for counter = 2:1:timestep
    e = C/delta + G;
    V_List(:,:,counter) = e\(C*V_List(:,:,counter-1)/delta +F_List(:,:,counter));
end

% Save our input and output voltages in new matrices for plotting
% Freq 1
V_Out_List(1,:) = V_List(5,1,:);
V_In_List(1,:) = V_List(1,1,:);

figure(10)
plot((1:timestep).*delta, V_Out_List(1,:))
hold on
plot((1:timestep).*delta, V_In_List(1,:))
title('Question 2 Gaussian Pulse Function')
xlabel('Time (s)')
ylabel('Voltage (v)')
grid on
legend('Output Voltage','Input Voltage')
hold off

% Fourier Transform
figure(11)
FF = abs(fftshift(fft(V_Out_List(1,:))));
plot(((1:length(FF))/timestep)-0.5,FF)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
xlim([-0.03 0.03])
title('Question 2 Gaussian Pulse Function Fourier Transform')

disp('Question 2: Changing the time step changes how accurate the simulation is. With a larger time step we have less accurate simulations since we are getting less detail overall.')
