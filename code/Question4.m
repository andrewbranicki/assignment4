%% ELEC 4700 Assignment 4 Question 4
% Andrew Branicki 100973961

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


time = 1;
timestep = 1000;
delta = time/timestep;

disp('QUESTION 4: In this case we can see that V is modelled by a cubic equation. We would need the values for beta and gamma, and then we can take the roots of this function and then use the values of the roots in our MNA matrices.')

