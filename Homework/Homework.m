clc
clear

Kp = input("Kp = ");

sq = sqrt(Kp);

s_1 = -sq + sq*sqrt((16/Kp) - 1);
s_2 = -sq - sq*sqrt((16/Kp) - 1);

disp(s_1)
disp(s_2)