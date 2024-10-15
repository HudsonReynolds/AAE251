%AAE 251 Fall 2024
%Homework 5, Question 2
%AAE251_HW5_2
%Authors: Hudson Reynolds & Preston Wright

%% Data:

% alpha list:

alpha = -2:2:18;

% cL and cD data:
cL = [0,0.2,0.42,0.63,0.85,1.08,1.28,1.43,1.56,1.62,1.57];

%cD data. Last two data points not given in charts.
cD = [0.0063, 0.0062, 0.0065, 0.007, 0.008, 0.0092, 0.0112, 0.0147, 0.0187, 0, 0];

L2D = cL ./ cD;

% perform a polynomial fitting to find the extrapolated values of cD.

% create a list of cL values with a finer resolution
cLList = linspace(0,2, 181);
alphaList = linspace(-2,18, 181);

% fit the cL v. cD to a fourth order polynomial
b = polyfit(cL(1:9),cD(1:9),4);

cdEquation = polyval(b, cLList);

L2DExtrapolated = cLList ./ cdEquation;

figure(1)
hold on
plot(cL(1:9),cD(1:9), '.', 'MarkerSize', 10, 'Color', 'Red')
plot(cLList, cdEquation, '--', 'linewidth', 1, 'Color', 'blue')
grid on;
title('Coefficient of Lift vs. Coefficient of Drag')
xlabel('Coefficient of Lift')
ylabel('Coefficient of Drag')
legend('Anderson Data', 'Best Fit Polynomial')
hold off

figure(2)
hold on
plot(alpha, L2D, '.', 'MarkerSize', 10)
plot(alphaList, L2DExtrapolated, '--')
xlim([-2.5,18])
xline(-2, '--', 'Zero Lift Angle of Attack', 'LabelVerticalAlignment','middle')
xline(16, '--', 'Critical Angle of Attack')
title('Lift to Drag Ratio as Function of Angle of Attack')
xlabel('Angle of Attack [deg]')
ylabel('Lift to Drag Ratio [-]')
legend('L/D Data', 'Best Fit Line', 'Location', 'Northwest')
grid on
hold off
