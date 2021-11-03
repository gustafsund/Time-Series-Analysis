%% Prophet for comparison
addpath('C:\Users\gusta\OneDrive\Dokument\Skola\Tidsserier\matlabTSA20\data')
addpath('C:\Users\gusta\OneDrive\Dokument\Skola\Tidsserier\matlabTSA20\matlab')
addpath('C:\Users\gusta\OneDrive\Dokument\Skola\Tidsserier\data20\data')
close all 
clear all
clc

load 'fjarrvarme89.dat'

load 'fjarrvarme90.dat'

start = 300;
N = 8*24*7;
Nv = 2*7*24;
Nt = 7*24;

%% Prophet validation, k = 8
close all
prophet_valid = readtable('C:\Users\gusta\OneDrive\Dokument\Skola\Tidsserier\Prophet_valid2.csv','HeaderLines',1);

pwr = fjarrvarme90(start:start+N+Nv-1,2);
prophet_pwr_valid = table2array(prophet_valid(:,3));
prophet_pwr_valid = prophet_pwr_valid(1:end-1);
pwr = pwr(2:end);
plot(pwr,'b')
hold on
plot(prophet_pwr_valid,'r')
xline(N);
legend('blue = real data','red = Prophet prediction','vertical = cut of between model data and validation')
title('Prophet 8-step prediction on validation data set')

res_prophet_valid= pwr(N:end) - prophet_pwr_valid(N:end);
figure
triplePlot(res_prophet_valid,20,false);
mean(res_prophet_valid)
var(res_prophet_valid)

%% prophet validation, k = 1 
close all
prophet_valid = readtable('C:\Users\gusta\OneDrive\Dokument\Skola\Tidsserier\Prophet_valid3.csv','HeaderLines',1);

pwr = fjarrvarme90(start:start+N+Nv-1,2);
prophet_pwr_valid = table2array(prophet_valid(:,3));
prophet_pwr_valid = prophet_pwr_valid(1:end-1);
pwr = pwr(2:end);
plot(pwr,'b')
hold on
plot(prophet_pwr_valid,'r')
xline(N);
legend('blue = real data','red = Prophet prediction','vertical = cut of between model data and validation')
title('Prophet 1-step prediction on validation data set')

res_prophet_valid= pwr(N:end) - prophet_pwr_valid(N:end);
figure
triplePlot(res_prophet_valid,20,false);
mean(res_prophet_valid)
var(res_prophet_valid)


%% prophet test1
close all
prophet_test1 = readtable('C:\Users\gusta\OneDrive\Dokument\Skola\Tidsserier\Prophet_test1.csv','HeaderLines',1);
pwr = fjarrvarme90(start:start+N+Nv+Nt-1,2);
prophet_pwr_test1 = table2array(prophet_test1(:,3));
prophet_pwr_test1 = prophet_pwr_test1(1:end-1);
pwr = pwr(2:end);
plot(pwr,'b')
hold on
plot(prophet_pwr_test1,'r')
xline(N+Nv);
legend('blue = real data','red = Prophet prediction','vertical = cut off between model + validation data and test set 1')
title('Prophet 6-step prediction on test set 1')
res_prophet_test1 = pwr(N+Nv:end) - prophet_pwr_test1(N+Nv:end);
figure
triplePlot(res_prophet_test1,20,false);
mean(res_prophet_test1)
var(res_prophet_test1)