%% Project c)
addpath('C:\Users\gusta\OneDrive\Dokument\Skola\Tidsserier\matlabTSA20\data')
addpath('C:\Users\gusta\OneDrive\Dokument\Skola\Tidsserier\matlabTSA20\matlab')
addpath('C:\Users\gusta\OneDrive\Dokument\Skola\Tidsserier\data20\data')
close all
clear all
clc
f1 = load('fjarrvarme89.dat');
f2 = load('fjarrvarme90.dat');
load 'MBJ.mat'

% csvwrite('fjarrvarme89.csv',f1);
% csvwrite('fjarrvarme90.csv',f2)
Nm = 8*7*24;
start = 300;
modelRange = (start:start+Nm-1);
Nv = 2*7*24;
validRange = (modelRange(end)+1:modelRange(end)+Nv);
Nt = 7*24;
test1Range = (start+Nm+Nv:start+Nm+Nv+Nt-1);
test2Range = (length(f2)-Nt:length(f2));

power = f2(:,2);
% power = power-mean(power);
airtemp = f2(:,3);
airtemp = airtemp - mean(airtemp);
water = f2(:,4);
water = water - mean(water);

%% Kalman setup
close all
y = power;

N = length(y);
temp =MboxJ.d((MboxJ.d~=0));
AB1 = conv(MboxJ.d,cell2mat(MboxJ.b(1)));
AB2 = conv(MboxJ.d,cell2mat(MboxJ.b(2)));
y_shifts = find(MboxJ.d~=0) -1;
y_shifts = y_shifts(2:end);
b1_shifts = find(AB1~=0) - 1;
b2_shifts = find(AB2 ~=0) -1;

xtt1 = [temp(2:end) AB1(AB1~=0) AB2(AB2~=0)]';
states = length(xtt1);


%% Kalman filter
close all
A = eye(states);

Re = 0.1*eye(states);
Rw = 0.01;

Rxx1 = 0.1*eye(states);


xsave = zeros(states,N);
k = 1;
ysave = zeros(N,1);

for t = 30:N-k
C = [-y(t-y_shifts)' water(t-b1_shifts)' airtemp(t-b2_shifts)'];
Ryy= C*Rxx1*C' + Rw;
Kt=Rxx1*C'/Ryy;
xtt= xtt1 + Kt*(y(t) - C*xtt1);
Rxx=(eye(length(Rxx1))-Kt*C)*Rxx1;

xsave(:,t) = xtt1;
Rxx1= A*Rxx*A' + Re;
xtt1= A*xtt;    
y_temp = y(1:t);
    
    for p = 1:k
    Ct = [-y(t+p-y_shifts)' water(t+p-b1_shifts)' airtemp(t+p-b2_shifts)'];    
    y_temp(t+p) = Ct*xtt;
    
    end
ysave(t+k) = Ct*xtt;


end
subplot(2,1,1)
plot(y(50:end),'-b')
hold on
plot(ysave(50:end),'-r')
xline(modelRange(1));
xline(modelRange(end));
xline(validRange(end));
xline(test1Range(end));
xline(test2Range(end));
subplot(2,1,2)
plot(xsave')
legend

%% Residuals
close all
res = y(validRange) - ysave(validRange);
% res = y-ysave;
triplePlot(res,100,false);
mean(res)
var(res)
if k == 1
    figure
    whitenessTest(res);
end

%% Kalman - eliminating parameters one by one to see how they affect residual variances
close all
y = power;

N = length(y);
temp =MboxJ.d((MboxJ.d~=0));
AB1 = conv(MboxJ.d,cell2mat(MboxJ.b(1)));
AB2 = conv(MboxJ.d,cell2mat(MboxJ.b(2)));
y_shifts = find(MboxJ.d~=0) -1;
y_shifts = y_shifts(2:end);
b1_shifts = find(AB1~=0) - 1;
b2_shifts = find(AB2 ~=0) -1;

states=29;
varvec=zeros(1,states);


%% Loop
close all
% Outer loop
for j=1:states
    A = eye(states-1);
    Re = 0.1*eye(states-1);
    Rw = 0.01;
    Rxx1 = 0.1*eye(states-1);
    
    
    xsave = zeros(states-1,N);
    k = 6;
    ysave = zeros(N,1);
    xtt1 = [temp(2:end) AB1(AB1~=0) AB2(AB2~=0)]';
    xtt1(j)=[];
    % Kalman loop
    for t = 30:N-k
        C = [-y(t-y_shifts)' water(t-b1_shifts)' airtemp(t-b2_shifts)'];
        C(j)=[];    
        
        Ryy= C*Rxx1*C' + Rw;
        Kt=Rxx1*C'/Ryy;
        xtt= xtt1 + Kt*(y(t) - C*xtt1);
        Rxx=(eye(length(Rxx1))-Kt*C)*Rxx1;
        %Save
        xsave(:,t) = xtt1;
        Rxx1= A*Rxx*A' + Re;
        xtt1= A*xtt;    
        y_temp = y(1:t);
            
            for p = 1:k
                Ct = [-y(t+p-y_shifts)' water(t+p-b1_shifts)' airtemp(t+p-b2_shifts)'];    
                Ct(j)=[]; 
                y_temp(t+p) = Ct*xtt;
            
            end
        ysave(t+k) = Ct*xtt;
        
    end
    % Variance vectors for when parameter 1,2,3... are removed
    %varvec(j)=var(ysave(validRange)-y(validRange))
    %varvec(j)=var(ysave(test1Range)-y(test1Range))
    varvec(j)=var(ysave(test2Range)-y(test2Range))
end

%% Kalman setup - Reduced model
close all
y = power;

N = length(y);
temp =MboxJ.d((MboxJ.d~=0));
AB1 = conv(MboxJ.d,cell2mat(MboxJ.b(1)));
AB2 = conv(MboxJ.d,cell2mat(MboxJ.b(2)));
y_shifts = find(MboxJ.d~=0) -1;
y_shifts = y_shifts(2:end);
b1_shifts = find(AB1~=0) - 1;
b2_shifts = find(AB2 ~=0) -1;

xtt1 = [temp(2:end) AB1(AB1~=0) AB2(AB2~=0)]';
xtt1 = [temp(2:end) AB1(AB1~=0) AB2(AB2~=0)]';
% Remove unnecessary parameters
index=[12 18 20 27];
xtt1(index)=[];
states = length(xtt1);


%% Kalman filter - Reduced model
close all
A = eye(states);

Re = 0.1*eye(states);
Rw = 0.01;

Rxx1 = 0.1*eye(states);


xsave = zeros(states,N);
k = 1;
ysave = zeros(N,1);

for t = 30:N-k
C = [-y(t-y_shifts)' water(t-b1_shifts)' airtemp(t-b2_shifts)'];
C(index)=[];
Ryy= C*Rxx1*C' + Rw;
Kt=Rxx1*C'/Ryy;
xtt= xtt1 + Kt*(y(t) - C*xtt1);
Rxx=(eye(length(Rxx1))-Kt*C)*Rxx1;

xsave(:,t) = xtt1;
Rxx1= A*Rxx*A' + Re;
xtt1= A*xtt;    
y_temp = y(1:t);
    
    for p = 1:k
    Ct = [-y(t+p-y_shifts)' water(t+p-b1_shifts)' airtemp(t+p-b2_shifts)'];    
   Ct(index)=[];
    y_temp(t+p) = Ct*xtt;
    
    end
ysave(t+k) = Ct*xtt;


end
subplot(2,1,1)
plot(y(50:end),'-b')
hold on
plot(ysave(50:end),'-r')
xline(modelRange(1));
xline(modelRange(end));
xline(validRange(end));
xline(test1Range(end));
xline(test2Range(end));
subplot(2,1,2)
plot(xsave')
legend
res = y(test1Range) - ysave(test1Range);
var(res)

%% Residuals
close all
res = y(test1Range) - ysave(test1Range);
% res = y-ysave;
triplePlot(res,100,false);
mean(res)
var(res)
if k == 1
    figure
    whitenessTest(res);
end