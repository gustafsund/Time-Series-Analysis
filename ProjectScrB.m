%% B
close all
clear all
f1 = load("fjarrvarme89.dat");
f2 = load("fjarrvarme90.dat");

%% Select data
start = 300;
N=8*7*24;
l=N;
margin=50;
mdata = f2(start:start+l,:);
vdata=f2(start+l+1-margin:start+l+24*7*2,:); %336 length
t1data =f2(start+l+24*7*2+1-margin:start+l+24*7*3,:); %168 length
startT2=length(f2)-24*7;
endT2=length(f2);
t2data =f2(startT2-margin:endT2,:);
maxLag=floor(N/4);
airtemp = mdata(:,3);
power = mdata(:,2);
water = mdata(:,4);
close all
plot(airtemp)
%% PW
airtemp=airtemp-mean(airtemp);
power=power-mean(power);
water=water-mean(water);
triplePlot(water,150)
%%
close all
CCF(power, water) % more correlated
%%
CCF(power, airtemp)
%%
triplePlot(water,50)

%% Identify A3 and C3
model_init = idpoly([1 zeros(1,25)],[],[1]); 
model_init.Structure.a.Free= [1 1 1 zeros(1,20) 1 0 1];
model_1 = pem(water,model_init);
Wpw = resid(model_1, water).y;
triplePlot(Wpw, 100)
whitenessTest(Wpw)
present(model_1)
A3=model_1.a;
C3=model_1.c;
% Passes 2
%% PW load
powerpw= resid(model_1, power).y;
triplePlot(powerpw, 50)
%% Compute CCF
CCF(powerpw, Wpw);
%% model B and A2 DETERMINE A2 AND B
A21 = [1];
B1 = [1 0 0];
Mi = idpoly ([1] ,[B1] ,[] ,[] ,[A21]); 
zpw = iddata(powerpw,Wpw);
Mba2 = pem(zpw,Mi); 
present(Mba2)
vhat = resid(Mba2,zpw);
v = vhat.y;
triplePlot(v,80);
A21 = Mba2.f;
B1 = Mba2.b;
%% CCF
CCF(Wpw,v) % not necc white
%% Model next input
z=iddata(power, water);
y2=resid(z,Mba2).y;
triplePlot(y2,80);

%% Identify A3 and C3
model_init = idpoly([1 0 0],[],[1]); 
ATdata =iddata(airtemp);
model_1 = pem(ATdata,model_init);
ATpw = resid(model_1, airtemp).y;
triplePlot(ATpw, 100)
whitenessTest(ATpw)
present(model_1)
A3=model_1.a;
C3=model_1.c;

%% PW load
y2pw= resid(model_1, y2).y;
triplePlot(y2pw, 50)
%% Compute CCF
CCF(y2pw, ATpw);
%% model B and A2 DETERMINE A2 AND B
A22 = [1];
B2 = [1];
Mi = idpoly ([1] ,[B2] ,[] ,[] ,[A22]); 
zpw = iddata(y2pw,ATpw);
Mba2 = pem(zpw,Mi); 
present(Mba2)
vhat = resid(Mba2,zpw);
v = vhat.y;
triplePlot(v,50);
A22 = Mba2.f;
B2 = Mba2.b;

%% CCF
CCF(ATpw,v)

%% etilde
z = iddata(y2, airtemp);
etilde = resid(Mba2,z).y;
triplePlot(etilde,50)

%% Find A1 and C
A1 = [1 zeros(1,25)];
C= [1];
model_init = idpoly([A1],[],[C]);
model_init.Structure.a.Free = [1 1 1 zeros(1,8) 1 1 zeros(1,10) 1 1 1];
data =iddata(etilde);
model_x = pem(data,model_init);
present(model_x)
res = resid(model_x, etilde).y;
triplePlot(res, 100)
whitenessTest(res)
A1=model_x.a;
C=model_x.c;


%% Estimate all
A21c=mat2cell(A21,1,1);
A22c=mat2cell(A22,1,1);
A2 =[A21c A22c];
B1c=mat2cell(B1,1,3);
B2c=mat2cell(B2,1,1);
B =[B1c B2c];

%%
Mi = idpoly (1 ,B,C,A1,A2);
Mi.Structure.d.Free = [1 1 1 zeros(1,8) 1 1 zeros(1,10) 1 1 1];
z = iddata(power,[water airtemp]);
MboxJ = pem(z ,Mi); 
present(MboxJ)
ehat = resid(MboxJ,z);
triplePlot(ehat.y,80)
whitenessTest(ehat.y)

%% Test model on validation data

valid=iddata(vdata(:,2)-mean(mdata(:,2)), [vdata(:,4)-mean(vdata(:,4)) vdata(:,3)-mean(vdata(:,3))]);
ehatv = resid(MboxJ,valid);
triplePlot(ehatv.y,80)
whitenessTest(ehatv.y)

%% Predict on validation data
A=MboxJ.d;
C=MboxJ.c;
B1=conv(cell2mat(MboxJ.b(1)), A);
B2=conv(cell2mat(MboxJ.b(2)), A);
k=1;
[Fk, Gk] = polydiv( C, A, k );
B1F=conv(B1, Fk);
B2F=conv(B2, Fk);

%%
%[Fku, Gku] = polydiv( BF, C, k );
%uhat=myFilter(Gku, C, u) ;
airtempV=valid.u(:,2);
waterV=valid.u(:,1);
airToPower=myFilter(B2F,C,airtempV);
waterToPower=myFilter(B1F,C,waterV);
close all
plot(airtempV(100:200), 'r')
hold on
plot(airToPower(100:200), 'b')
%%
close all
plot(waterV(100:200), 'r')
hold on
plot(waterToPower(100:200), 'b')
%%
powerV=valid.y;
powerToPower=myFilter(Gk, C, powerV);
%%
close all
plot(waterToPower, 'b')
hold on
plot(powerToPower, 'r')
plot(airToPower, 'g')
%%
prediction=powerToPower+airToPower+waterToPower;
prediction=prediction(margin+1:end);
powerVNoMarg=powerV(margin+1:end);
close all
plot(prediction(1:199), 'b')
hold on
plot(powerVNoMarg(1:199), 'r-')
%% Addmean
predm=prediction+mean(mdata(:,2));
powerVm=powerVNoMarg+mean(mdata(:,2));
close all
plot(predm(1:199), 'b')
hold on
plot(powerVm(1:199), 'r-')
%% Prediction error
res=powerVNoMarg-prediction;
triplePlot(res, 80)
whitenessTest(res)
var(res)
% pred k=8 ser inte ut som en MA(7). Vad göra?

%% Predict on test1
test1=iddata(t1data(:,2)-mean(mdata(:,2)), [t1data(:,4)-mean(t1data(:,4)) t1data(:,3)-mean(t1data(:,3))]);
k=6;
[Fk, Gk] = polydiv( C, A, k );
B1F = conv(B1, Fk);
B2F = conv(B2, Fk);

%%
%[Fku, Gku] = polydiv( BF, C, k );
%uhat=myFilter(Gku, C, u) ;
airtempT=test1.u(:,2);
waterT=test1.u(:,1);
airToPower=myFilter(B2F,C,airtempT);
waterToPower=myFilter(B1F,C,waterT);
close all
plot(airtempT(50:149), 'r')
hold on
plot(airToPower(50:149), 'b')
%%
close all
plot(waterT(50:149), 'r')
hold on
plot(waterToPower(50:149), 'b')
%%
powerT=test1.y;
powerToPower=myFilter(Gk, C, powerT);
%%
close all
plot(waterToPower, 'b')
hold on
plot(powerToPower, 'r')
plot(airToPower, 'g')
%%
prediction=powerToPower+airToPower+waterToPower;
prediction=prediction(margin+1:end);
powerTNoMarg=powerT(margin+1:end);
close all
plot(prediction(50:149), 'b')
hold on
plot(powerTNoMarg(50:149), 'r-')
hold off
%% Addmean
predm=prediction+mean(mdata(:,2));
powerTm=powerTNoMarg+mean(mdata(:,2));
plot(predm(50:149), 'b')
hold on
plot(powerTm(50:149), 'r-')
%% Prediction error
res=powerTNoMarg-prediction;
triplePlot(res, 80)
whitenessTest(res)
var(res)
% pred k=6fukt. K=1 ger white!!

%% Predict on test2
test2=iddata(t2data(:,2)-mean(mdata(:,2)), [t2data(:,4)-mean(t2data(:,4)) t2data(:,3)-mean(t2data(:,3))]);
k=6;
[Fk, Gk] = polydiv( C, A, k );
B1F = conv(B1, Fk);
B2F = conv(B2, Fk);
%%
%[Fku, Gku] = polydiv( BF, C, k );
%uhat=myFilter(Gku, C, u) ;
airtempT=test2.u(:,2);
waterT=test2.u(:,1);
airToPower=myFilter(B2F,C,airtempT);
waterToPower=myFilter(B1F,C,waterT);
close all
plot(airtempT(50:149), 'r')
hold on
plot(airToPower(50:149), 'b')
%%
close all
plot(waterT(50:149), 'r')
hold on
plot(waterToPower(50:149), 'b')
%%
powerT=test2.y;
powerToPower=myFilter(Gk, C, powerT);
%%
close all
plot(waterToPower, 'b')
hold on
plot(powerToPower, 'r')
plot(airToPower, 'g')
%%
prediction=powerToPower+airToPower+waterToPower;
prediction=prediction(margin+1:end);
powerTNoMarg=powerT(margin+1:end);
close all
plot(prediction(50:149), 'b')
hold on
plot(powerTNoMarg(50:149), 'r-')
hold off
%% Addmean
predm=prediction+mean(mdata(:,2));
powerTm=powerTNoMarg+mean(mdata(:,2));
plot(predm(50:149), 'b')
hold on
plot(powerTm(50:149), 'r-')
%% Prediction error
res=powerTNoMarg-prediction;
triplePlot(res, 80)
whitenessTest(res)
var(res)
% pred k=6 fukt. K=1 ger white!!



