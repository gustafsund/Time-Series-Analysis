%% A
close all
clear all
f1 = load("fjarrvarme89.dat");
f2 = load("fjarrvarme90.dat");
%% ;
subplot(3,2,1)
plot(f1(:,2))
%% hej
subplot(3,2,2)
plot(f2(:,2))
subplot(3,2,3)
plot(f1(:,3))
subplot(3,2,4)
plot(f2(:,3))
subplot(3,2,5)
plot(f1(:,4))
subplot(3,2,6)
plot(f2(:,4))
%% Select data
start = 300;
N=8*7*24;
l=N;
margin=30;
mdata = f2(start:start+l,:);
vdata=f2(start+l+1-margin:start+l+24*7*2,:); %336 length
t1data =f2(start+l+24*7*2+1-margin:start+l+24*7*3,:); %168 length
startT2=length(f2)-24*7;
endT2=length(f2);
t2data =f2(startT2-margin:endT2,:);
maxLag=floor(N/4);
airtemp = mdata(:,3);
power = mdata(:,2);
close all
plot(airtemp)
%% PW
airtemp=airtemp-mean(airtemp);
power=power-mean(power);
triplePlot(airtemp,100)
z=iddata(power,airtemp);

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
powerpw= resid(model_1, power).y;
triplePlot(powerpw, 50)
%% Compute CCF
CCF(powerpw, ATpw);
% fråga om denna!!!

%% model B and A2 DETERMINE A2 AND B
A2 = [1];
B = [1 0 0];
Mi = idpoly ([1] ,[B] ,[] ,[] ,[A2]); 
zpw = iddata(powerpw,ATpw);
Mba2 = pem(zpw,Mi); 
present(Mba2)
vhat = resid(Mba2,zpw);
v = vhat.y;
triplePlot(v,50);

%% CCF
CCF(ATpw,v)

%% Whiteness
%whitenessTest(v)
%% X
xhat = resid(Mba2,z);
x = xhat.y;
triplePlot(x,50)

%% Find A1 and C
A1 = [1 zeros(1,25)];
C= [1];
model_init = idpoly([A1],[],[C]);
model_init.Structure.a.Free = [1 1 1 zeros(1,20) 1 1 1];
data =iddata(x);
model_x = pem(data,model_init);
res = resid(model_x, x).y;
triplePlot(res, 100)
whitenessTest(res)

%% Estimate all
A1=[1 zeros(1,25)];
A2=[1];
B = [1 0 0];
C=[1];
Mi = idpoly (1 ,B,C,A1,A2);
Mi.Structure.d.Free = [1 1 1 zeros(1,20) 1 1 1];
z_ = iddata(power,airtemp);
MboxJ = pem(z_ ,Mi); 
present(MboxJ)
ehat = resid(MboxJ,z_);
triplePlot(ehat.y,80)
whitenessTest(ehat.y)
% Passes 1
%% Compute CCF for e and u
CCF(airtemp,ehat.y)

%% Test model on validation data

valid=iddata(vdata(:,2)-mean(mdata(:,2)), vdata(:,3)-mean(vdata(:,3)));
ehatv = resid(MboxJ,valid);
triplePlot(ehatv.y,80)
whitenessTest(ehatv.y)

%% Predict on validation
A=MboxJ.d;
C=MboxJ.c;
B=conv(MboxJ.b, A);
k=1;
[Fk, Gk] = polydiv( C, A, k );
BF = conv(B,Fk);

%%
%[Fku, Gku] = polydiv( BF, C, k );
%uhat=myFilter(Gku, C, u) ;
airtempV=valid.u;
uhat=myFilter(BF,C,airtempV);
close all
plot(airtempV(100:200), 'r')
hold on
plot(uhat(100:200), 'b')
%%
powerV=valid.y;
yhat=myFilter(Gk, C, powerV);
% Fotnot 3: initiate predictions well before the beginning...?
%%
close all
plot(uhat, 'b')
hold on
plot(yhat, 'r')
%%
ytilde=yhat+uhat;
ytilde=ytilde(margin+1:end);
powerV=powerV(margin+1:end);
close all
plot(ytilde(1:199), 'b')
hold on
plot(powerV(1:199), 'r-')
%% Prediction error
res=powerV-ytilde;
triplePlot(res, 80)
whitenessTest(res)
var(res)

%% naive predictor - validation data
close all
mpower = mean(mdata(:,2));
powerV = powerV+mpower
naive_valid = [mpower*ones(k,1); powerV(1:end-k)];
plot(powerV(margin+1:end),'b')
hold on 
plot(naive_valid(margin+1:end),'r')
legend('blue = realization','red = Naive prediction, validation data')

res = naive_valid(margin+1:end) - powerV(margin+1:end);
figure 
triplePlot(res,100,false);
if k == 1
    figure
    whitenessTest(res)
end
mean(res)
var(res)

%% Predict on test 1
test1=iddata(t1data(:,2)-mean(mdata(:,2)), t1data(:,3)-mean(t1data(:,3)));
k=6;
[Fk, Gk] = polydiv( C, A, k );
BF = conv(B,Fk);

%%
airtempT1=test1.u;
uhat=myFilter(BF,1,airtempT1);
close all
plot(airtempT1(50:149), 'r')
hold on
plot(uhat(50:149), 'b')
%%
powerT1=test1.y;
yhat=myFilter(Gk, C, powerT1);
% Fotnot 3: initiate predictions well before the beginning...?
%%
close all
plot(uhat, 'b')
hold on
plot(yhat, 'r')
%%
ytilde=yhat+uhat;
ytilde=ytilde(margin+1:end);
powerT1=powerT1(margin+1:end);
close all
plot(ytilde, 'r')
hold on
plot(powerT1, 'b-')
%% Prediction error
res=powerT1-ytilde;
triplePlot(res, 80)
var(res)

%% naive predictor - test data 1
close all
powerT1 = powerT1 + mpower
naive_test1= [mpower*ones(k,1); powerT1(1:end-k)];
plot(powerT1(margin+1:end),'b')
hold on 
plot(naive_test1(margin+1:end),'r')
legend('blue = realization','red = Naive prediction, validation data')

res = naive_test1(margin+1:end) - powerT1(margin+1:end);
figure 
triplePlot(res,100,false);
if k == 1
    figure
    whitenessTest(res)
end
mean(res)
var(res)



%% Predict on test 2
test2=iddata(t2data(:,2)-mean(mdata(:,2)), t2data(:,3)-mean(t2data(:,3)));
k=6;
[Fk, Gk] = polydiv( C, A, k );
BF = conv(B,Fk);

%%
airtempT2=test2.u;
uhat=myFilter(BF,1,airtempT2);

%%
powerT2=test2.y;
yhat=myFilter(Gk, C, powerT2);
%%
ytilde=yhat+uhat;
ytilde=ytilde(margin+1:end);
powerT2=powerT2(margin+1:end);
close all
plot(ytilde, 'r')
hold on
plot(powerT2, 'b-')
%% Prediction error
res=powerT2-ytilde;
triplePlot(res, 80)
var(res)



%% naive predictor - test data 1
close all
powerT2 = powerT2 + mpower
naive_test2= [mpower*ones(k,1); powerT2(1:end-k)];
plot(powerT2(margin+1:end),'b')
hold on 
plot(naive_test2(margin+1:end),'r')
legend('blue = realization','red = Naive prediction, validation data')

res = naive_test2(margin+1:end) - powerT2(margin+1:end);
figure 
triplePlot(res,100,false);
if k == 1
    figure
    whitenessTest(res)
end
mean(res)
var(res)


