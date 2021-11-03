from fbprophet import Prophet
import pandas as pd
import matplotlib.pyplot as plt



df = pd.read_csv(r'c:\users\gusta\onedrive\dokument\skola\tidsserier\fjarrvarme90.csv') 
df['ds'] = pd.to_datetime(df['ds'])
start = 300
N = 8*7*24
Nv = 2*7*24
Nt = 7*24
mdf = df[start:start+N-1]
vdf = df[start:start + N + Nv -1]
#plt.plot(df['ds'],df['y'])

#dfm = df[300:300+8*7*24][:]

plt.plot(df['ds'], df['y'])
k = 1
predictions = df[start:start+N-1]

for i in range(start + N -k,start+N + Nv + 1 - k):
    m = Prophet()
    cdf = df[start:i]
    m.fit(cdf)
    future = m.make_future_dataframe(k,freq='H')
    forecast = m.predict(future)
    ind = len(forecast)-1-k-start
    predictions = predictions.append({'ds':forecast.at[ind,'ds'],'y':forecast.at[ind,'yhat']},ignore_index = True)

predictions.to_csv(r'C:\Users\gusta\OneDrive\Dokument\Skola\Tidsserier\Prophet_valid3.csv')
