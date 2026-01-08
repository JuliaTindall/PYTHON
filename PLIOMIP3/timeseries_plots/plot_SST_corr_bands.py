#python program
#
# this program will plot the salinty from different bands and see what leads
# what
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import sys

expt = 'xpsig'
band1='90S-60S'
band2='60N-90N'

corr_yearstart =100
corr_yearend = 1800

# Function to read data from a file, skipping the first line
def read_data(filename,bandreq):
    band_index={'30S-30N':1,'30N-60N':2,'60N-90N':3,'30N-90N':4,
                '60S-30S':5,'90S-60S':6,'90S-30S':7}
    ix=band_index.get(bandreq)

    x_vals = []
    y_vals = []
    with open(filename, 'r') as file:
        next(file)  # Skip the title line
        for line in file:
            try:
                values = line.strip().split(',')
                x_vals.append(float(values[0]))
                y_vals.append(float(values[ix]))
             
            except ValueError:
                print(f"Skipping invalid line in {filename}: {line.strip()}")
    return x_vals, y_vals

# Function to compute running mean**
def running_mean(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')


##################################################################
plt.figure(figsize=(12, 9))
period = {'xpsid':'LP','xpsij':'LP490','xpsie':'EP400','xpsig':'EP',
          'xpsic':'PI'}
axis_label = {'salinity':'psu','AABW':'SV'}

#plot data from all the files
if expt == 'xpsij':
    startyear=1991
else:
    startyear=12

labels=[]
filename = ('/home/earjcti/um/' + expt + '/timeseries/SST_bands_' +
            expt + '_' + str(startyear) + '_2999.tex')

xband1, yband1 = read_data(filename,band1)
xband2, yband2= read_data(filename,band2)

# Plot the data
plt.subplot(211)

line1=plt.plot(xband1, [value - yband1[0] for value in yband1], label=band1, color='blue' )
line2=plt.plot(xband2, [value - yband2[0] for value in yband2], label=band2, color='orange' )
        
plt.title('Temperature trend: '+expt + ': ' +  period.get(expt))
plt.xlabel('year')
plt.ylabel('degC')
plt.legend(loc='lower right')
plt.grid(True)


#do a 30 year running mean
window=100
yband1_mean=running_mean(yband1,window)
xband1_mean=xband1[window-1:]
yband2_mean=running_mean(yband2,window)
xband2_mean=xband2[window-1:]

# Plot the data
plt.subplot(212)

plt.plot(xband1_mean, [value - yband1_mean[20] for value in yband1_mean], label=band1, color='blue' )
plt.plot(xband2_mean, [value - yband2_mean[20] for value in yband2_mean], label=band2,color='orange' )
        
plt.title(expt + '(30 year mean) : ' +  period.get(expt))

plt.xlabel('year')
plt.ylabel('degC')
plt.legend(loc='lower right')
plt.grid(True)
plt.xlim(corr_yearstart,corr_yearend)
plt.grid(True)
#plt.show()
plt.close()

##################################################################
# do a lag correlation on the raw data.  To avoid too high correlations
# only use xvalues between year start and yearend
print('LAG CORRELATIONS RAW DATA')
print('--------------------------------------')

xband1_arr = np.array(xband1)
mask=(xband1_arr >= corr_yearstart) & (xband1_arr <=corr_yearend)
yband1_sub = np.array(yband1)[mask]
yband2_sub=np.array(yband2)[mask]

yband1_series = pd.Series(yband1_sub)
yband2_series = pd.Series(yband2_sub)

# Choose the maximum lag to test (in time steps)
max_lag = 30  # 30 years

# Compute correlations for lags from -max_lag to +max_lag
lags = range(-max_lag, max_lag + 1)
lag_corrs = [yband1_series.corr(yband2_series.shift(lag)) for lag in lags]

# Print results
for lag, corr in zip(lags, lag_corrs):
    if lag < 0:
        relation = f"{band1} leads {band2} by {lag} steps"
    elif lag > 0:
        relation = f"{band2} leads {band1} by {abs(lag)} steps"
    else:
        relation = "No lead/lag (simultaneous)"
        
    r_squared = corr**2 if corr is not None else np.nan
    if r_squared > 0.2:
        print(f"{relation}: r = {corr:.3f}, R² = {r_squared:.3f}")


##################################################################
# do a lag correlation on the 30 year averaged data
print('LAG CORRELATIONS 30 yr avg DATA ')
print('--------------------------------------')
xband1_arr = np.array(xband1_mean)
mask=(xband1_arr >= corr_yearstart) & (xband1_arr <=corr_yearend)
yband1_sub = np.array(yband1_mean)[mask]
yband2_sub=np.array(yband2_mean)[mask]
x_sub = xband1_arr[mask]
print(np.shape(yband1_sub))
print(np.shape(yband1))

yband1_series = pd.Series(yband1_sub)
yband2_series = pd.Series(yband2_sub)

# Choose the maximum lag to test (in time steps)
max_lag = 50  # 30 years


# Compute correlations for lags from -max_lag to +max_lag on 30 year mean data
lags = range(-max_lag, max_lag + 1)
lag_corrs = [yband1_series.corr(yband2_series.shift(lag)) for lag in lags]
lag_r2 = [corr**2 if corr is not None else np.nan for corr in lag_corrs]

lagged_AABW = np.array(yband2_series.shift(-2))

# Print results
for lag, corr, r2 in zip(lags, lag_corrs, lag_r2):
    if lag < 0:
        relation = f"{band1} leads {band2} by {lag} steps"
    elif lag > 0:
        relation = f"{band2} leads {band1} by {abs(lag)} steps"
    else:
        relation = "No lead/lag (simultaneous)"
        
    if r2 > 0.2:
        print(f"{relation}: r = {corr:.3f}, R² = {r2:.3f}")

# --- Create subplots: top = data, bottom = lead–lag correlation ---
fig, (ax_data, ax_corr) = plt.subplots(2, 1, figsize=(8, 8), sharex=False)

# --- Top panel: time series data ---

ax_data.plot(x_sub, [value - yband1_sub[20] for value in yband1_sub], label=band1, color='blue')
ax_data.plot(x_sub, [value - yband2_sub[20] for value in yband2_sub], label=band2, color='orange')
ax_data.set_ylabel("Value")
ax_data.set_title("Time Series: SST")
ax_data.legend(loc="lower right")
ax_data.legend()
ax_data.grid(True)

# --- Bottom panel: lead–lag correlation ---
ax_corr.plot(lags, lag_corrs, marker='o', label='Correlation (r)')
#ax_corr.plot(lags, lag_r2, marker='s', linestyle='--', label='R²')

# Highlight the max correlation
max_idx = np.nanargmax(np.abs(lag_corrs))
ax_corr.axvline(lags[max_idx], color='red', linestyle=':', 
                label=f"Peak at lag {lags[max_idx]}")

# Shade regions for clarity
ax_corr.axvspan(min(lags), -1, color='lightblue', alpha=0.2, label=f'{band1} leads')
ax_corr.axvspan(1, max(lags), color='lightgreen', alpha=0.2, label=f'{band2} leads')

ax_corr.axhline(0, color='black', linewidth=0.8)
ax_corr.set_xlabel(f"Lag (positive = {band2} leads)")
ax_corr.set_ylabel("Correlation")
ax_corr.set_title("Lead–Lag Correlation")
ax_corr.legend()
ax_corr.grid(True)

plt.tight_layout()
plt.show()
#plt.savefig(expt + '_lead_lag_corr_salin_and_AABW.png')

