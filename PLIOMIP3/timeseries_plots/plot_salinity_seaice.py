#python program
#
# this program will plot  SH salinity and sea ice on the same plot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

expt = 'xpsig'
corr_yearstart = 1500
corr_yearend = 1800

# Function to read data from a file, skipping the first line
def read_data(filename):
    x_vals = []
    y_vals = []
    with open(filename, 'r') as file:
        next(file)  # Skip the title line
        for line in file:
            try:
                x, y = map(float, line.strip().split(','))
                x_vals.append(x)
                y_vals.append(y)
            except ValueError:
                print(f"Skipping invalid line in {filename}: {line.strip()}")
    return x_vals, y_vals

# Function to compute running mean**
def running_mean(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')


##################################################################
plt.figure(figsize=(12, 9))
period = {'xpsid':'LP','xpsij':'LP490','xpsie':'EP400','xpsig':'EP'}
axis_label = {'salinity':'psu','seaice_area':'km2'}

#plot data from all the files
if expt == 'xpsij':
    startyear=1991
else:
    startyear=12

labels=[]
salfile = ('salinity_' + expt + '_60.0S-90S_' + str(startyear) + '_2999.tex')
cicefile = ('seaice_area_' + expt + '_SH_' + str(startyear) + '_2999.tex')
   
xsal, ysal = read_data(salfile)
xcice, ycice= read_data(cicefile)

# Plot the data
plt.subplot(211)
ax1=plt.gca()
ax2=ax1.twinx()

line1=ax1.plot(xsal, ysal, label='Salinity', color='blue' )
line2=ax2.plot(xcice, ycice, label='Sea ice area',color='orange' )
        
plt.title(expt + ': ' +  period.get(expt))

ax1.set_xlabel('year')
ax1.set_ylabel('psu',color='blue')
ax2.set_ylabel('km2',color='orange')

lines=line1 + line2
labels = [l.get_label() for l in lines]

ax1.legend(lines, labels, loc='lower left')

plt.grid(True)



#do a 30 year running mean
window=30
ysal_mean=running_mean(ysal,window)
xsal_mean=xsal[window-1:]
ycice_mean=running_mean(ycice,window)
xcice_mean=xcice[window-1:]

# Plot the data
plt.subplot(212)
ax1=plt.gca()
ax2=ax1.twinx()

line1=ax1.plot(xsal_mean, ysal_mean, label='Salinity', color='blue' )
line2=ax2.plot(xcice_mean, ycice_mean, label='Sea ice area',color='orange' )
        
plt.title(expt + '(30 year mean) : ' +  period.get(expt))

ax1.set_xlabel('year')
ax1.set_ylabel('psu',color='blue')
ax2.set_ylabel('km2',color='orange')

lines=line1 + line2
labels = [l.get_label() for l in lines]

ax1.legend(lines, labels, loc='lower left')

#plt.grid(True)
#plt.xlim(1500,2000)
ax1.grid(True)


# Save the plot to a file
plt.savefig(expt + '_salinity_and_seaice.png')
print("Plot saved to 'combined_plot.png'.")
plt.close()

##################################################################
# do a lag correlation on the raw data.  To avoid too high correlations
# only use xvalues between year start and yearend
print('LAG CORRELATIONS RAW DATA')
print('--------------------------------------')

xsal_arr = np.array(xsal)
mask=(xsal_arr >= corr_yearstart) & (xsal_arr <=corr_yearend)
ysal_sub = np.array(ysal)[mask]
ycice_sub=np.array(ycice)[mask]
print(np.shape(ysal_sub))
print(np.shape(ysal))

ysal_series = pd.Series(ysal_sub)
ycice_series = pd.Series(ycice_sub)

# Choose the maximum lag to test (in time steps)
max_lag = 30  # 30 years

# Compute correlations for lags from -max_lag to +max_lag
lags = range(-max_lag, max_lag + 1)
lag_corrs = [ysal_series.corr(ycice_series.shift(lag)) for lag in lags]

# Print results
for lag, corr in zip(lags, lag_corrs):
    if lag < 0:
        relation = f"Salinity leads Sea Ice by {lag} steps"
    elif lag > 0:
        relation = f"Sea Ice leads Salinity by {abs(lag)} steps"
    else:
        relation = "No lead/lag (simultaneous)"
        
    r_squared = corr**2 if corr is not None else np.nan
    if r_squared > 0.2:
        print(f"{relation}: r = {corr:.3f}, R² = {r_squared:.3f}")


##################################################################
# do a lag correlation on the 30 year averaged data
print('LAG CORRELATIONS 30 yr avg DATA ')
print('--------------------------------------')
xsal_arr = np.array(xsal_mean)
mask=(xsal_arr >= corr_yearstart) & (xsal_arr <=corr_yearend)
ysal_sub = np.array(ysal_mean)[mask]
ycice_sub=np.array(ycice_mean)[mask]
x_sub = xsal_arr[mask]
print(np.shape(ysal_sub))
print(np.shape(ysal))

ysal_series = pd.Series(ysal_sub)
ycice_series = pd.Series(ycice_sub)

# Choose the maximum lag to test (in time steps)
max_lag = 25  # 30 years


# Compute correlations for lags from -max_lag to +max_lag on 30 year mean data
lags = range(-max_lag, max_lag + 1)
lag_corrs = [ysal_series.corr(ycice_series.shift(lag)) for lag in lags]
lag_r2 = [corr**2 if corr is not None else np.nan for corr in lag_corrs]

lagged_sea_ice = np.array(ycice_series.shift(-2))

# Print results
for lag, corr, r2 in zip(lags, lag_corrs, lag_r2):
    if lag < 0:
        relation = f"Salinity leads Sea Ice by {lag} steps"
    elif lag > 0:
        relation = f"Sea Ice leads Salinity by {abs(lag)} steps"
    else:
        relation = "No lead/lag (simultaneous)"
        
    if r2 > 0.2:
        print(f"{relation}: r = {corr:.3f}, R² = {r2:.3f}")

# --- Create subplots: top = data, bottom = lead–lag correlation ---
fig, (ax_data, ax_corr) = plt.subplots(2, 1, figsize=(8, 8), sharex=False)

# --- Top panel: time series data ---
ax2=ax_data.twinx()
ax_data.plot(x_sub, ysal_sub, label='Salinity', color='blue')
ax2.plot(x_sub, ycice_sub, label='Sea Ice', color='orange')
#ax2.plot(x_sub, lagged_sea_ice, label='Sea Ice lagged', color='green')
ax_data.set_ylabel("Value")
ax_data.set_title("Time Series: Salinity and Sea Ice")
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
ax_corr.axvspan(min(lags), -1, color='lightblue', alpha=0.2, label='Salinity leads')
ax_corr.axvspan(1, max(lags), color='lightgreen', alpha=0.2, label='Seaice leads')

ax_corr.axhline(0, color='black', linewidth=0.8)
ax_corr.set_xlabel("Lag (positive = sea ice leads)")
ax_corr.set_ylabel("Correlation")
ax_corr.set_title("Lead–Lag Correlation")
ax_corr.legend()
ax_corr.grid(True)

plt.tight_layout()
plt.savefig(expt + '_lead_lag_corr_salinity_and_seaice.png')
