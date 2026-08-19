#! /bin/python3
import matplotlib as mpl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

plt.style.use("seaborn-v0_8-bright")                                            
plt.rcParams.update({                                                           
     "font.family": "serif",  # use serif/main font for text elements            
     "text.usetex": True,     # use inline math for ticks                        
     "pgf.rcfonts": True     # don't setup fonts from rc parameters              
         }) 

z_spec = 0.5
quantity = "T" 

dataFileName = f"{quantity}_Z{z_spec}_Ideal.csv"
cmaxFileName = f"{quantity}_Z{z_spec}_Cmax.csv"
hmaxFileName = f"{quantity}_Z{z_spec}_hmax.csv"
hminFileName = f"{quantity}_Z{z_spec}_hmin.csv"

outFileName = f"{quantity}_Z{z_spec}_normalized.csv"
res = 200

def read_points(path, quantity):
    return pd.read_csv(path, comment='#', header=None, names=[quantity, 'h', 'C'])


# Read main data
data = read_points(dataFileName, quantity)
n_total = len(data)
data_raw = data.copy()

# Read boundary data point sets
cmax_pts = read_points(cmaxFileName, quantity)
hmax_pts = read_points(hmaxFileName, quantity)
hmin_pts = read_points(hminFileName, quantity)

#######################################################################
# Normalize C using Cmax(h) from the boundary points
#######################################################################

h_env = cmax_pts['h'].values
Cmax_env = cmax_pts['C'].values
order = np.argsort(h_env)
h_env, Cmax_env = h_env[order], Cmax_env[order]
data['CNorm'] = data['C'] / np.interp(data['h'], h_env, Cmax_env)

#######################################################################
# Normalize enthalpy using hmin(CNorm) and hmax(CNorm)
#######################################################################

def to_CNorm_h(pts): 
    CNorm = pts['C'].values / np.interp(pts['h'].values, h_env, Cmax_env)
    h = pts['h'].values
    order = np.argsort(CNorm)
    return CNorm[order], h[order]

CNorm_hmin, hmin_line = to_CNorm_h(hmin_pts)
CNorm_hmax, hmax_line = to_CNorm_h(hmax_pts)

h_lo = np.interp(data['CNorm'], CNorm_hmin, hmin_line)
h_hi = np.interp(data['CNorm'], CNorm_hmax, hmax_line)
data['hNorm'] = ((data['h'] - h_lo) / (h_hi - h_lo)).clip(lower=0)

# Drop points outside the normalized bounds
data = data[(data['CNorm'] >= 0) & (data['CNorm'] <= 1) & (data['hNorm'] <= 1)]

#######################################################################
# Fill CNorm = 0 boundary
#######################################################################

edge_hNorm = np.linspace(0, 1, res)
if quantity == "T": 
    edge_point_val = 300.0
    colormap = "inferno"
    edge_points = pd.DataFrame({
        quantity: edge_point_val,
        'CNorm': 0.0,
        'hNorm': edge_hNorm,
    })
    data = pd.concat([data, edge_points], ignore_index=True)
else:
    colormap = "viridis"

#######################################################################
# Interpolate onto a cartesian grid in (hNorm, CNorm) space
#######################################################################

grid_coords = np.linspace(0, 1, res)
grid_h, grid_c = np.meshgrid(grid_coords, grid_coords)

points = data[['hNorm', 'CNorm']].values
values = data[quantity].values

# Linear interpolation
grid_q = griddata(points, values, (grid_h, grid_c), method='linear')

# Fill nan points with nearest-neighbor extrapolation
nan_mask = np.isnan(grid_q)
if np.any(nan_mask):
    grid_q_nearest = griddata(points, values, (grid_h, grid_c), method='nearest')
    grid_q[nan_mask] = grid_q_nearest[nan_mask]

#######################################################################
# Save the interpolated grid
#######################################################################

out = pd.DataFrame({
    'hNorm': grid_h.ravel(),
    'CNorm': grid_c.ravel(),
    quantity: grid_q.ravel(),
})
out.to_csv(outFileName, index=False)

#######################################################################
# Plot
#######################################################################
fig, axes = plt.subplots(1, 2, figsize=(12, 5))


cf = axes[0].contourf(grid_c, grid_h, grid_q, levels=40, cmap=mpl.colormaps[colormap])
axes[0].scatter(data['CNorm'], data['hNorm'], c=data[quantity], cmap=colormap,
                 s=6, edgecolors='k', linewidths=0.3)


axes[0].set_xlabel(r'$C_{norm}$')
axes[0].set_ylabel(r'$h_{norm}$')
fig.colorbar(cf, ax=axes[0], label=quantity)
axes[1].scatter(data_raw['C'], data_raw['h'], s=6, color='0.6', label='Rohdaten')
axes[1].scatter(Cmax_env, h_env, color='red', marker='^', label='$C_{max}(h)$')
hmin_sorted = hmin_pts.sort_values('C')
hmax_sorted = hmax_pts.sort_values('C')
axes[1].scatter(hmin_sorted['C'], hmin_sorted['h'], color='blue', marker='o', label=r'$h_{min}$')
axes[1].scatter(hmax_sorted['C'], hmax_sorted['h'], color='green', marker='D', label=r'$h_{max}$')
axes[1].set_xlabel('C')
axes[1].set_ylabel('h in J/kg')
axes[1].legend(fontsize=8)
title = f"Normierte und interpolierte Daten vs. Rohdaten bei " + r"$z = $" + f" {z_spec}"
fig.suptitle(title)
plt.show()

