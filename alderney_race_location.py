import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import ConnectionPatch

# --- Tes coordonnées ---
lat_point = 49.712
lon_point = -2.010

# --- Configuration de la figure ---
fig = plt.figure(figsize=(10, 10)) 
main_ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# --- 1. La carte principale ---
main_ax.set_extent([-5, 10, 41, 52], crs=ccrs.PlateCarree())
main_ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='#f2f2f2')
main_ax.add_feature(cfeature.OCEAN.with_scale('10m'), facecolor='#b3d9ff')
main_ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
main_ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':', alpha=0.5)

# --- MODIFICATION ICI : Point plus gros sur la carte principale ---
# markersize=8 (au lieu de 3), + bordure blanche pour le faire ressortir
main_ax.plot(lon_point, lat_point, 'ro', markersize=8, markeredgecolor='white', markeredgewidth=2, transform=ccrs.PlateCarree())

# --- 2. Le Zoom en HAUT à DROITE ---
zoom_radius = 0.4 
zoom_extent = [lon_point - zoom_radius, lon_point + zoom_radius, 
               lat_point - zoom_radius, lat_point + zoom_radius]

zoom_ax = fig.add_axes([0.58, 0.58, 0.35, 0.35], projection=ccrs.PlateCarree())
zoom_ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())

zoom_ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='#f2f2f2')
zoom_ax.add_feature(cfeature.OCEAN.with_scale('10m'), facecolor='#b3d9ff')
zoom_ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
zoom_ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':', alpha=0.5)

# Point rouge sur le zoom
zoom_ax.plot(lon_point, lat_point, 'ro', markersize=8, markeredgecolor='white', markeredgewidth=2, transform=ccrs.PlateCarree())

# Cadre rouge autour du zoom
zoom_ax.spines['geo'].set_edgecolor('red')
zoom_ax.spines['geo'].set_linewidth(2)

# --- 3. Ligne de liaison plus visible ---
con = ConnectionPatch(xyA=(lon_point, lat_point), coordsA=main_ax.transData,
                      xyB=(lon_point - zoom_radius, lat_point), coordsB=zoom_ax.transData,
                      arrowstyle="-", 
                      linestyle="--", 
                      color="#333333",  # Gris plus foncé (presque noir)
                      linewidth=2.5,    # Épaisseur augmentée (2.5 au lieu de 1 par défaut)
                      alpha=1.0)        # Opacité maximale
fig.add_artist(con)

plt.title("Coordinates :    ", fontsize=22.5)

plt.savefig('C:/these_docs/mon_papier/covXtreme/figs/carte_zoom.png', dpi=300, bbox_inches='tight')
print("Image générée ! (Point gros et trait épais)")
plt.show()