import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rasterio
import mplcursors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader
from matplotlib.colors import LightSource
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

# =========================
# 1. 读取DEM
# =========================

dem_file = "GRAY_HR_SR_W.tif"

with rasterio.open(dem_file) as src:
    dem = src.read(1)
    bounds = src.bounds
    extent = [
        bounds.left,
        bounds.right,
        bounds.bottom,
        bounds.top
    ]

# =========================
# 2. 东亚范围
# =========================

lon_min = 80
lon_max = 140
lat_min = 20
lat_max = 45

# =========================
# 3. 裁剪DEM
# =========================

height, width = dem.shape
lons = np.linspace(extent[0], extent[1], width)
lats = np.linspace(extent[3], extent[2], height)

lon_mask = (lons >= lon_min) & (lons <= lon_max)
lat_mask = (lats >= lat_min) & (lats <= lat_max)

dem_crop = dem[np.ix_(lat_mask, lon_mask)]
# 把海洋（高程 <= 0）设为 NaN，生成阴影时将变为透明
dem_crop = np.where(dem_crop < 0, np.nan, dem_crop)

crop_extent = [
    lons[lon_mask].min(),
    lons[lon_mask].max(),
    lats[lat_mask].min(),
    lats[lat_mask].max()
]
# =========================
# 4. 自定义灰度地形
# =========================

terrain_colors = [
    "#ececec",
    "#f0f0f0",
    "#d9d9d9",
    "#bdbdbd",
    "#969696",
    "#5F5959"
]

custom_cmap = LinearSegmentedColormap.from_list(
    "gray_terrain",
    terrain_colors,
    N=256
)

# =========================
# 5. hillshade阴影
# =========================

ls = LightSource(azdeg=315, altdeg=45)

# 用 0 临时填补 NaN，避免 hillshade 全白
dem_for_shade = np.nan_to_num(dem_crop, nan=0)

rgb = ls.shade(
    dem_for_shade,
    cmap=custom_cmap,
    blend_mode='soft',
    vert_exag=1.2
)

# 把原来是 NaN 的区域设为透明
rgb[..., -1] = np.where(np.isnan(dem_crop), 0, 1)
# =========================
# 6. 创建地图
# =========================

fig = plt.figure(figsize=(12, 10))
ax = plt.axes(projection=ccrs.PlateCarree())

ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

# =========================
# 7. 绘制地形
# =========================

ax.imshow(
    rgb,
    origin='upper',
    extent=crop_extent,
    transform=ccrs.PlateCarree(),
    zorder=1
)

# =========================
# 8. 添加国界线、河流
# =========================

ax.add_feature(
    cfeature.BORDERS,
    edgecolor='white',
    linewidth=1.5,
    zorder=5
)
# =========================
# 8.1 添加中国省界，并高亮云南/广东/江苏
# =========================
import cartopy.io.shapereader as shapereader

province_shp = shapereader.natural_earth(
    resolution='10m',
    category='cultural',
    name='admin_1_states_provinces'
)

province_reader = shapereader.Reader(province_shp)

# 需要高亮的省份及颜色
target_provinces = {
    "Yunnan": "yellow",
    "Guangdong": "green",
    "Jiangsu": "red"
}

china_count = 0

for record in province_reader.records():
    attrs = record.attributes

    # 判断是否属于中国
    is_china = (
        attrs.get('adm0_name') == 'China'
        or attrs.get('admin') == 'China'
        or attrs.get('geonunit') == 'China'
        or attrs.get('sovereignt') == 'China'
        or attrs.get('adm0_a3') == 'CHN'
        or attrs.get('gu_a3') == 'CHN'
        or attrs.get('iso_a2') == 'CN'
        or str(attrs.get('iso_3166_2', '')).startswith('CN-')
    )

    if not is_china:
        continue

    china_count += 1

    # 尝试获取省名（Natural Earth 不同版本字段可能不同）
    prov_name = (
        attrs.get('name')
        or attrs.get('name_en')
        or attrs.get('gn_name')
        or attrs.get('postal')
    )

    # 默认：其他省份不填充，仅画黑色边界
    facecolor = 'none'

    # 对目标省份填色
    if prov_name in target_provinces:
        facecolor = target_provinces[prov_name]

    ax.add_geometries(
        [record.geometry],
        crs=ccrs.PlateCarree(),
        facecolor=facecolor,
        edgecolor='black',
        linewidth=0.8,
        zorder=30
    )

print("中国省级行政区记录数：", china_count)
# =========================
# 8.2 添加长江和黄河
# =========================
ax.add_feature(
    cfeature.RIVERS,
    edgecolor='#3182bd',
    linewidth=1,
    zorder=6
)
df = pd.read_csv("species_points07.csv")


# =========================
# 10. 所有采集点绘制成同一个颜色
# =========================

point_color = "#d73027"   # 你可以改成任意颜色，例如 "black"、"#4575b4"

sc = ax.scatter(
    df["lon"],
    df["lat"],
    s=5,
    c=point_color,
    edgecolor='black',
    linewidth=0.4,
    transform=ccrs.PlateCarree(),
    label="Sampling points",
    zorder=10
)


# =========================
# 11. 添加经纬度边框标注
# =========================

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker

# 显示地图边框
ax.spines['geo'].set_visible(True)
ax.spines['geo'].set_linewidth(1.2)
ax.spines['geo'].set_edgecolor('black')

# 添加经纬度网格线
gl = ax.gridlines(
    crs=ccrs.PlateCarree(),
    draw_labels=True,
    linewidth=0.5,
    color='gray',
    alpha=0.6,
    linestyle='--',
    zorder=30
)

# 只在左侧和下侧显示经纬度标签
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = True
gl.left_labels = True

# 根据地图范围自动设置经纬度间隔
lon_span = lon_max - lon_min
lat_span = lat_max - lat_min

if lon_span <= 20:
    lon_step = 2
elif lon_span <= 60:
    lon_step = 5
else:
    lon_step = 10

if lat_span <= 20:
    lat_step = 2
elif lat_span <= 60:
    lat_step = 5
else:
    lat_step = 10

gl.xlocator = mticker.FixedLocator(
    np.arange(np.floor(lon_min), np.ceil(lon_max) + lon_step, lon_step)
)

gl.ylocator = mticker.FixedLocator(
    np.arange(np.floor(lat_min), np.ceil(lat_max) + lat_step, lat_step)
)

# 设置经纬度格式
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()

# 设置标签样式
gl.xlabel_style = {
    'size': 10,
    'color': 'black'
}

gl.ylabel_style = {
    'size': 10,
    'color': 'black'
}

# =========================
# 12. 海拔颜色图例
# =========================

norm = Normalize(
    vmin=0,
    vmax=9000
)

sm = ScalarMappable(
    cmap=custom_cmap,
    norm=norm
)

sm.set_array([])

cbar = plt.colorbar(
    sm,
    ax=ax,
    orientation='vertical',
    shrink=0.6,
    pad=0.03
)

cbar.set_label(
    'Elevation (m)',
    fontsize=12
)

cbar.outline.set_visible(False)

plt.legend(frameon=False, loc='lower left')

# =========================
# 13. 比例尺
# =========================

scalebar_length = 1000  # km

# 比例尺位置：这里仍然放在中国附近
# 如果你画全球图，可以根据需要改位置
lon0 = 120
lat0 = 25

scalebar_deg = scalebar_length / 111

ax.plot(
    [lon0, lon0 + scalebar_deg],
    [lat0, lat0],
    transform=ccrs.PlateCarree(),
    color='black',
    linewidth=3,
    zorder=20
)

ax.plot(
    [lon0, lon0],
    [lat0 - 0.3, lat0 + 0.3],
    color='black',
    transform=ccrs.PlateCarree(),
    linewidth=2,
    zorder=20
)

ax.plot(
    [lon0 + scalebar_deg, lon0 + scalebar_deg],
    [lat0 - 0.3, lat0 + 0.3],
    color='black',
    transform=ccrs.PlateCarree(),
    linewidth=2,
    zorder=20
)

ax.text(
    lon0 + scalebar_deg / 2,
    lat0 + 1,
    f'{scalebar_length} km',
    transform=ccrs.PlateCarree(),
    ha='center',
    fontsize=11,
    zorder=20
)

# =========================
# 14. 交互光标
# =========================

cursor = mplcursors.cursor(sc, hover=True)

@cursor.connect("add")
def on_add(sel):
    index = sel.index
    row = df.iloc[index]

    text = f"lon: {row['lon']:.3f}\nlat: {row['lat']:.3f}"

    if "species" in df.columns:
        text = f"{row['species']}\n" + text

    sel.annotation.set_text(text)

# =========================
# 15. 保存和显示
# 注意：保存要放在比例尺、散点、图例全部画完之后
# =========================

plt.savefig(
    "Word04.png",
    dpi=600,
    bbox_inches='tight',
    facecolor='white'
)
plt.savefig(
    "Word04.pdf",
    bbox_inches='tight',
    facecolor='white'
)
plt.show()