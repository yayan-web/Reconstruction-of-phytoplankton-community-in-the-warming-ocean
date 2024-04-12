#!/usr/bin/env python
# coding: utf-8

# ## 大亚湾附近遥感数据探索

# In[30]:


import xarray as xr
import numpy as np
import pandas as pd
from plotnine import *
import geopandas as gpd
from pyproj import CRS
from scipy.interpolate import griddata
import os
import ctd
from pathlib import Path
from wodpy import wod
from os.path import relpath
import netCDF4
from datetime import datetime, timedelta
import warnings
from concurrent.futures import ThreadPoolExecutor
import psutil
from scipy.stats import pearsonr
import patchworklib as pw
warnings.filterwarnings('ignore')


# In[2]:


sst_files = os.listdir("sst/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr")
chl_files = os.listdir("chla/pub/socd1/mecb/coastwatch/viirs/science/L3/global/chlora/dineof")
lat_min, lat_max = 22, 22.6 # 从上到下....
lon_min, lon_max = 114.5, 115.2


# In[3]:


map_data = gpd.read_file("./2021-4/shi/CN-shi-A.shp")
daya_data = map_data[map_data["CityNameC"].isin(["惠州市", "深圳市"])]
daya_data_1 = daya_data.to_crs(CRS("EPSG:4326"))


# In[4]:


Stations = {
    'Station': ["A", "B", "C", "D", "E", "F"],
    'lat': [22.125, 22.125, 22.125, 22.375, 22.375,22.375],
    'lon': [114.625, 114.875, 115.125, 114.625, 114.875, 115.125]
}

Stations = pd.DataFrame(Stations)


# In[34]:


p1 = (ggplot() +
 theme_bw() +
 geom_map(data=daya_data_1, fill = "#8a7967") +
 geom_point(Stations, aes(x = "lon", y = "lat"), size = 7, color = "#ff4f81") +
 coord_fixed() +
 geom_label(Stations, aes(x = "lon", y = "lat-0.05", label = "Station"), size = 10) +
 xlim(lon_min, lon_max) +
 ylim(lat_min, lat_max) +
 labs(x = "Longitude", y = "Latitude") +
 theme(axis_text = element_text(size = 10, family = "Arial"),
         axis_title = element_text(size = 12, family = "Arial", face = "bold"),
         strip_text = element_text(size = 10, family = "Arial", face = "bold"),
         legend_text = element_text(size = 10, family = "Arial"),
         legend_title = element_text(size = 12, family = "Arial", face = "bold"),
         legend_position = "bottom"))
p1


# In[142]:


def convert_to_date(input_str):
    # 解析输入字符串，提取年份和第几天
    year_str, day_of_year_str = input_str[:4], input_str[4:]
    # 将字符串转换为整数
    year = int(year_str)
    day_of_year = int(day_of_year_str)
    # 计算日期
    base_date = datetime(year, 1, 1)  # 每年的第一天
    target_date = base_date + timedelta(days=day_of_year - 1)
    # 格式化日期为字符串
    result_str = target_date.strftime("%Y%m%d")
    return result_str


# In[ ]:


sst_files = [item for item in sst_files if "html" not in item]
i = 0
columns = ['lat', 'lon', 'sst', 'date']
sst_data = pd.DataFrame(columns=columns)
for dir in sst_files:
    dir = "sst/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/" + dir
    files  = os.listdir(dir)
    files =  [item for item in files if "html" not in item]
    for sst_file in files:
        i = i + 1
        file_path = dir + "/" + sst_file
        dataset = xr.open_dataset(file_path, engine = "h5netcdf")
        subset = dataset.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        df = subset.to_dataframe().reset_index()
        data_part = df[["lat", "lon", "sst"]]
        data_part = data_part.dropna()
        date = sst_file.split(".")[1]
        data_part["date"] = date
        sst_data = pd.concat([sst_data, data_part], ignore_index= True)
        if i % 365 == 0:
            print(dir, sst_file, i)


# In[ ]:


sst_data_copy = sst_data.copy() # 备份数据


# In[ ]:


sst_data["date"] = sst_data["date"].str.replace("_preliminary", "")


# In[ ]:


sst_data = pd.merge(sst_data, Stations, on = ["lat", "lon"])


# In[ ]:


sst_data['Year'] = sst_data['date'].str[:4].astype(int)
sst_data['Month'] = sst_data['date'].str[4:6].astype(int)
sst_data['Day'] = sst_data['date'].str[6:8].astype(int)


# In[ ]:


sst_data.to_csv("sst_data_1981-2023.csv")


# In[6]:


sst_data = pd.read_csv("sst_data_1981-2023.csv", header = 0, index_col = 0)


# In[100]:


sst_data.head()


# In[111]:


data = sst_data.groupby(["Year", "Month"])["sst"].mean().reset_index()


# In[114]:


data["Index"] = range(1, len(data) + 1)
data.head()


# In[126]:


(ggplot(data[data["Year"] < 2024], 
        aes(x = "Year", y = "sst", color = "sst"))
 + theme_bw()
 + geom_smooth(method = "lm")
 + geom_vline(xintercept = 1994)
 + geom_point(size = 2, shape = ".")
 + facet_wrap("~Month", scales = "free_y")
 + theme(axis_text = element_text(size = 6),
        strip_text = element_text(size = 8),
        figure_size = (8, 5)))


# In[7]:


sst_data_2 = sst_data[(sst_data["Month"]==4)&(sst_data["Day"].isin([5,6,7]))]
sst_data_2 = sst_data_2.drop("date", axis = 1)
# sst_data_2 = sst_data_2.groupby(["Year", "Month", "Station"]).mean().reset_index()


# ## 18年到23年几个点的平均温度

# In[39]:


sst_data_2.head()


# In[74]:


sst_data_3 = sst_data_2[sst_data_2["Year"].isin([2019, 2020, 2021, 2022, 2023])]
sst_data_3["Year"] = sst_data_3["Year"].astype("str")


# ## 这几天的平均温度

# In[48]:


sst_data_3["sst"].mean() # 19年到23年这些点的平均温度


# In[55]:


sst_data_3["Day"] = sst_data_3["Day"].replace(5, "April 5th")
sst_data_3["Day"] = sst_data_3["Day"].replace(6, "April 6th")
sst_data_3["Day"] = sst_data_3["Day"].replace(7, "April 7th")


# In[61]:


sst_data_3[(sst_data_3["Year"] == "2021") & (sst_data_3["Day"] == "April 6th")]["sst"].mean()


# In[63]:


p2 = (ggplot(sst_data_3, aes(x = "Year", y = "sst", fill = "Day", group = "Year"))
 + theme_bw()
 + geom_point(shape = ".", size = 2)
 + geom_boxplot(alpha = 0.4)
 + facet_wrap("~Day")
 + labs(x = "Year", y = "Temperature [°C]")
 + theme(axis_text_x = element_text(size = 10, family = "Arial", angle = 45, vjust = 1),
         axis_text_y = element_text(size = 10, family = "Arial"),
         axis_title = element_text(size = 12, family = "Arial", face = "bold"),
         strip_text = element_text(size = 10, family = "Arial", face = "bold"),
         legend_text = element_text(size = 10, family = "Arial"),
         legend_title = element_text(size = 12, family = "Arial", face = "bold"),
         figure_size = (8, 4),
         legend_position = "none")
)
p2


# In[65]:


p2.save("mean_sst.tiff",  width=8, height=4, units='in', dpi = 300)


# In[75]:


sst_data_3 = sst_data_3[sst_data_3["Year"]=="2021"]
sst_data_3 = sst_data_3.drop("Year", axis = 1)
sst_data_3 = sst_data_3.groupby(["Station"]).mean().reset_index()


# In[76]:


sst_data_3 = sst_data_3.drop(["Month", "Day"], axis = 1)


# In[77]:


sst_data_3


# In[78]:


grid_x, grid_y = np.mgrid[lon_min:lon_max:100j, lat_min:lat_max:100j]
grid_z = griddata((sst_data_3["lon"], sst_data_3["lat"]), sst_data_3["sst"], 
                  (grid_x, grid_y), method='cubic')
plot_df = pd.DataFrame({
    'x': grid_x.flatten(),
    'y': grid_y.flatten(),
    'z': grid_z.flatten()
})
plot_df.columns = ["lon", "lat", "sst"]


# In[79]:


plot_df = plot_df.dropna()
plot_df.tail()


# In[97]:


p5 = (ggplot() +
 theme_bw() +
 geom_map(data=daya_data_1, fill = "#8a7967") +
 geom_tile(plot_df, aes(x = "lon", y = "lat", fill = "sst")) +
 scale_fill_gradient(high = "#ff9933", low = "#0099cc", breaks = [24.6, 24.8, 25], name = "Temperature [°C]") +
 coord_fixed() +
 xlim(114.5, 115.15) +
 ylim(22.1, 22.6) +
 labs(x = "Longitude", y = "Latitude") + 
 theme(axis_text = element_text(size = 10, family = "Arial"),
         axis_title = element_text(size = 12, family = "Arial", face = "bold"),
         strip_text = element_text(size = 10, family = "Arial", face = "bold"),
         legend_text = element_text(size = 10, family = "Arial"),
         legend_title = element_text(size = 12, family = "Arial", face = "bold"),
      legend_position = "bottom")
)
p5


# ## 2021年这几天的平均温度的分布

# In[98]:


p5.save("sst_map.tiff", width=8, height=6, units='in', dpi=300)


# ## 获取各点位18-23年的区域叶绿素数据

# In[151]:


file_idx = 0
chl_files =  [item for item in chl_files if "html" not in item]
columns = ['lat', 'lon', "chlor_a", 'date']
chl_data = pd.DataFrame(columns=columns)
for dir in chl_files:
    dir = "chla/pub/socd1/mecb/coastwatch/viirs/science/L3/global/chlora/dineof/" + dir
    files  = os.listdir(dir)
    files =  [item for item in files if "html" not in item]
    for chl_file in files:
        file_idx += 1  
        file_path = os.path.join(dir, chl_file)  
        dataset = xr.open_dataset(file_path, engine="h5netcdf")
        subset = dataset.sel(lat=slice(lat_max, lat_min),
                             lon=slice(lon_min, lon_max))
        tmp_df = subset.to_dataframe().reset_index() 
        data_part = tmp_df[["lat", "lon", "chlor_a"]].dropna()
        data_part = data_part.drop_duplicates()
        date = chl_file.split("_")[0].replace("V", "")
        data_part["date"] = date
        chl_data = pd.concat([chl_data, data_part], ignore_index=True)
        if file_idx % 100 == 0:
            print(dir, chl_file)


# In[152]:


len(data_part)


# ## 获取各点位18-23年的Station平均叶绿素数据

# In[153]:


chl_data["date"] = chl_data["date"].apply(convert_to_date)
chl_data.head()


# In[154]:


chl_data.to_csv("chl_all_2018-2023.csv")


# In[ ]:





# In[ ]:


file_idx = 0
chl_files =  [item for item in chl_files if "html" not in item]
columns = ['lat', 'lon', "chlor_a", 'Station', 'date']
chl_data = pd.DataFrame(columns=columns)
for dir in chl_files:
    dir = "chla/pub/socd1/mecb/coastwatch/viirs/science/L3/global/chlora/dineof/" + dir
    files  = os.listdir(dir)
    files =  [item for item in files if "html" not in item]
    for chl_file in files:
        file_idx += 1  
        file_path = os.path.join(dir, chl_file)  
        dataset = xr.open_dataset(file_path, engine="h5netcdf")
        tmp_data = pd.DataFrame(columns=columns) 
        for station_idx in range(len(Stations)):  
            Station = Stations["Station"][station_idx]
            latitude = Stations["lat"][station_idx]
            longitude = Stations["lon"][station_idx]
            subset = dataset.sel(lat=slice(latitude + 0.5, latitude - 0.5),
                                 lon=slice(longitude - 0.5, longitude + 0.5))
            tmp_df = subset.to_dataframe().reset_index() 
            data_part = tmp_df[["lat", "lon", "chlor_a"]].dropna() 
            chl_mean = pd.DataFrame(data_part.mean()).transpose()
            chl_mean["Station"] = Station
            tmp_data = pd.concat([tmp_data, chl_mean], ignore_index=True)  
        date = chl_file.split("_")[0].replace("V", "")
        tmp_data["date"] = date  # Adding date to temp data
        chl_data = pd.concat([chl_data, tmp_data], ignore_index=True)
        if file_idx % 100 == 0:
            print(dir, chl_file)


# In[ ]:


chl_data.head(12)


# In[ ]:


chl_data_2 = chl_data.drop("lon", axis = 1)
chl_data_2 = chl_data_2.drop("lat", axis = 1)


# In[ ]:


chl_data_2.head()


# In[ ]:


chl_data_2["date"] = chl_data_2["date"].apply(convert_to_date)


# In[ ]:


df = pd.merge(chl_data_2, sst_data, on = ["date", "Station"])


# In[ ]:


df.head()


# In[ ]:


df.to_csv("chl_sst_data_2018-2023.csv")


# In[18]:


df = pd.read_csv("chl_sst_data_2018-2023.csv", header = 0, index_col = 0)


# ## 2020-2022，3年春季叶绿素与温度的关系

# In[19]:


df_2 = df[df.loc[:, "Year"].isin([2019, 2020, 2021,2022, 2023])]
df_2 = df_2[df_2.loc[:, "Month"].isin([3, 4, 5])]
df_2["Month"] = df_2["Month"].astype(str)
df_2["Sample"] = df_2["Year"].astype(str) + "_" + df_2["Station"].astype(str)


# In[20]:


(ggplot(df_2, 
        aes(x = "sst", y = "chlor_a", color = "Day"))
 + theme_bw()
 + geom_vline(xintercept = 25.5, linetype = "dashed", color = "red")
 + geom_point(size = 2, shape = ".")
 + xlim(23, 28)
 + geom_smooth(method='gls', 
             formula='y ~ I((x - 25.5) * (x > 25.5))', se=False, color='blue')
 + facet_wrap("~Sample", scales = "free_y", nrow = 5)
 + theme(axis_text = element_text(size = 6),
        strip_text = element_text(size = 6),
        legend_position = "bottom",
        figure_size = (8, 6)))


# ## 降水量

# In[21]:


pre_data = pd.read_csv("2站点2021.csv", header = 0)


# In[22]:


pre_data.columns = ["City", "Year", "Month", "Day", "Precipitation"]
pre_data = pre_data.groupby(["Year", "Month", "Day"])["Precipitation"].mean().reset_index()
pre_data = pre_data[pre_data["Month"].isin([2, 3, 4, 5])]


# In[23]:


(ggplot(pre_data, 
        aes(x = "Day", y = "Precipitation"))
 + theme_bw()
 + geom_point(size = 2, shape = ".")
 + facet_wrap("Month", nrow = 2)
 + theme(axis_text = element_text(size = 6),
        strip_text = element_text(size = 8),
        legend_position = "bottom",
        figure_size = (8, 5)))


# In[137]:


df_3 = df_2.copy()
df_3 = df_3[df_3["Year"] == 2021]


# 这部分按照降水量的先舍弃吧
# ```python
# condition1 = (df_3["Month"] == "3") & (df_3["Day"].isin(range(1, 20)))
# condition2 = (df_3["Month"] == "4") & (df_3["Day"].isin(range(25, 31)))
# condition3 = (df_3["Month"] == "5") & (df_3["Day"].isin(range(1, 7)))
# condition4 = (df_3["Month"] == "5") & (df_3["Day"].isin(range(24, 32)))
# condition5 = (df_3["Month"] == "4") & (df_3["Day"].isin(range(1, 14)))
# df_3 = df_3[~(condition1 | condition2 | condition3 | condition4 | condition5)]
# ```

# ### 为了绘图美观，改变ABC和DEF的顺序

# In[138]:


replace_dict = {"A": "D", "B": "E", "C": "F", "D": "A", "E": "B", "F": "C"}
df_3["Station"] = df_3["Station"].replace(replace_dict)


# In[139]:


p3 = (ggplot(df_3, 
        aes(x = "sst", y = "chlor_a", color = "Day"))
 + theme_bw()
 + geom_vline(xintercept = 25.5, linetype = "dashed", color = "red")
 + geom_vline(xintercept = 28, linetype = "dashed", color = "red")
 + geom_point(aes(shape = "Month"), size = 1.5, alpha = 0.5)
 + geom_smooth(method = 'rlm', 
               formula = 'y ~ I((x - 25.6) * (x > 25.6)) + I((x - 28) * (x > 28))', 
               se = False,
               color = '#ff4c4c',
               fill = "#fb8a2e")
 + labs(x = "SST [°C]", y = "Chlorophyll a [mg·m$^{-3}$]")
 + facet_wrap("~Station", scales = "free_y")
 + theme(axis_text = element_text(size = 10, family = "Arial"),
         axis_title = element_text(size = 12, family = "Arial", face = "bold"),
         strip_text = element_text(size = 10, family = "Arial", face = "bold"),
         legend_text = element_text(size = 10, family = "Arial"),
         legend_title = element_text(size = 12, family = "Arial", face = "bold"),
         figure_size = (8, 6),
         legend_position = "bottom"))
p3


# In[99]:


p1.save("Station.tiff", width=8, height=6, units='in', dpi=300)


# In[136]:


p3.save("Chl_sst.tiff", units='in', dpi=300)


# ## 三个温区的平均叶绿素

# In[ ]:





# ## 可能需要对调ABCDEF以保持和站位的上下一致

# ## 获取每天，叶绿素和温度的相关性，然后和平均温度的相关性
# 不好用，舍弃掉

# In[26]:


df_4 = df_2.copy()


# In[27]:


df_4 = df_4.drop(["Station", "Sample", "lat", "lon"], axis = 1)
df_4 = df_4.dropna()
df_4.tail()


# In[28]:


grouped_df = df_4.groupby(['date', "Year", 'Month', "Day"])


# In[29]:


correlation_df = pd.DataFrame(columns=['date', "Year", 'Month', "Day",
                                       'correlation', 'p_value', 'sst', 'chlor_a'])
for group_name, group_data in grouped_df:
    # 计算相关性和 p-value
    correlation, p_value = pearsonr(group_data['sst'], group_data['chlor_a'])
    
    sst = group_data['sst'].mean()
    chl = group_data['chlor_a'].mean()
    # 将结果添加到数据框
    correlation_df = pd.concat([correlation_df, 
                                pd.DataFrame({'date': group_name[0],
                                              "Year": group_name[1],
                                              'Month': group_name[2],
                                              "Day": group_name[3],
                                              'correlation': [correlation], 
                                              'p_value': [p_value],
                                              'sst': [sst],
                                              'chlor_a': [chl]})])

# 重置索引
correlation_df.head()


# In[30]:


result_df_2 = pd.melt(correlation_df, 
                      id_vars = ['date', "Year", 'Month', "Day", "sst"])
result_df_2["Day"] = result_df_2["Day"].astype("int")
result_df_2["Month"] = result_df_2["Month"].astype("int")
result_df_2.loc[result_df_2["Month"] >= 4, "Day"] += 31
result_df_2.loc[result_df_2["Month"] >= 5, "Day"] += 30
result_df_2["variable"] = pd.Categorical(result_df_2["variable"],
                                         categories=['sst', 'chlor_a', 'correlation', 'p_value'],
                                         ordered=True)
result_df_2.tail()


# In[71]:


result_df_3 = result_df_2.copy()
result_df_3 = result_df_3[result_df_3["variable"] == "chlor_a"]


# In[72]:


result_df_3 = result_df_3.groupby(["Year", "Month", "Day", "sst"])["value"].mean().reset_index()


# In[86]:


(ggplot(result_df_3,
        aes(x = "sst",  y = "value", color = "Day"))
 + theme_bw()
 + geom_vline(xintercept = 25.5)
 + facet_wrap("~Year", scales = "free_y")
 + geom_point(shape = ".", size = 2)
 + geom_smooth(method = "lm",
               formula='y ~ I((x - 25.5) * (x > 25.5)) + I((x - 28) * (x > 28))')
 + theme(figure_size = (8, 5),
        legend_position = "left"))


# In[ ]:


memory_info = psutil.virtual_memory()
print(memory_info)


# 探索降雨量和温度的关系 (太复杂了，暂舍弃)

# In[ ]:


print(df_2.dtypes)
df_2.tail()


# In[ ]:


df_6 = df[df.loc[:, "Year"].isin([2021])]
df_6 = df_6[df_6.loc[:, "Month"].isin([3, 4, 5])]
df_6["Month"] = df_6["Month"].astype(str)
df_6["Sample"] = df_6["Year"].astype(str) + "_" + df_6["Station"].astype(str)


# In[ ]:


df_6["Month"] = df_6["Month"].astype("int")
df_6 = pd.merge(df_6, pre_data, on = ["Year", "Month", "Day"])


# In[ ]:


# df_6.loc[df_6["Month"] >= 3, "Day"] += 30
df_6.loc[df_6["Month"] >= 4, "Day"] += 31
df_6.loc[df_6["Month"] >= 5, "Day"] += 30
df_6.tail()


# In[ ]:


df_7 = df_6[["Year", "Month", "Day", "sst", "Station", "Precipitation"]]
df_7 = df_7.groupby(["Year", "Month", "Day"])[["Precipitation", "sst"]].mean().reset_index()
df_7["Precipitation"] = df_7["Precipitation"].shift(3)
df_7 = df_7[df_7["Precipitation"] > 0]
df_7.head()


# In[ ]:


(ggplot(df_7, 
        aes(x = "sst", y = "Precipitation"))
 + theme_bw()
 # + geom_vline(xintercept = 21, linetype = "dashed", color = "red")
 + geom_point(size = 2, shape = ".")
 + theme(axis_text = element_text(size = 6),
        strip_text = element_text(size = 8),
        legend_position = "bottom",
        figure_size = (8, 5)))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




