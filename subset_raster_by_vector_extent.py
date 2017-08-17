# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 10:51:30 2017

@author: braatenj
"""

import os
from osgeo import ogr
import subprocess


######################################################################################################################
inShape = '/vol/v2/stem/multi_lt_test/spatial/conus_tiles_northeast_epsg5070.geojson'
inRaster = '/vol/v1/general_files/datasets/spatial_data/nlcd/nlcd_2001_v2/nlcd_2001_landcover_3x3_equal.tif'
outRaster = '/vol/v2/stem/multi_lt_test/spatial/nlcd_2001_landcover_3x3_equal_ne.tif'
clipRaster = 'true'
burnValue = 0
######################################################################################################################

# get the driver from the inShape file ext
ext = str.lower(os.path.splitext(inShape)[-1])
drivers = {'.shp'    :'ESRI Shapefile', 
           '.geojson': 'GeoJSON'}           
driver = ogr.GetDriverByName(drivers[ext])

# read in the inShape file and get the extent
inDataSource = driver.Open(inShape, 0)
extent = inDataSource.GetLayer().GetExtent()

# format the exent as -projwin arguments for gdal translate
projwin = '-projwin {} {} {} {} '.format(extent[0], extent[3], extent[1], extent[2])  

# make gdal_translate cmd and run it as subprocess
cmd = 'gdal_translate -of GTiff -tr 30 30 ' + projwin + inRaster +' '+ outRaster
subprocess.call(cmd, shell=True)

# if values outside the inShape should be chagned to some 
if clipRaster.lower() == 'true':
  cmd = 'gdal_rasterize -i -burn '+str(burnValue)+' '+inShape+' '+outRaster
  subprocess.call(cmd, shell=True)

