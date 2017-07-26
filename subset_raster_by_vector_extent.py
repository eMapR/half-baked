# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 10:51:30 2017

@author: braatenj
"""

from osgeo import ogr
import subprocess


######################################################################################################################
inShape = '/vol/v2/stem/multi_lt_test/spatial/conus_tiles_northeast_epsg5070.geojson'
inRaster = '/vol/v1/general_files/datasets/spatial_data/nlcd/nlcd_2001_v2/nlcd_2001_landcover_3x3_equal.tif'
outRaster = '/vol/v2/stem/multi_lt_test/spatial/nlcd_2001_landcover_3x3_equal_ne.tif'
clipRaster = 'true'
######################################################################################################################



inDriver = ogr.GetDriverByName('GeoJSON') # this needs to be smarter so that a person can provide something besides geojson
inDataSource = inDriver.Open(inShape, 0)
extent = inDataSource.GetLayer().GetExtent()

projwin = '-projwin {} {} {} {} '.format(extent[0], extent[3], extent[1], extent[2])  
cmd = 'gdal_translate -of GTiff -tr 30 30 ' + projwin + inRaster +' '+ outRaster
subprocess.call(cmd, shell=True)

if clipRaster.lower() == 'true':
  cmd = 'gdal_rasterize -i -burn 0 '+inShape+' '+outRaster
  subprocess.call(cmd, shell=True)

