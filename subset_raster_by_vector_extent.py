# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 10:51:30 2017

@author: braatenj
"""

from osgeo import ogr
import subprocess


# Get a Layer's Extent
inShape = '/vol/v2/stem/multi_lt_test/spatial/conus_tiles_northeast_epsg5070.geojson'
inRaster = '/vol/v1/general_files/datasets/spatial_data/nlcd/nlcd_2001_v2/nlcd_2001_landcover_3x3_equal.tif'
outRaster = '/vol/v2/stem/multi_lt_test/spatial/nlcd_2001_landcover_3x3_equal_ne.tif'

inDriver = ogr.GetDriverByName('GeoJSON')
inDataSource = inDriver.Open(inFile, 0)
inLayer = inDataSource.GetLayer()
extent = inLayer.GetExtent()

#projwin = '-projwin {} {} {} {} '.format(ulxy[0], ulxy[1], lrxy[0], lrxy[1])  
projwin = '-projwin {} {} {} {} '.format(extent[0], extent[3], extent[1], extent[2])  
cmd = 'gdal_translate -of GTiff -tr 30 30 ' + projwin + inRaster +' '+ outRaster
subprocess.call(cmd, shell=True)
