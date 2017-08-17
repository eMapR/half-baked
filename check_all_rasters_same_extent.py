# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 12:17:52 2017

@author: braatenj
"""

from osgeo import gdal
from glob import glob
import os


def get_dims(fileName):
  src = gdal.Open(fileName)
  ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
  lrx = ulx + (src.RasterXSize * xres)
  lry = uly + (src.RasterYSize * yres)
  return [ulx,uly,lrx,lry]

##################################################################################################
# inputs
searchDir = '/vol/v1/general_files/user_files/justin/for_others/david/raster/mosaics_label/'
searchGlob = '*.bsq'
##################################################################################################

# normalize the searchDir
searchDir = os.path.normpath(searchDir)

# find the files that match the search term in the searchDir
files = glob(os.path.join(searchDir, searchGlob))

# get extent of 1st raster
standard = get_dims(files[0])

# loop through files to check if each extent is the same as the 1st raster
for fn in files[1:]:
  if not get_dims(fn) == standard:
    print("Not all rasters have the same extent")
    break