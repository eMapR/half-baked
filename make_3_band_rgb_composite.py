# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:02:26 2017

@author: braatenj
"""

from osgeo import gdal
import numpy as np
import math
import subprocess
import os
from glob import glob



def get_dims(fileName):
  src = gdal.Open(fileName)
  ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
  sizeX = src.RasterXSize
  sizeY = src.RasterYSize
  lrx = ulx + (sizeX * xres)
  lry = uly + (sizeY * yres)
  return [ulx,uly,lrx,lry,xres,-yres,sizeX,sizeY]

def make_geo_trans(fileName, trgtDim):
  src   = gdal.Open(fileName)
  ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
  return((trgtDim[0], xres, xskew, trgtDim[1], yskew, yres))

def get_intersec(files):
  ulxAll=[]
  ulyAll=[]
  lrxAll=[]
  lryAll=[]
  for fn in files:
    dim = get_dims(fn)
    ulxAll.append(dim[0])
    ulyAll.append(dim[1])
    lrxAll.append(dim[2])
    lryAll.append(dim[3])
  return([max(ulxAll),min(ulyAll),min(lrxAll),max(lryAll)])

def get_offsets(fileName, trgtDim):
  dim = get_dims(fileName)
  xoff = math.floor(abs(dim[0]-trgtDim[0])/dim[4])
  yoff = math.ceil(abs(dim[1]-trgtDim[1])/dim[4])
  xsize = abs(trgtDim[0]-trgtDim[2])/dim[4]
  ysize = abs(trgtDim[1]-trgtDim[3])/dim[4]
  return([int(i) for i in [xoff, yoff, xsize, ysize]])

def get_band(fileName, trgtDim, band):
  offsets = get_offsets(fileName, trgtDim)
  src = gdal.Open(fileName)
  band = src.GetRasterBand(band)
  array = band.ReadAsArray(
            offsets[0],
            offsets[1],
            offsets[2],
            offsets[3])
  return(array)

def write_img(outFile, refImg, trgtDim, nBands, dataType, of):
  convertDT = {
    'uint8': 1,
    'int8': 1,
    'uint16': 2,
    'int16': 3,
    'uint32': 4,
    'int32': 5,
    'float32': 6,
    'float64': 7,
    'complex64': 10,
    'complex128': 11
  }
  dataType = convertDT[dataType]
  geoTrans = make_geo_trans(refImg, trgtDim)
  proj = gdal.Open(refImg).GetProjection()
  dims = get_offsets(refImg, trgtDim)
  driver = gdal.GetDriverByName(of)
  driver.Register()
  outImg = driver.Create(outFile, dims[2], dims[3], nBands, dataType) # file, col, row, nBands, dataTypeCode
  outImg.SetGeoTransform(geoTrans)
  outImg.SetProjection(proj)
  return(outImg)

def scale_to_8bit(img, stretchMin, stretchMax):
  img[np.where(img < stretchMin)] = stretchMin
  img[np.where(img > stretchMax)] = stretchMax
  img = np.round(((img-stretchMin)/(stretchMax-stretchMin+0.0))*255)     
  return img

def stdev_stretch_limits(img, n_stdev):
  goods = img[np.where(img != 0)]  
  mean = np.mean(goods)
  stdev = np.std(goods)
  imin = mean-(stdev * n_stdev)
  imax = mean+(stdev * n_stdev)
  return [imin, imax]

def write_bands(r, g, b, outFile, ref, trgtDim):
  outImg = write_img(outFile, ref, trgtDim, 3, 'int8', 'GTIFF')
  outBand = outImg.GetRasterBand(1) 
  outBand.WriteArray(r)
  outBand = outImg.GetRasterBand(2) 
  outBand.WriteArray(g)
  outBand = outImg.GetRasterBand(3) 
  outBand.WriteArray(b)
  outImg = None




#############################################################################################################
####  INPUTS  ###############################################################################################
#############################################################################################################

inFile = '/vol/v2/conus_tiles/vrts/lt_ee_conus_nbr_20170417_ftv_nbr.vrt'
outFile = '/vol/v1/proj/field_trips/nccn_sept_2017/raster/images_for_tiling/mora_z6-z13_2.tif'
redBand   = 1
greenBand = 17
blueBand  = 33
source_srs = 'EPSG:5070' # what is the EPSG of the source imagery
target_srs = 'EPSG:3857' # what EPSG should the output images be web mercator: EPSG:3857 | wgs84 lat lon: EPSG:4326 | iGIS: EPSG:900913
stretchStdev = 2 # if this is greater than 0, then it will be stdev stretch, else stretchMin and stretchMax need to be correcty defined
stretchMin = 0
stretchMax = 0
useFullImage = 0 # if this is 1 the whole image will be used, if 0, then fill our the next three variables
te_srs = 'EPSG:4326'# what EPSG are the following bounding coordinages in - target extent srs
ulxy = [-122.4810791015625,47.20557536955536] # corvallis demo: [-123.51001739501952, 44.669629526920104]
lrxy = [-121.38519287109375,46.65414950711959]  # corvallis demo: [-123.15673828124999, 44.48744370505992]           
make_tiles = 0
zMin = 7
zMax = 16

#############################################################################################################
#############################################################################################################
#############################################################################################################

# make sure the outdir exists
if not os.path.isdir(os.path.dirname(outFile)):
  os.mkdir(os.path.dirname(outFile))

# band and spatial subset the image
tempoutFileTrans = os.path.splitext(outFile)[0]+'_temp_translate.tif'
tempoutFile = tempoutFileTrans # make a copy of this temp file name for later use and possible overwriting

if useFullImage:
  cmd = 'gdal_translate -of GTiff -b '+redBand+' -b '+greenBand+' -b '+blueBand+' '+inFile+' '+tempoutFileTrans  
else:
  te = '{} {} {} {}'.format(ulxy[0]-1, ulxy[1]+0.5, lrxy[0]+1, lrxy[1]-0.5)
  cmd = 'gdal_translate -of GTiff -b '+str(redBand)+' -b '+str(greenBand)+' -b '+str(blueBand)+' -projwin '+te+' -projwin_srs '+te_srs+' '+inFile+' '+tempoutFile

subprocess.call(cmd, shell=True)

# reproject the image if requested

if source_srs != target_srs:
  tempoutFileWarp = os.path.splitext(outFile)[0]+'_temp_warp.tif'  
  te = '-te {} {} {} {} '.format(ulxy[0], lrxy[1], lrxy[0], ulxy[1])  
  cmd = 'gdalwarp -s_srs '+source_srs+' -t_srs '+target_srs+' -of GTiff -r near -te_srs '+te_srs+' '+te+tempoutFileTrans+' '+tempoutFileWarp  
  #cmd = 'gdalwarp -s_srs '+source_srs+' -t_srs '+target_srs+' -of GTiff -r near -overwrite '+tempoutFileTrans+' '+tempoutFileWarp
  subprocess.call(cmd, shell=True)
  tempoutFile = tempoutFileWarp

# get the image dimensions for image reading 
trgtDim = get_intersec([tempoutFile])

# if a standard deviation stretch is requested, then find the limits
if stretchStdev > 0:
  stretchMin, stretchMax = stdev_stretch_limits(get_band(tempoutFile, trgtDim, redBand), stretchStdev)

# read in image bands and scale them to 8-bit range using 2 stdev - do this for R, G, and B bands
r = scale_to_8bit(get_band(tempoutFile, trgtDim, 1), stretchMin, stretchMax)
g = scale_to_8bit(get_band(tempoutFile, trgtDim, 2), stretchMin, stretchMax)
b = scale_to_8bit(get_band(tempoutFile, trgtDim, 3), stretchMin, stretchMax)

# write out the temp 8-bit RGB composite 
outFile = os.path.splitext(outFile)[0]+'.tif' # make sure the file is a tif
write_bands(r, g, b, outFile, tempoutFile, trgtDim)


# get rid of the temp file
deleteThese = glob(os.path.dirname(outFile)+'/*temp*.tif*')
for deleteThis in deleteThese:
  os.remove(deleteThis)

# possibly make tiles
if make_tiles:
  cmd = 'gdal2tiles.py -s '+target_srs+' -z '+str(zMin)+'-'+str(zMax)+' '+outFile+' '+os.path.dirname(outFile)+'/tiles'
  subprocess.call(cmd, shell=True)

