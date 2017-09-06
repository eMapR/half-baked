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

def scale_to_8bit(img, n_stdev):
  mean = np.mean(img)
  stdev = np.std(img)
  imin = mean-(stdev*n_stdev)
  imax = mean+(stdev*n_stdev)
  img[np.where(img < imin)] = imin
  img[np.where(img > imax)] = imax
  img = np.round(((img-imin)/(imax-imin+0.0))*255)     
  return img

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

inFile = '/vol/v1/proj/nccn/2017/raster/mora/ltee_nccn_mora_06010930_20170806/ltee_nccn_mora_06010930_20170806_ftv_nbr.bsq'
outfile = '/vol/v1/general_files/user_files/justin/temp/nbr_rgb_tiles/nbr_rgb.tif'
redBand   = 1
greenBand = 17
blueBand  = 33
source_srs = 'EPSG:5070' # what is the EPSG of the source imagery
target_srs = 'EPSG:3857' # what EPSG should the output images be
make_tiles = 1


#############################################################################################################
#############################################################################################################
#############################################################################################################

# get the image dimensions for image reading 
trgtDim = get_intersec([inFile])

# read in image bands and scale them to 8-bit range using 2 stdev - do this for R, G, and B bands
r = scale_to_8bit(get_band(inFile, trgtDim, redBand), 2)
g = scale_to_8bit(get_band(inFile, trgtDim, greenBand), 2)
b = scale_to_8bit(get_band(inFile, trgtDim, blueBand), 2)

# make a temp file name
tempOutFile = os.path.splitext(outfile)[0]+'_temp.tif'

# write out the temp 8-bit RGB composite 
write_bands(r, g, b, tempOutFile, inFile, trgtDim)

# write out the final reprojected RGB composite
cmd = 'gdalwarp -s_srs '+source_srs+' -t_srs '+target_srs+' -of GTiff -r near '+tempOutFile+' '+outfile
subprocess.call(cmd, shell=True)

# get rid of the temp file
os.remove(tempOutFile)

# possibly make tiles
if make_tiles:
  cmd = 'gdal2tiles.py -z 0-13 '+outfile+' '+os.path.dirname(outfile)+'/tiles'
  subprocess.call(cmd, shell=True)

