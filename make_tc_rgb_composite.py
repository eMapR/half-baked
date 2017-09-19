# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 15:15:47 2017

@author: braatenj
"""


from osgeo import gdal
import numpy as np
import math
import subprocess
from glob import glob
import os
import shutil


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


def scale_to_8bit(img):
  mean = np.mean(img)
  stdev = np.std(img)
  n_stdev = 2
  imin = mean-(stdev*n_stdev)
  imax = mean+(stdev*n_stdev)
  if imin < 0:
    imin = 0
  img[np.where(img < imin)] = imin
  img[np.where(img > imax)] = imax
  img = np.round(((img-imin)/(imax-imin+0.0))*255)     
  return img


def scale_to_8bit_tc(img, tc):
  # standard TC stretch SR * 10000  
  n_stdev = 2  
  if tc == 'b':  
    imin = 3098-(1247*n_stdev)
    imax = 3098+(1247*n_stdev)
  if tc == 'g':
    imin = 1549-(799*n_stdev)
    imax = 1549+(799*n_stdev)
  if tc == 'w':  
    imin = -701-(772*n_stdev)
    imax = -701+(772*n_stdev)  
  
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
  
def color_map(inFile, colorTable, outFile):
  trgtDim = get_intersec([inFile])
  img = get_band(inFile, trgtDim, 1)
  r = np.copy(img)
  g = np.copy(img)
  b = np.copy(img)
  
  l = colorTable.shape[0]
  for i in range(colorTable.shape[0]):
    print('working on class: '+str(i+1)+'/'+str(l))
    these = np.where(img == colorTable.ix[i,'value'])
    r[these] = colorTable.ix[i,'r']
    g[these] = colorTable.ix[i,'g']
    b[these] = colorTable.ix[i,'b']
  
  write_bands(r, g, b, outFile, inFile, trgtDim)


def make_even(value):
  if (value % 2) != 0:
    value -= 1
  return value



#############################################################################################################
#############################################################################################################
#############################################################################################################


tcb = '/vol/v2/conus_tiles/vrts/lt_ee_conus_nbr_20170417_ftv_tcb.vrt'
tcg = '/vol/v2/conus_tiles/vrts/lt_ee_conus_nbr_20170417_ftv_tcg.vrt'
tcw = '/vol/v2/conus_tiles/vrts/lt_ee_conus_nbr_20170417_ftv_tcw.vrt'
band = 31
outFile = '/vol/v1/proj/field_trips/nccn_sept_2017/raster/images_for_tiling/mora_tc_z6-z16.tif'
clipFile = 0 #'/vol/v1/proj/nccn/2017/vector/LPa01_LEWI_MORA_NOCA_OLYM.shp' 
useFullImage = 0

source_srs = 'EPSG:5070' # what is the EPSG of the source imagery
target_srs = 'EPSG:3857' # what EPSG should the output images be - use EPSG
te_srs = 'EPSG:4326'# what EPSG are the following bounding coordinages in

ulxy = [-121.97296142578124,47.033630599264086] 
lrxy = [-121.35910034179688,46.68336307047754]

make_tiles = 0
zMin = 7
zMax = 14

#############################################################################################################


# make sure the outdir exists
if not os.path.isdir(os.path.dirname(outFile)):
  os.mkdir(os.path.dirname(outFile))


tempOutFiles = [os.path.splitext(outFile)[0] + '_'+tci + '_temp.tif' for fn, tci in zip([tcb, tcg, tcw],['tcb', 'tcg', 'tcw'])]
for inFile, tOutFile in zip([tcb, tcg, tcw], tempOutFiles):
  if useFullImage:
    cmd = 'gdalwarp -s_srs '+source_srs+' -t_srs '+target_srs+' -of GTiff -r near '+inFile+' '+tOutFile
  else:
    te = '-te {} {} {} {} '.format(ulxy[0], lrxy[1], lrxy[0], ulxy[1])  
    cmd = 'gdalwarp -s_srs '+source_srs+' -t_srs '+target_srs+' -of GTiff -r near -te_srs '+te_srs+' '+te+inFile+' '+tOutFile

  subprocess.call(cmd, shell=True)

trgtDim = get_intersec(tempOutFiles)

r = scale_to_8bit_tc(get_band(tempOutFiles[0], trgtDim, band), 'b')
g = scale_to_8bit_tc(get_band(tempOutFiles[1], trgtDim, band), 'g')
b = scale_to_8bit_tc(get_band(tempOutFiles[2], trgtDim, band), 'w')

outFile = os.path.splitext(outFile)[0]+'.tif' # make sure the file is a tif
write_bands(r, g, b, outFile, tempOutFiles[0], trgtDim)

# get rid of the temp file
deleteThese = glob(os.path.dirname(outFile)+'/*_temp.tif*')
for deleteThis in deleteThese:
  os.remove(deleteThis)
  
# possibly make tiles
if make_tiles:
  cmd = 'gdal2tiles.py -s '+target_srs+' -z '+str(zMin)+'-'+str(zMax)+' '+outFile+' '+os.path.dirname(outFile)+'/tiles'
  subprocess.call(cmd, shell=True)