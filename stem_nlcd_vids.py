# -*- coding: utf-8 -*-
"""
Created on Fri May 12 13:39:27 2017

@author: braatenj
"""

from osgeo import gdal
import numpy as np
import math
import subprocess
from glob import glob
import os
import pandas as pd
import fnmatch
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


def get_subset_bounds(ulxy, urxy, targetWidth, targetHeight):
    ulxAdj = math.ceil(ulxy[0] / 30.0) * 30
    if ulxAdj-ulxy[0] >= 15: 
        ulxAdj -= 30 + 15
    else:
        ulxAdj -= 15
    ulyAdj = round(ulxy[1] / 30.0) * 30 + 15
    ratio = targetHeight/(targetWidth + 0.0)   
    geoWidth = round((urxy[0] - ulxy[0]) / 30.0) * 30
    geoHeight = geoWidth * ratio
    lrxAdj = ulxAdj + geoWidth
    lryAdj = ulyAdj - geoHeight
                       
    return ([ulxAdj, ulyAdj],[lrxAdj, lryAdj]) #"{0} {1} {2} {3}".format(ulxAdj, lryAdj, lrxAdj, ulyAdj);  



#############################################################################################################
#############################################################################################################
#############################################################################################################


nlcd_stem_dir = '/vol/v2/stem/conus_time_series/' # where to search for annual nldc files
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/dc_dulles' # where files will be written
colorMap = '/vol/v1/general_files/datasets/spatial_data/nlcd/nlcd_lc_color_map.csv' # csv colorMap file path

#robert vids
targetWidthSmall = 800
targetHeightSmall = 600


#DC - suburbs NW of Dulles International
ulxy = [1572411, 1935368]
urxy = [1588603, 1935305]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/dc_dulles'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#Orlando 
ulxy = [1391473, 738269]
urxy = [1457645, 739354]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/orlando'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#houston
ulxy = [499,803458]
urxy = [117982, 794527]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/houston'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#denver
ulxy = [-805607,1932341]
urxy = [-712113,1930027]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/denver'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#sacramento
ulxy = [-2183836,2055652]
urxy = [-2102402,2053351]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/sacramento'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#salt lake city
ulxy = [-1517900,2212730]
urxy = [-1264730,2218230]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/saltlakecity'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#phoenix
ulxy = [-1544566,1331888]
urxy = [-1394045,1330808]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/phoenix'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#portland
ulxy = [-2082049, 2816328]
urxy = [-2025392, 2815473]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/portland'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#reno
ulxy = [-2035314, 2108561]
urxy = [-1975306, 2108452]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/reno'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#las vegas
ulxy = [-1719719, 1635926]
urxy = [-1683180, 1635960]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/lasvegas'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#detroit
ulxy = [1020426, 2247734]
urxy = [1078629, 2247734]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/detroit'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#new orleans
ulxy = [572980, 787381]
urxy = [582252, 787408]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/neworleans'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)

#sioux falls
ulxy = [-69993, 2291375]
urxy = [-48032, 2291315]
outDir = '/vol/v1/general_files/user_files/justin/temp/stem/for_robert/siouxfalls'
ulxy, lrxy = get_subset_bounds(ulxy, urxy, targetWidthSmall, targetHeightSmall)


#############################################################################################################

# read in the color nlcd table
colorTable = pd.read_csv(colorMap)

# scan the annual nlcd images
maps = []  
for root, dirnames, filenames in os.walk(nlcd_stem_dir):
  for filename in fnmatch.filter(filenames, '*vote.tif'):
    if('_pct_' not in filename):    
      maps.append(os.path.join(root, filename))

# get a list of the years      
years = [os.path.basename(os.path.dirname(fn)) for fn in maps]

goods = [i for i, year in enumerate(years) if (int(year) >= 1990) and (int(year) % 2 == 0)]
maps = [maps[good] for good in goods]
years = [years[good] for good in goods]

# make the width even - needs even width and height for making mp4 vid
if (targetWidthSmall % 2) != 0:
  targetWidthSmall -= 1



# make sure the outDir exists
outDir = os.path.join(outDir)
if not os.path.exists(outDir):
    os.makedirs(outDir)
   
for fn, year in zip(maps, years):  
  outFile = os.path.join(outDir,'nlcd_'+year+'.tif')
  outFileTemp = outFile.replace('.tif', '_temp.tif')
  projwin = '-projwin {} {} {} {} '.format(ulxy[0], ulxy[1], lrxy[0], lrxy[1])  
  cmd = 'gdal_translate -of GTiff -r near -outsize '+str(targetWidthSmall)+' 0 '+projwin+fn+' '+outFileTemp
  subprocess.call(cmd, shell=True)
  color_map(outFileTemp, colorTable, outFile)
  os.remove(outFileTemp)


# add year label to each image
print('adding year label to images')
rgbFilesSmall = glob(outDir+'/*.tif')
for fn, year in zip(rgbFilesSmall,years):
  newFile = fn.replace('.tif', '.jpg') 
  cmd = 'convert -extent 800x600 '+fn+' -pointsize 72 -fill black -annotate +25+70 '+'"'+str(year)+'" '+newFile
  subprocess.call(cmd, shell=True)  

jpgFiles = glob(outDir+'/*.jpg')
for jpgFile in jpgFiles:
  for i in range(4):  
    newFile = jpgFile.replace('.jpg', '_'+str(i)+'.jpg')
    shutil.copyfile(jpgFile, newFile)
  os.remove(jpgFile)

  
# make webm and mp4 videos 
print('making videos')
search = outDir+'/*.jpg'
outWebm = outDir+'/vid.webm'
outMp4 = outWebm.replace('webm', 'mp4')

webmCmd = 'ffmpeg -framerate 25 -pattern_type glob -i "'+search+'" -c:v libvpx-vp9 -b:v 1M  '+outWebm #25
mp4Cmd = 'ffmpeg -framerate 25 -pattern_type glob -i "'+search+'" -vcodec libx264 -pix_fmt yuv420p '+outMp4 

subprocess.call(webmCmd, shell=True)  
subprocess.call(mp4Cmd, shell=True)   
  
