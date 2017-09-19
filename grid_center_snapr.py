# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 15:36:09 2017

@author: braatenj
"""

import sys
import gdal
import subprocess
import math



def nearestExt(ulx, uly, sizeX, sizeY):        
    print("Using nearest extent method")
    ulxAdj = math.ceil(ulx / 30.0) * 30
    if ulxAdj-ulx >= 15: 
        ulxAdj -= 30 + 15
    else:
        ulxAdj -= 15
    ulyAdj = round(uly / 30.0) * 30 + 15
    lrxAdj = ulxAdj + sizeX + 30
    lryAdj = ulyAdj + sizeY - 30 # yres is negative so need to + ulyAdj instead of -
                       
    return "{0} {1} {2} {3}".format(ulxAdj, lryAdj, lrxAdj, ulyAdj);
    

    
def maxExt(ulx, uly, lrx, lry):        
    print("Using max extent method")
    ulxAdj = math.floor(ulx / 30.0) * 30 - 15                                         
    ulyAdj = math.ceil(uly / 30.0) * 30 + 15
    lrxAdj = math.ceil(lrx / 30.0) * 30 + 15                                          
    lryAdj = math.floor(lry / 30) * 30 - 15                                               
    
    return "{0} {1} {2} {3}".format(ulxAdj, lryAdj, lrxAdj, ulyAdj);


    
def main(inFile, outFile, method, r):
    src = gdal.Open(inFile)
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    sizeX = src.RasterXSize * xres
    sizeY = src.RasterYSize * yres
    lrx = ulx + sizeX
    lry = uly + sizeY

    if method == "nearest":
        ext = nearestExt(ulx, uly, sizeX, sizeY)
    elif method == "max":
        ext = maxExt(ulx, uly, lrx, lry)
    
    print('...with ' + r + ' resampling')
    cmd = 'gdalwarp -r ' + r + ' -tr 30 30 -te ' + ext + ' ' + inFile + ' ' + outFile   
    subprocess.call(cmd)


    
if __name__ == "__main__":  
    args = sys.argv
    inFile = args[1]
    outFile = args[2]
    method = args[3]
    r = args[4]
    
    main(inFile, outFile, method, r)