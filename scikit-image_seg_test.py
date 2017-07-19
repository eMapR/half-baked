# -*- coding: utf-8 -*-
"""
Created on Thu May 25 16:24:16 2017

@author: braatenj
"""


import numpy as np
from skimage import segmentation, color, io
from skimage.future import graph
from skimage.segmentation import mark_boundaries
from matplotlib import pyplot as plt



def weight_mean_color(graph, src, dst, n):
  """Callback to handle merging nodes by recomputing mean color.

  The method expects that the mean color of `dst` is already computed.

  Parameters
  ----------
  graph : RAG
      The graph under consideration.
  src, dst : int
      The vertices in `graph` to be merged.
  n : int
      A neighbor of `src` or `dst` or both.

  Returns
  -------
  data : dict
      A dictionary with the `"weight"` attribute set as the absolute
      difference of the mean color between node `dst` and `n`.
  """

  diff = graph.node[dst]['mean color'] - graph.node[n]['mean color']
  diff = np.linalg.norm(diff)
  return {'weight': diff}


def merge_mean_color(graph, src, dst):
  """Callback called before merging two nodes of a mean color distance graph.

  This method computes the mean color of `dst`.

  Parameters
  ----------
  graph : RAG
      The graph under consideration.
  src, dst : int
      The vertices in `graph` to be merged.
  """
  graph.node[dst]['total color'] += graph.node[src]['total color']
  graph.node[dst]['pixel count'] += graph.node[src]['pixel count']
  graph.node[dst]['mean color'] = (graph.node[dst]['total color'] /
                                   graph.node[dst]['pixel count'])


def make_display_image(img, bgPixels):
  """Function to apply 2 stdev stretch and compress to image to 8-bit for
  higher contrast viewing

  Parameters
  ----------
  img : np.array
      The image converted to np.array
  bgPixels : np.array
      np.array of indexes for background pixels
  """
  img[bgPixels] = np.nan
  imgMean = np.nanmean(img)
  imgStd = np.nanstd(img)
  imgMinOrig = np.nanmin(img)
  imgMinCalc = imgMean-(imgStd*2)
  if imgMinCalc < imgMinOrig:
    imgMin = imgMinOrig
  else:
    imgMin = imgMinCalc
  img[bgPixels] = imgMin
  imgMax = imgMean+(imgStd*2)
  img[img < imgMin] = imgMin
  img[img > imgMax] = imgMax
  img = np.round( ( (img-imgMin) / (imgMax-imgMin+0.0) ) *255).astype(np.uint8)
  return img

# activate the gdal plugin
io.use_plugin('gdal')

# define the image filename
fn = '/vol/v1/proj/cmonster/lidar/lt_lidar_comparison/aggregated/rasters/hja_lidar_biomass_reproject_nearest_masked_1x1.bsq'

# set background value
bgValue = -9999

# read in the image
img = io.imread(fn)

# get the index of the background pixels
bgPixels = np.where(img == bgValue)

# create the display image
displyImg = make_display_image(img, bgPixels)
plt.imshow(displyImg, cmap='gray')

# create the image to segment
img[bgPixels] = -100
img = np.round( ( (img-np.min(img)) / (np.max(img)-np.min(img)+0.0) ) * 600)

# use k-means to segment the image
seg1 = segmentation.slic(img, compactness=110, n_segments=4000)

# compare the original image to the initial segmentation
labels1Ave = color.label2rgb(seg1, img, kind='avg')
fig, (ax0, ax1) = plt.subplots(1, 2 ,figsize=(20,10),dpi=72)
ax0.imshow(img, cmap='gray')
ax1.imshow(labels1Ave)
ax0.axis('off')
ax1.axis('off')


#
g = graph.rag_mean_color(img, seg1)

labels2 = graph.merge_hierarchical(seg1, g, thresh=150, rag_copy=False,  #100
                                   in_place_merge=True,
                                   merge_func=merge_mean_color,
                                   weight_func=weight_mean_color)

labels2[bgPixels] = labels2[0,0]


plt.figure(1,figsize=(10,10),dpi=72)
plt.imshow(mark_boundaries(displyImg, labels2))





