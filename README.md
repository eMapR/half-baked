# half-baked
A collection of one-off, messy, or hard-coded scripts that could be useful down the line.

When a new script is added, please alter this readme to include the script name under a broad category and provide a brief description of its purpose and any particularly useful tidbits that would help someone else run it or what could be changed to make it full-baked.

[CONUS VRT Manipulation](#conVRTman)

[Image Spatial Segmentation](#spatialSeg)

[Videos](#videos)

[Raster manipulation](#rasterManipulation)

[LT Processing in IDL](#ltProcessingIDL)

**To add scripts to the repository:**

1. Make sure script file is in the /vol/v1/general_files/script_library/half-baked/ dir
2. In a terminal, `cd` to the half-baked script dir

`cd /vol/v1/general_files/script_library/half-baked/`

3. Add the file to the local GIT repository

`git add the_name_of_the_file.py`

4. Commit the file to the local repository

`git commit -m 'added file' the_name_of_the_file.py`

5. Push the local repository changes to the GitHub origin

`git push origin master`

**Commit changes to the scripts in the repository**

1. Make changes to a file
2. In a terminal, `cd` to the half-baked script dir

`cd /vol/v1/general_files/script_library/half-baked/`

3. Commit the file to the local repository

`git commit -m 'note about the change' the_name_of_the_file_that_was_changed.py`

4. Push the local repository changes to the GitHub origin

`git push origin master`


### <a id="conVRTman"></a>CONUS VRT Manipulation 

**spatial_subset_conus_vrt_by_geojson.py**

*Author:* Justin

*Date added:* 2017/07/19

*Description:* This script will subset/clip CONUS VRT files by a geojson file. The parameters are hard-coded. It is expecting certain VRT files to exist and is intended to work only on a pre-defined set of VRTs. For these reasons it is not a very flexible script, but could be made more generic by giving it glob-type search terms for the VRTs you want subset or provide full paths and output names.

### <a id="spatialSeg"></a>Image Spatial Segmentation 

**scikit-image_seg_test.py**

*Author:* Justin

*Date added:* 2017/07/19

*Description:* This script will spatially segment an image using the scikit-image library. **NEED TO FILL OUT MORE**

### <a id="videos"></a>Videos

**stem_nlcd_vids.py**

*Author:* Justin

*Date added:* 2017/07/19

*Description:* This script will create videos of STEM NLCD for even years 1990 to 2000. **NEED TO FILL OUT MORE**

### <a id="rasterManipulation"></a>Raster manipulation

**subset_raster_by_vector_extent.py**

*Author:* Justin

*Date added:* 2017/07/26

*Description:* This script will subset a raster by a vector extent using OGR and gdal_translate. It will also burn in a nodata value for pixels outside the bounds of the provided vector. All parameter arguments are hardcoded. Right now it only accepts geojson features - driver is hardcoded. **NEED TO FILL OUT MORE**

### <a id="ltProcessingIDL"></a>LT Processing in IDL

**execute_idl_runfiles_in_parallel.py**

*Author:* Justin

*Date added:* 2017/08/02

*Description:* This script will run IDL 'run files' or 'batchfiles' in parallel. The input is a list of files paths to the IDL .pro run files you want executed. This list is entered directly into the script.



