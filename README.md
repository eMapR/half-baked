# half-baked
A collection of one-off, messy, or hard-coded scripts that could be useful down the line.

When a new script is added, please alter this readme to include the script name under a broad category and provide a brief description of its purpose and any particularly useful tidbits that would help someone else run it or what could be changed to make it full-baked.

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


#### CONUS VRT Manipulation 

**spatial_subset_conus_vrt_by_geojson.py**

*Author:* Justin

*Date added:* 2017/07/19

*Description:* This script will subset/clip CONUS VRT files by a geojson file. The parameters are hard-coded. It is expecting certain VRT files to exist and is intended to work only on a pre-defined set of VRTs. For these reasons it is not a very flexible script, but could be made more generic by giving it glob-type search terms for the VRTs you want subset or provide full paths and output names.


