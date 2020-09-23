#!/usr/bin/env python3
# Ross Gray
# ross.gray12@imperial.ac.uk
#
#  imageanalysis4.py
#

""" Script for importing orthomosaic tiff files, conducting individual
tree crown delineation and analysing the ITC segments for pixel data. """


#################
##Loading modules
#################

import sys #For module writing
import subprocess #For bash and R operations
import os #For directory access
import glob #For file selection
import re #For manipulation of filenames

import cv2 #For image manipulation
from PIL import Image #For image analysis
import tifffile as tiff #For 16-bit tiff files

from osgeo import gdal, osr, ogr #For gdal spatial manipulation
# import fiona #For GIS manipulation
# import rasterio #For raster manipulaiton
# from rasterio.mask import mask #For raster manipulaiton
# from rasterio import Affine # or from affine import Affine
import geopandas as gpd #For shapefile manipulation
from shapely.geometry import mapping, point #For changing shapely geometries

import numpy as np #For image analysis and data manipulation
import csv #For saving the results table
import pandas as pd #For saving the results table and conversion to tables
import statistics #For calculating standard deviation

import time #For timing processing
from tqdm import tqdm #Produces progress bar

#################
###Georeferencing
#################

#~ ###Batch Georeferencing
#~ #Unable to do so mannaul georeferencing
#~ for f in glob.glob("../Data/Test/"+"*.tif"):

	#~ ##Filenames
	#~ filename = os.path.split(f)[1]
	#~ temp = "../Data/tmp/" + filename
	#~ out_file = "../Data/Geotest/" + filename

	#~ ##Translating original raster to new GCPs
	#~ #Command line
	#~ cmd1 = ["gdal_translate", "-of", "GTiff",
			#~ #GCPs   Source X, Source Y, Dest X, Dest Y
			#~ "-gcp", "566055", "-522336", "566058", "522336",
			#~ "-gcp", "566074", "-522269", "566077", "522269",
			#~ "-gcp", "566269", "-522273", "566271", "522273",
			#~ f, temp]
	#~ subprocess.call(cmd1)

	#~ ##Warping to create new raster
	#~ #Command line
	#~ cmd2 = ["gdalwarp", "-r", "near", "-tps", "-co", "COMPRESS=LZW", temp, out_file]
	#~ subprocess.call(cmd2)


#----------------------------------------------------------------------#
#----------------------------------------------------------------------#

###################################################
### TREE CROWN IDENTIFICAION AND PIXEL ANALYSIS ###
###################################################

def main(argv):

	###################################################
	## RASTERIZING LIDAR POLYGONS TO 16-BIT TEMPLATE ##
	###################################################

	## Timing the program
	starttime = time.time()

	#Update
	print("Starting image analysis.")

	##Interate through Orthomosaics
	def image_analysis(f):
	# for f in tqdm(glob.glob("../Data/LongTermStudy/Orthomosaics/Georeferenced/"+"*.tif")):
	## Translating RGB image into 16-bit template
		rgb_img = gdal.Open(f)

		#Getting filenames
		fullname = os.path.split(f)[1]
		filename = os.path.splitext(fullname)[0]

		#Output to new format
		template = gdal.Translate("../Data/LongTermStudy/Templates/NewGeoref/" + filename + "_ITC.tif", rgb_img,
								  format="GTiff", outputType= gdal.GDT_UInt16,
								  creationOptions=['COMPRESS=PACKBITS'])

		#Properly close the datasets to flush to disk
		rgb_img = None
		template = None

		##Rasterizing tree crown polygons onto template
		#Open RGB image, raster template and polygons to burn
		bgr_img = cv2.imread(f, cv2.IMREAD_COLOR)

		#Getting dimension of rgb image
		[height, width, dim] = bgr_img.shape

		#Burn tree crown polygons onto template
		if argv[2] == 'Manual':
			#For manual delineation
			rasterizeOptions = gdal.RasterizeOptions(format = "GTiff", width = width,
													 height = height, attribute = "Tree_ID",
													 outputType = gdal.GDT_UInt16)
			rasterpoly = gdal.Rasterize("../Data/LongTermStudy/Templates/NewGeoref/" + filename + "_ITC.tif",
								        "../Data/LongTermStudy/ITCSegmentation/Manual/ManualDelineation.shp",
								        options = rasterizeOptions)
		elif argv[2] == 'VTree':
			#For VTrees
			rasterizeOptions = gdal.RasterizeOptions(format = "GTiff", width = width,
													 height = height, attribute = "id",
													 outputType = gdal.GDT_UInt16)
			rasterpoly = gdal.Rasterize("../Data/LongTermStudy/Templates/NewGeoref/" + filename + "_ITC.tif",
								        "../Data/LongTermStudy/ITCSegmentation/VTreeDelineation.shp",
								        options = rasterizeOptions)
		elif argv[2] == '2DSFM':
			#For 3D Point Cloud delineation
			rasterizeOptions = gdal.RasterizeOptions(format = "GTiff", width = width,
													 height = height, attribute = "ID",
													 outputType = gdal.GDT_UInt16)
			rasterpoly = gdal.Rasterize("../Data/LongTermStudy/Templates/NewGeoref/" + filename + "_ITC.tif",
								        "../Data/LongTermStudy/ITCSegmentation/itc2DSFM/2DSFMDelineation.shp",
								        options = rasterizeOptions)
		elif argv[2] == 'LiDAR':
			#For LiDAR delineation
			rasterizeOptions = gdal.RasterizeOptions(format = "GTiff", width = width,
													 height = height, attribute = "ID",
													 outputType = gdal.GDT_UInt16)
			rasterpoly = gdal.Rasterize("../Data/LongTermStudy/Templates/NewGeoref/" + filename + "_ITC.tif",
								        "../Data/LongTermStudy/ITCSegmentation/itcLiDAR/LiDARDelineation.shp",
								        options = rasterizeOptions)

		#Properly close the datasets to flush to disk
		rasterpoly = None
		#~ rgb_img = None
		bgr_img = None

		print("Tree crown delineation successful.")

#----------------------------------------------------------------------#
#----------------------------------------------------------------------#

		#####################################
		## EXTRACTING FILENAME INFORMATION ##
		#####################################

		##Create empty arrays for the data
		Date = []
		Start_Time = []
		FlightNo = []
		Location = []
		Software = []
		Light = []
		Cloud = []
		Wind = []
		Height = []

		##Splitting path from filename
		fullname = os.path.split(f)[1]
		filename = os.path.splitext(fullname)[0]

		##When the image was taken
		#Date
		find_date = re.compile(r"_([0-9]+-[0-9]+-[0-9]+)_")
		date = find_date.search(filename).group(1)
		Date.append(date)

		#Start Time
		find_time = re.compile(r"_([0-9]+.[0-9]+)_")
		start_time = find_time.search(filename).group(1)
		start_time.replace(".", ":")
		Start_Time.append(start_time)

		##Where the image was taken
		find_where = re.compile(r"^([a-zA-Z]+\d+){1}_")
		where = find_where.search(filename).group(1)
		where = re.split('(\d+)',where)
		FlightNo.append(where[1])
		Location.append(where[0])

		##Software used to stitch
		find_software = re.compile(r"_([a-z]+)_")
		software = find_software.search(filename).group(1)
		if software[0] == "o":
			Software.append("OpenDroneMap")
		if software[0] == "p":
			Software.append("Pix4D")
		if software[0] == "d":
			Software.append("DroneDeploy")

		##Conditions it was taken in
		find_cond = re.compile(r"_([A-Z]+)")
		cond = find_cond.search(filename).group(1)
		#Explaning conditions
		if cond[0] == "D":
			Light.append("Dull")
		if cond[0] == "B":
			Light.append("Bright")
		if cond[0] == "S":
			Light.append("Sunny")
		if cond[1] == "N":
			Cloud.append("None")
		if cond[1] == "S":
			Cloud.append("Some")
		if cond[1] == "C":
			Cloud.append("Cloudy")
		if cond[1] == "O":
			Cloud.append("Overcast")
		if cond[2] == "N":
			Wind.append("None")
		if cond[2] == "L":
			Wind.append("Light")
		if cond[2] == "M":
			Wind.append("Medium")
		if cond[2] == "H":
			Wind.append("High")

		##() Height of flight in ft or m
		find_height = re.compile(r"_(\d+[a-zA-Z]+)_")
		height = find_height.search(filename).group(1)
		Height.append(height)

		print("Variables have been added.")

#----------------------------------------------------------------------#
#----------------------------------------------------------------------#
		####################################
		## PIXEL ANALYSIS OF ITC SEGMENTS ##
		####################################

		##Reading 16 bit ITC tiff file and bgr image
		treecrowns = tiff.imread("../Data/LongTermStudy/Templates/NewGeoref/" + filename + "_ITC.tif")
		img = cv2.imread(f)
		#Converting bgr to rgb
		b,g,r = cv2.split(img)
		img = cv2.merge([r,g,b])

		##Empty lists to populate
		Trees = []
		R_mean = []
		G_mean = []
		B_mean = []
		R_SD = []
		G_SD = []
		B_SD = []
		RCC = []
		GCC = []
		BCC = []
		ExG = []

		##Pixel Analysis
		print("Starting pixel analysis.")
		def pixel_analysis(i):
		# for i in tqdm(range(1, np.amax(treecrowns+1))):

			#Individual Tree Number
			Trees.append(i)

			#Finding tree crowns in array
			locations = np.where(treecrowns == i)

			#Extract based on tree crown cells
			values = img[locations]
			values = np.ma.masked_equal(values, 0)
			values_table = pd.DataFrame(values, columns = ["R", "G", "B"]) #Edit to speed up

			#Calculating mean for each colour channel
			Rmean = values_table["R"].mean()
			Gmean = values_table["G"].mean()
			Bmean = values_table["B"].mean()

			#Calculating standard deviation for each colour channel
			Rsd = values_table["R"].std()
			Gsd = values_table["G"].std()
			Bsd = values_table["B"].std()

			#Appending results
			R_mean.append(Rmean)
			G_mean.append(Gmean)
			B_mean.append(Bmean)
			R_SD.append(Rsd)
			G_SD.append(Gsd)
			B_SD.append(Bsd)

			#Calculating overall brightness
			rgb = (Rmean + Gmean + Bmean)

			#Calculating chromatic coordinates for each channel
			rcc = Rmean/rgb
			gcc = Gmean/rgb
			bcc = Bmean/rgb
			exg = (2*Gmean)/(Rmean+Bmean)

			#Appending chromatic coordinates to lists
			GCC.append(gcc)
			RCC.append(rcc)
			BCC.append(bcc)
			ExG.append(exg)
		# *map(pixel_analysis, tqdm(range(1, np.amax(treecrowns+1)))),
		[pixel_analysis(i) for i in tqdm(range(1, np.amax(treecrowns+1)))]

		print("Pixel Analysis Completed.")
#~ #----------------------------------------------------------------------#
#~ #----------------------------------------------------------------------#

		#########################
		## CREATING DATAFRAMES ##
		#########################

		##Converting results table to dataframe
		pixels_df = pd.DataFrame({"Tree_Crown_ID": Trees,
								  "R_Mean": R_mean, "G_Mean": G_mean, "B_Mean": B_mean,
								  "R_StDev": R_SD, "G_StDev": G_SD, "B_StDev": B_SD,
								  "RCC": RCC, "GCC": GCC, "BCC": BCC, "ExG": ExG})
		variables_df = pd.DataFrame({"Date" : Date, "Start_Time": Start_Time,
									 "Flight_Number": FlightNo, "Location": Location,
									 "Software": Software, "Light": Light, "Cloud": Cloud,
									 "Wind": Wind, "Height": Height})

		##Matching dataframe lenghts
		repeat_variables_df = pd.concat([variables_df]*len(pixels_df), ignore_index=True)

		##Combining dataframe
		combined_df = pd.concat([repeat_variables_df, pixels_df], axis = 1)

		##Rearranging dataframe
		#~ results_table = results_table[["Number", "Date", "Location", "Software", "Start_Time", "Height", "Light", "Cloud", "Wind", "Mean_Green", "GCC", "Mean_Red", "RCC", "Mean_Blue", "BCC", "ExG"]]

		## Saving dataframe to new csv or existing csv
		#Command line arguments to name the csv file
		if not os.path.isfile(argv[1]):
			combined_df.to_csv(argv[1], index=False)
		else:
			with open(argv[1], "a") as f:
				combined_df.to_csv(f, header = False, index=False)

	##Using map or list comprehension
	# *map(image_analysis, tqdm(glob.glob("../Data/LongTermStudy/Orthomosaics/Georeferenced/"+"*.tif"))),
	[image_analysis(f) for f in tqdm(glob.glob("../Data/LongTermStudy/Orthomosaics/Georeferenced/"+"*.tif"))]

	##Calculating time elapsed
	endtime = time.time()
	hours, rem = divmod(endtime-starttime, 3600)
	minutes, seconds = divmod(rem, 60)
	print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

#----------------------------------------------------------------------#
#----------------------------------------------------------------------#

if(__name__ == "__main__"):
	status = main(sys.argv)

#End
