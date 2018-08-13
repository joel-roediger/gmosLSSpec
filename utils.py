import time
from pyraf import iraf
from iraf import system
from iraf import images, tv
from iraf import gemini, gmos
from astropy.io import fits

# function to ...
def examIm(im, **kwargs):
	
	# TODO: check whether ds9 open
	iraf.imexamine(im, **kwargs)

	return

# function to extract paths to individual images listed in a file
def getImPaths(imListPath):

	imList = open(imListPath, "r")
	scratch = imList.readlines()
	imPaths = [imPath.strip() for imPath in scratch]
	imList.close()

	return imPaths

# function to match images based on their central wavelengths
def matchCentWave(im1, im1Path, im2List, im2Path):

	# grab central wavelength of first image
	centWave1 = fits.open(im1Path + im1)[0].header["CENTWAVE"]

	# find best match amongst second image set via central wavelength
	im2Names = getImPaths(im2List)
	for im2 in im2Names:
		centWave2 = fits.open(im2Path + im2)[0].header["CENTWAVE"]
		if centWave2 == centWave1:
			break

	# return the best match
	return im2

# function to display a group of images (as a slowly-paced movie)
#def viewIms(imList, frameNo=1, bias="no", rawPath="./", z1=0., z2=0., sat="no"):
def viewIms(imList, **kwargs):

	# TODO: reformat loop to only use python functions (i.e. replace iraf.type)
	#images = iraf.type(imList, Stdout=1)
	imPaths = getImPaths(imList)
	for i, image in enumerate(imPaths):
		#image = rawPath + image.strip()
		print ("\nDISPLAYING " + image)

		# display each image for 3 s
		#iraf.gdisplay(image, frame=frameNo, fl_bias=bias, rawpath=rawPath, \
		#	z1=z1, z2=z2, fl_sat=sat)
		iraf.gdisplay(image, **kwargs)
		time.sleep(3)

	return
