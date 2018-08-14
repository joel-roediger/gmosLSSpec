import os
import utils
from pyraf import iraf
from iraf import proto
from iraf import system
from iraf import fitsutil
from iraf import gemini, gmos
from iraf import images, imutil
from iraf import stsdas, analysis, gasp
from astropy.io import fits

class ReduceIfuObj():
	def __init__(self, imList, steps, rawPath="./", calPath="", mdf="", \
		biasIm="", flatList="", arcDir="", arcList="", logRoot=""):
		
		self.steps = steps

		self.rawPath = rawPath
		self.calPath = calPath

		self.bias = self.calPath + biasIm
		self.mdf = mdf
		self.imList = imList

		self.arcDir = arcDir
		self.arcList = arcList

		self.flats = utils.getImPaths(flatList)

		self.logRoot = logRoot

		return

	# function to mask cosmetics in exposures
	def adjustDQ(self, inIm, ext, txtMask):

		# allow user to dynamically declare the extension
		#ext = str(raw_input("\nWhich extension needs masking? (EXTNAME,EXTVER): "))
		#ext = "[" + ext + "]"

		print "\nADJUSTING MASK FOR " + inIm + ext.upper()
		
		# grab dimensions of mask from header of image
		# TODO: does expression for num need to be generalized to one-slit case?
		hdu = fits.open(inIm)
		num = 3 * (int(ext[-2]) - 1) + 2
		xDim = hdu[num].header["NAXIS1"]
		yDim = hdu[num].header["NAXIS2"]
		hdu.close()

		# display input image
		iraf.display(inIm + ext)

		# allow user to adjust mask iteratively
		while True:
			
			# ask user to input coordinate ranges to mask
			fixDQ = bool(input("\nDoes mask need improvement? (True/False): "))
			if fixDQ:
				os.system("rm -iv " + txtMask)
				ranges = raw_input("Coord ranges to be masked (space-separated): ")
				cmd = "echo '" + ranges + "' >> " + txtMask
				#print cmd
				os.system(cmd)
			else:
				break

			# transform text mask to pixel list format
			plMask = txtMask[:txtMask.find(".")] + ".pl"
			iraf.delete(plMask)
			iraf.text2mask(txtMask, plMask, xDim, yDim)

			# interpolate over bad pixels and mark them in DQ plane
			images = "tmp_" + inIm + "," + "x" + inIm + "," + "tmpdq_" + inIm
			iraf.imdelete(images)
			iraf.copy(inIm, "tmp_" + inIm)
			iraf.fixpix("tmp_" + inIm + ext, plMask, linterp="1,2,3,4")
			iraf.copy("tmp_" + inIm, "x" + inIm)
			iraf.imarith(plMask, "+", "x" + inIm + "[DQ," + ext[-2:], \
				"tmpdq_" + inIm)
			iraf.imcopy("tmpdq_" + inIm + "[0]", "x" + inIm + "[DQ," + \
				ext[-2:-1] + ",overwrite]")

			# display result
			iraf.display("x" + inIm + ext)

		return

	# function to clean cosmic ray hits from exposures
	def cleanCRs(self, inIm, pref="x", **kwargs):

		print "\nCLEANING COSMIC RAYS FROM " + inIm

		iraf.imdelete(pref + inIm)
		iraf.gemcrspec(inIm, pref + inIm, **kwargs)

		return

	# function to compute sensitivity function b/o observed flux standard
	def compSensFunc(self, inIm, **kwargs):

		print "\nCOMPUTING SENSITIVITY FUNCTION FROM " + inIm

		suffStart = inIm.find(".")
		date = inIm[suffStart - 13:suffStart - 5]

		hdu = fits.open(inIm)
		tgtName = hdu[0].header["OBJECT"]
		centWave = str(int(hdu[0].header["CENTWAVE"]))
		hdu.close()

		outFlux = tgtName + "_" + centWave + "_" + date + "_flux.txt"
		sFunc = tgtName + "_" + centWave + "_" + date + "_sFunc"

		iraf.gsstandard(inIm, sfile=outFlux, sfunction=sFunc, **kwargs)

		return

	# function to correct exposures for QE changes between CCD chips
	def correctQE(self, inIm, **kwargs):

		print "\nCORRECTING QE IN " + inIm

		# delete previous correction data and QE-corrected exposure
		iraf.imdelete("q" + inIm)
		iraf.imdelete("qecorr" + kwargs["refimages"])

		iraf.gqecorr(inIm, **kwargs)

		return

	# function to run the pipeline following the steps in K Labrie's tutorial
	def run(self):

		imPaths = utils.getImPaths(self.imList)
		
		for i, im in enumerate(imPaths):
			# attach MDF, subtract the bias+overscan, and trim the overscan
			if self.steps["reduceIm"]:
				logFile = self.logRoot + "_reduceIm.log"
				# TODO: remove hard-coding of slits and mdfdir parms
				self.reduceIm(im, rawpath=self.rawPath, slits="both", \
					fl_over="yes", fl_trim="yes", fl_gscrrej="no", \
					fl_wavtran="no", fl_skysub="no", fl_extract="no", \
					fl_fluxcal="no", fl_inter="no", fl_vardq="yes", \
					bias=self.bias, mdffile=self.mdf, mdfdir="./", \
					logfile=logFile, verbose="no")

			# model and subtract the scattered light
			if self.steps["subScatLgt"]:
				self.subtractScat("rg" + im, self.calPath + "blkMask_" + \
					self.flats[i], prefix = "b", fl_inter="yes", cross="yes")

			# clean the cosmic rays
			if self.steps["cleanCosRays"]:
				logFile = self.logRoot + "_cleanCosRays.log"
				self.cleanCRs("brg" + im, fl_vardq="yes", xorder=9, yorder=-1, \
					sigclip=4.5, sigfrac=0.5, objlim=1.0, niter=4, \
					key_ron="RDNOISE", key_gain="GAIN", logfile=logFile, \
					verbose="no")

			# correct for QE changes
			if self.steps["correctQE"]:
				logFile = self.logRoot + "_correctQE.log"

				# retrieve arc of matching central wavelength
				arc = utils.matchCentWave(im, self.rawPath, self.arcList, "arc/")
				refIm = "erg" + arc
				os.system("mv " + self.calPath + refIm + " ./")

				self.correctQE("xbrg" + im, refimages=refIm, fl_correct="yes", \
					fl_vardq="yes", logfile=logFile, verbose="no")

				# return arc from whence it came
				os.system("mv " + refIm + " " + self.calPath)

			# extract the spectra
			if self.steps["extractSpec"]:
				logFile = self.logRoot + "_extractSpec.log"

				# retrieve flat of matching central wavelength
				# TODO: set flat through call to matchCentWave?
				flat = self.flats[i]
				refIm = "eqbrg" + flat
				os.system("mv " + self.calPath + refIm + " ./")

				respIm = self.calPath + flat[:flat.find(".")] + "_resp.fits"

				self.extractSpec("qxbrg" + im, reference=refIm, \
					response=respIm, fl_inter="no", fl_vardq="yes", \
					recenter="no", trace="no", weights="none", logfile=logFile, \
					verbose="no")
				
				# return flat from whence it came
				os.system("mv " + refIm + " " + self.calPath)

			# TODO: place this in extractSpec
			# view extracted spectra in detector plane
			#for j in range(1, 3):
			#	print "\nDISPLAYING eqxbrg" + im + "[SCI," + str(j) + "]"
			#	iraf.display("eqxbrg" + im + "[SCI," + str(j) + "]")
			#	continue

			# adjust the mask to cover cosmetics
			if self.steps["maskSpec"]:
				hdu = fits.open("eqxbrg" + im)

				# create a separate mask for each science extension
				for j in range(1, hdu[-1].header["EXTVER"] + 1):
					txtMask = "mask_" + im[:im.find(".")] + "_" + str(j) + ".txt"
					ext = "[SCI," + str(j) + "]"
					self.adjustDQ("eqxbrg" + im, ext, txtMask)

			# rectify the extracted spectra
			if self.steps["rectifySpec"]:
				logFile = self.logRoot + "_rectifySpec.log"

				# retrieve arc of matching central wavelength
				arc = utils.matchCentWave(im, self.rawPath, self.arcList, "arc/")
				refIm = "erg" + arc
				os.system("mv " + self.calPath + refIm + " ./")

				# use first spectra to inform choice of final wavelength sampling
				print "\nRECTIFYING SPECTRA IN", "xeqxbrg" + im
				if i == 0:
					self.rectifySpec("xeqxbrg" + im, wavtraname=refIm, \
						fl_vardq="no", dw="INDEF", logfile=logFile, verbose="no")

				dw = input("Desired wavelength sampling: ")
				self.rectifySpec("xeqxbrg" + im, wavtraname=refIm, \
					fl_vardq="yes", dw=dw, logfile=logFile, verbose="no")

				# return arc from whence it came
				os.system("mv " + refIm + " " + self.calPath)

			# subtract the sky from the object fibers
			if self.steps["subtractSky"]:
				logFile = self.logRoot + "_subtractSky.log"
				self.subtractSky("txeqxbrg" + im, fl_inter="no", \
					logfile=logFile, verbose="yes")
			
			# ...


			continue

		return

	# function to extract spectra from the fibers and build a datacube
	def extractSpec(self, inIm, **kwargs):

		print "\nEXTRACTING SPECTRA FROM", inIm
		
		# destroy previous extracted spectra
		iraf.imdelete("e" + inIm)

		# extract the spectra and view the result
		iraf.gfextract(inIm, **kwargs)
		self.viewCube("e" + inIm, frame=1, z1=0., z2=0., extname="SCI", \
			version="*")

		return

	# function to rectify IFU spectra using established wavelength calibration
	def rectifySpec(self, inIm, **kwargs):

		# remove previous rectified spectra
		iraf.imdelete("t" + inIm)

		# transform the spectra and view the result
		iraf.gftransform(inIm, **kwargs)
		iraf.display("t" + inIm + "[SCI]")

		return

	# function to perform a tailored set of reductions on an IFU exposure
	def reduceIm(self, inIm, **kwargs):

		print "\nREDUCING", inIm

		# remove previously reduced exposures
		inPref = inIm[:inIm.find(".") - 14]
		if inPref == "":
			outPref = "g"
			images = outPref + inIm

			outPref = "r" + outPref
			images += "," + outPref + inIm
		else:
			outPref = ""
			images = ""

		if "fl_gscrrej" in kwargs.keys() and kwargs["fl_gscrrej"] == "yes":
			outPref = "x" + outPref
			images += "," + outPref + inIm
		else:
			pass

		if "fl_scatsub" in kwargs.keys() and kwargs["fl_scatsub"] == "yes":
			outPref = "b" + outPref
			images += "," + outPref + inIm
		else:
			pass

		if "fl_qecorr" in kwargs.keys() and kwargs["fl_qecorr"] == "yes":
			outPref = "q" + outPref
			images += "," + outPref + inIm
		else:
			pass

		if "fl_extract" in kwargs.keys() and kwargs["fl_extract"] == "no":
			pass
		else:
			outPref = "e" + outPref
			images += "," + outPref + inIm

		iraf.imdelete(images)

		# reduce the image
		iraf.gfreduce(inIm, **kwargs)

		return

	# function to subtract scattered light b/o bundle gaps identified from flats
	def subtractScat(self, inIm, mask, **kwargs):

		print "\nSUBTRACTING SCATTERED LIGHT FROM", inIm

		pref = kwargs["prefix"]

		# allow user to iterate over orders used to model scattered light
		while True:
			fixScat = bool(input("\nDoes subtraction of scattered light need improvement? (True/False): "))
			if fixScat:
				# request orders from user
				xOrders = raw_input("X orders (csv): ")
				yOrders = raw_input("Y orders (csv): ")

				# model and subtract the light, and view the result
				iraf.imdelete(pref + inIm)
				iraf.gfscatsub(inIm, mask, xorder=xOrders, yorder=yOrders, \
					**kwargs)

				utils.examIm(pref + inIm + "[SCI]", frame=1)
			else:
				break

		return

	# function to subtract the sky component from the object fibers
	def subtractSky(self, inIm, **kwargs):

		print "\nSUBTRACTING SKY FROM " + inIm

		# delete previous sky-subtracted spectra
		iraf.imdelete("s" + inIm)

		# subtract the sky and view the result (as image and datacube)
		iraf.gfskysub(inIm, **kwargs)
		iraf.display("s" + inIm + "[SCI]")
		self.viewCube("s" + inIm, extname="SCI", version="1")

		return

	# function to sum all the spectra from a datacube into one
	def sumFibers(self, inIm, **kwargs):

		print "\nSUMMING FIBERS FROM " + inIm

		# delete ...
		# TODO: add attribute to class to turn imdelete verifications on/off
		iraf.imdelete("a" + inIm)

		# sum the fibers and view the final spectrum
		iraf.gfapsum(inIm, **kwargs)
		iraf.splot("a" + inIm + "[SCI,1]")

		return

	# function to view an extracted datacube and its individual spectra
	def viewCube(self, cubeFile, **kwargs):

		iraf.gfdisplay(cubeFile, **kwargs)

		return


if __name__ == "__main__":

	# set root filename of logs
	logRoot = "v1895_sci"

	# declare which steps of pipeline to execute
	pipeSteps = {}
	pipeSteps["reduceIm"] = False
	pipeSteps["subScatLgt"] = False
	pipeSteps["cleanCosRays"] = False
	pipeSteps["correctQE"] = False
	pipeSteps["extractSpec"] = False
	pipeSteps["maskSpec"] = False
	pipeSteps["rectifySpec"] = False
	pipeSteps["subtractSky"] = False

	# setup object for VCC1895 and run pipeline 
	v1895_sci = ReduceIfuObj("targetFiles.txt", steps=pipeSteps, \
		rawPath="target/", calPath="calibrations/", mdf="gsifu_slits_mdf.fits", \
		biasIm="v1895_masterBias.fits", flatList="flatFiles.txt", \
		arcDir="arc/", arcList="arcFiles.txt", logRoot=logRoot)
	v1895_sci.run()
