import os
import shutil
import utils
from pyraf import iraf
from iraf import system
from iraf import ttools
from iraf import gemini, gmos
from iraf import images, imutil

class ReduceIfuCals():

	# TODO: add attr to turn verify on/off for file deletions
	def __init__(self, target="", steps={}, rawMdf="", mdf="", mdfPath="./", \
		biasList="", biasPath="bias/", flatList="", flatPath="flat/", \
		logRoot=""):

		self.target = target
		self.logRoot = logRoot

		self.steps = steps

		self.rawMdf = rawMdf
		self.mdf = mdf
		self.mdfPath = mdfPath

		self.biasList = biasList
		self.biasPath = biasPath
		self.masterB = self.logRoot + "_masterBias.fits"

		self.flatList = flatList
		self.flatPath = flatPath

		return

	# function to construct a model for the wavelength solution b/o arc lines
	def calWave(self, arcList, **kwargs):

		iraf.gswavelength(arcList, **kwargs)

		# TODO: include call to gftransform to rectify arcs?
		
		return

	# function to check the viability of the raw MDF by extracting flat spectra
	def checkMdf(self, rawMdf, flatList, **kwargs):
		
		# copy raw MDF to local directory
		#mdf = rawMdf
		if not os.path.exists(self.mdf):
			print "\nFETCHING RAW MDF"
			rawPath = "/Users/roedigerj/anaconda/envs/astroconda/iraf_extern/gemini/gmos/data/"
			rawMdf = rawPath + rawMdf
			shutil.copy(rawMdf, self.mdfPath + self.mdf)
		
		# TODO: check that flatList is not empty

		# grab the name of a single flat to reduce
		im = utils.getImPaths(flatList)[0]

		print "\nCHECKING MDF WITH " + im

		# repeat the following steps until the MDF deemed correct
		while True:

			# do basic reduction of selected exposure
			_ = self.reduceIm(im, fl_inter="no", fl_gscrrej="no", \
				fl_extract="no", fl_wavtran="no", fl_skysub="no", \
				fl_fluxcal="no", mdffile=self.mdf, mdfdir=self.mdfPath, \
				verbose="no", **kwargs)

			# extract spectrum to check match between MDF and IFU flat
			logFile = ""
			if "logfile" in kwargs.keys():
				logFile = kwargs["logfile"]
			iraf.delete("database/aperg" + im[:im.find(".")] + "*", verify="yes")
			self.extractSpec("rg" + im, fl_inter="yes", logfile=logFile, \
				verbose="no")

			# TODO: remove any extracted files

			# call tcalc to fix MDF (if necessary)
			fixMdf = bool(input("\nDoes MDF need fixing? (True/False): "))
			if fixMdf:
				fiber = input("Fiber to fix: ")
				flag = input("Flag value: ")
				iraf.tcalc(self.mdf, "BEAM", "if NO == " + str(fiber) + \
					"then " + str(flag) + " else BEAM")
			else:
				break
		
		return

	# function to compute response function of GMOS IFU
	def computeResp(self, inIm, **kwargs):

		print "\nCOMPUTING RESPONSE FUNCTION FOR", inIm
		
		# find endpoint of image prefix
		prefEnd = inIm.find(".") - 14

		# form title of output image and compute it with gfresponse
		outIm = inIm[prefEnd:inIm.find(".")] + '_resp'
		iraf.imdelete(outIm, verify="yes")
		
		# compute and view response function and view it
		iraf.gfresponse(inIm, outIm, **kwargs)
		self.viewCube(outIm)

		return

	# function to ...
	def extractSpec(self, inIm, **kwargs):
		
		# remove previous spectra and extract new ones
		iraf.imdelete("e" + inIm, verify="yes")
		iraf.gfextract(inIm, **kwargs)

		return

	# function to produce a high-SNR master bias from a list of GMOS images
	def reduceBias(self, imList, outIm, **kwargs):

		print "\nREDUCING BIASES"
		
		# set up input redirection
		imList = "@" + imList

		# remove previous master
		iraf.imdelete(outIm, verify="yes")

		# reduce bias exposures and view master
		iraf.gbias(imList, outIm, **kwargs)
		util.viewIm(outIm)

		# remove prepared images
		iraf.imdelete("g" + imList)

		return

	# function to reduce flat exposures, incl. subtraction of scattered light
	def reduceFlats(self, flatList, inPref="", arcList="", arcPath="", \
		subscat="no", **kwargs):

		print "\nREDUCING FLATS"

		# retrieve filenames of individual exposures
		flatPaths = utils.getImPaths(flatList)

		# process flats one at a time
		for flat in flatPaths:

			# make additional preparations if QE correction requested
			if "fl_qecorr" in kwargs.keys() and kwargs["fl_qecorr"] == "yes":

				# move extracted arcs to working directory
				refIm = "erg" + utils.matchCentWave(inPref + flat, "./", \
					arcList, "arc/")
				os.system("mv " + arcPath + refIm + " ./")
				
				# add reference image to kwargs
				kwargs["qe_refim"] = refIm

			# model and remove the scattered light, if requested; 
			# otherwise, reduce the image
			if subscat == "yes":
				# TODO: display all science extensions

				# find file holding extracted flat spectra
				# TODO: determination of extFlat probably needs improvement
				dirList = os.listdir("./")
				descendList = [im for im in dirList if flat in im]
				extFlat = [im for im in descendList if "erg" in im][0]
				print extFlat

				# find gaps between fiber bundles on detector
				gapFile = "blkMask_" + flat[:flat.find(".")] + ".txt"
				iraf.delete(gapFile)
				iraf.gffindblocks(extFlat[1:], extFlat, gapFile)

	    		# model and remove scattered light, if necessary 
    			# (and repeat until satisfied)
				utils.examIm(extFlat[1:] + "[SCI]", frame=1)
				while True:
					fixScat = bool(input("Does subtraction of scattered light need improvement? (True/False): "))
					if fixScat:
						xOrders = raw_input("X orders (csv): ")
						yOrders = raw_input("Y orders (csv): ")

						iraf.imdelete("b" + extFlat[1:])
						iraf.gfscatsub(extFlat[1:], gapFile, fl_inter="yes", \
							prefix="b", xorder=xOrders, yorder=yOrders, \
							cross="yes")

						utils.examIm("b" + extFlat[1:] + "[SCI]", frame=1)
					else:
						break
			else:
				# reduce the image
				outPref = self.reduceIm(inPref + flat, **kwargs)

				# TODO: for this to work, need to know the proper prefix ...
				# (have reduceIm return the prefix?)
				self.viewCube(outPref + inPref + flat)
				pass

			# return extracted arcs to their rightful place
			if "fl_qecorr" in kwargs.keys() and kwargs["fl_qecorr"] == "yes":
				os.system("mv " + refIm + " " + arcPath)

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

		iraf.imdelete(images, verify="yes")

		# reduce the image
		iraf.gfreduce(inIm, **kwargs)

		return outPref

	# function to run the pipeline following the steps in K Labrie's tutorial
	def run(self):

		# view bias exposures
		if self.steps["viewBias"]:
			utils.viewIms(self.biasList, sat="yes")	

		# produce master bias
		if self.steps["reduceBias"]:
			logFile = self.logRoot + "_reduceBias.log"
			self.reduceBias(self.biasList, self.masterB, rawpath=self.biasPath, \
				fl_inter="no", fl_vardq="yes", logfile=logFile, verbose="no")

		# view flat exposures
		if self.steps["viewFlat"]:
			utils.viewIms(self.flatList, sat="yes")

		# check MDF (and modify, if necessary)
		if self.steps["checkMdf"]:
			#rawMdf = 
			logFile = self.logRoot + "_checkMdf.log"
			self.checkMdf(self.flatList, slits="both", rawpath=self.flatPath, \
				bias=self.masterB, logfile=logFile)

		# start reducing flats and create fiber trace for each central wavelength
		if self.steps["traceFibe"]:
			#mdf = "gsifu_slits_mdf.fits"
			logFile = self.logRoot + "_traceFibe.log"

			#rawFlat = "S20080405S0066"
			#v1895_cal.reduceIm(rawFlat, slits="both", rawpath="flat/", \
			#	fl_gscrrej="no", fl_wavtran="no", fl_skysub="no", \
			#	fl_fluxcal="no", fl_inter="no", fl_vardq="yes", \
			#	mdffile=self.mdf, mdfdir="./", bias=self.masterB, \
			#	logfile=logFile, verbose="no")

			#rawFlat = "S20080405S0070"
			#v1895_cal.reduceIm(rawFlat, slits="both", rawpath="flat/", \
			#	fl_gscrrej="no", fl_wavtran="no", fl_skysub="no", \
			#	fl_fluxcal="no", fl_inter="no", fl_vardq="yes", \
			#	mdffile=self.mdf, mdfdir="./", bias=self.masterB, \
			#	logfile=logFile, verbose="no")

			self.reduceFlats(flatList, slits="both", rawpath=self.flatPath, \
				fl_gscrrej="no", fl_wavtran="no", fl_skysub="no", \
				fl_fluxcal="no", fl_inter="no", fl_vardq="yes", \
				mdffile=self.mdf, mdfdir=self.mdfPath, bias=self.masterB, \
				logfile=logFile, verbose="no")

		# reduce the arcs
		if self.steps["reducsArcs"]:
			logFile = self.logRoot + "_reduceArcs.log"

			for arc in self.arcFiles:
				flat = utils.matchCentWave(arc, self.arcPath, self.flatList, \
					self.flatPath)

				_ = self.reduceIm(arc, slits="both", rawpath=self.arcPath, \
					fl_bias="no", fl_gscrrej="no", fl_wavtran="no", \
					fl_skysub="no", fl_fluxcal="no", fl_inter="no", \
					mdffile=self.mdf, mdfdir=self.mdfPath, \
					reference="erg" + flat, recenter="no", trace="no", \
					logfile=logFile, verbose="no")

		# determine wavelength solution for each central wavelength
		if self.steps["calWave"]:
			logFile = self.logRoot + "_calWave.log"

			self.calWave("erg@" + self.arcList, \
				coordlist="gmos$data/GCALcuar.dat", fl_inter="yes", \
				threshold=25., nlost=10, ntarget=15, logfile=logFile, \
				verbose="no")

		# reduce the lamp flat, including removal of scattered light
		if self.steps["subScatLgt"]:
			logFile = self.logRoot + "_subScatLgt.log"

			self.reduceFlats(self.flatList, subscat="yes", slits="both", \
				rawpath=self.flatPath, fl_gscrrej="no", fl_extract="no", \
				fl_wavtran="no", fl_skysub="no", fl_fluxcal="no", \
				fl_inter="no", fl_vardq="yes", mdffile=self.mdf, \
				mdfdir=self.mdfPath, bias=self.masterB, logfile=logFile, \
				verbose="no")
		
		# remove QE effects and view the results
		if self.steps["correctQE"]:
			logFile = self.logRoot + "_correctQE.log"

			self.reduceFlats(self.flatList, inPref="brg", arcList=self.arcList, \
				arcPath=self.arcPath, slits="both", fl_addmdf="no", \
				fl_bias="no", fl_over="no", fl_trim="no", fl_qecorr="yes", \
				fl_gscrrej="no", fl_extract="yes", fl_wavtran="no", \
				fl_skysub="no", fl_fluxcal="no", fl_inter="no", fl_vardq="yes", \
				logfile=logFile, verbose="no")

			#arc = "S20080405S0109"
			#imList = "flatFiles_473.txt"
			#imPaths = utils.getImPaths(imList)
			#for im in imPaths:
			#	v1895_cal.reduceIm("brg" + im, slits="both", fl_addmdf="no", \
			#		fl_bias="no", fl_over="no", fl_trim="no", fl_qecorr="yes", \
			#		fl_gscrrej="no", fl_wavtran="no", fl_skysub="no", \
			#		fl_extract="yes", fl_fluxcal="no", fl_inter="no", \
			#		fl_vardq="yes", qe_refim="erg" + arc, mdffile=self.mdf, \
			#		mdfdir=self.mdfPath, logfile=logFile, verbose="no")

			#v1895_cal.viewCube("eqbrg" + im)

		arc = "S20080405S0110"
		imList = "flatFiles_478.txt"
		imPaths = utils.getImPaths(imList)
		for im in imPaths:
			#v1895_cal.reduceIm("brg" + im, slits="both", fl_addmdf="no", \
			#	fl_bias="no", fl_over="no", fl_trim="no", fl_qecorr="yes", \
			#	fl_gscrrej="no", fl_wavtran="no", fl_skysub="no", \
			#	fl_extract="yes", fl_fluxcal="no", fl_inter="no", \
			#	fl_vardq="yes", qe_refim="erg" + arc, mdffile=mdf, mdfdir='./', \
			#	logfile=logFile, verbose="no")
			
			#v1895_cal.viewCube("eqbrg" + im)
			continue

		# calculate response function for each flat exposure
		imList = "flatFiles.txt"
		logFile = "v1895_respFunc.log"
		imPaths = utils.getImPaths(imList)
		for im in imPaths:
			v1895_cal.compResp("eqbrg" + im, skyimage="", fl_inter="no", \
				fl_fit="yes", function="spline3", order=45, sample="*", \
				logfile=logFile, verbose="no")
			continue

		return

	# function to view an extracted datacube and its individual spectra
	def viewCube(self, cubeFile, frameNo=1, z1=0., z2=0., ext="SCI", \
		version="*"):

		iraf.gfdisplay(cubeFile, frame=frameNo, z1=z1, z2=z2, extname=ext, \
			version=version)

		return


if __name__ == "__main__":

	# set root filename of logs
	target = "v1087"
	logRoot = target + "_cal"

	# setup names of raw and refined MDFs
	rawMdf = "gsifu_slits_mdf.fits"
	mdf = "mdf_cal.fits"
	mdfPath = "./"

	# setup ...
	biasList = "biasFiles.txt"

	flatList = "flatFiles.txt"

	# declare which steps of pipeline to execute
	pipeSteps = {}
	pipeSteps["viewBias"] = True
	pipeSteps["reduceBias"] = True
	pipeSteps["viewFlat"] = True
	pipeSteps["checkMdf"] = True
	pipeSteps["traceFibe"] = True

	# execute pipeline for calibrations
	calObj = ReduceIfuCals(target=target, steps=pipeSteps, rawMdf=rawMdf, \
		mdf=mdf, biasList=biasList, flatList=flatList, logRoot=logRoot)

