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
	def __init__(self, biasList="", flatList=""):
		
		#self.logFile = logFile

		self.biasList = biasList
		self.flatList = flatList

		return

	# function to construct a model for the wavelength solution b/o arc lines
	def calWave(self, arcList, **kwargs):

		iraf.gswavelength(arcList, **kwargs)

		# TODO: include call to gftransform to rectify arcs?
		
		return

	# function to check the viability of the raw MDF by extracting flat spectra
	def checkMDF(self, rawMDF, flatList, **kwargs):
		
		# copy raw MDF to local directory
		mdf = rawMDF
		if not os.path.exists(mdf):
			print "\nFETCHING RAW MDF"
			mdfPath = "/Users/roedigerj/anaconda/envs/astroconda/iraf_extern/gemini/gmos/data/"
			rawMDF = mdfPath + rawMDF
			shutil.copy(rawMDF, "./")
		
		# TODO: check that flatList is not empty

		# grab the name of a single flat to reduce
		im = utils.getImPaths(flatList)[0]

		print "\nCHECKING MDF WITH " + im

		# repeat the following steps until the MDF deemed correct
		while True:

			# do basic reduction of selected exposure
			_ = self.reduceIm(im, fl_inter="no", fl_gscrrej="no", \
				fl_extract="no", fl_wavtran="no", fl_skysub="no", \
				fl_fluxcal="no", mdffile=mdf, mdfdir="./", verbose="no", \
				**kwargs)

			# extract spectrum to check match between MDF and IFU flat
			logFile = ""
			if "logfile" in kwargs.keys():
				logFile = kwargs["logfile"]
			iraf.delete("database/aperg" + im[:im.find(".")] + "*", verify="yes")
			self.extractSpec("rg" + im, fl_inter="yes", logfile=logFile, \
				verbose="no")

			# TODO: remove any extracted files

			# call tcalc to fix MDF (if necessary)
			fixMDF = bool(input("\nDoes MDF need fixing? (True/False): "))
			if fixMDF:
				fiber = input("Fiber to fix: ")
				flag = input("Flag value: ")
				iraf.tcalc(mdf, "BEAM", "if NO == " + str(fiber) + "then " + str(flag) + " else BEAM")
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

	# function to produce a high-SNR master from a list of GMOS bias images
	def reduceBias(self, imList, outIm, **kwargs):

		print "\nREDUCING BIASES"
		
		# append prefix to file list
		imList = "@" + imList

		# remove previous master bias
		iraf.imdelete(outIm)

		iraf.gbias(imList, outIm, **kwargs)

		# remove prepared images
		junk = "g" + imList
		iraf.imdelete(junk)

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

	def run(self):

		# view bias exposures
		#utils.viewIms(self.biasList, sat="yes")	

		# produce master bias
		logFile = "v1895_reduceBias.log"
		masterBias = "v1895_masterBias.fits"
		#self.reduceBias(self.biasList, masterBias, rawpath="bias/", \
		#	fl_inter="no", fl_vardq="yes", logfile=logFile, verbose="no")
		#util.viewIm(masterBias)


		# view flat exposures
		#utils.viewIms(self.flatList, sat="yes")

		# check MDF (and modify, if necessary)
		rawMDF = "gsifu_slits_mdf.fits"
		logFile = "v1895_checkMDF.log"
		#self.checkMDF(rawMDF, self.flatList, slits="both", rawpath="flat/", \
		#	bias=masterBias, logfile=logFile)
		
		
		# create fiber trace reference (for each central wavelength)
		# TODO: replace two calls to reduceIm with single call to reduceFlats
		mdf = "gsifu_slits_mdf.fits"
		logFile = "v1895_fiberTrace.log"

		rawFlat = "S20080405S0066"
		#v1895_cal.reduceIm(rawFlat, slits="both", rawpath="flat/", \
		#	fl_gscrrej="no", fl_wavtran="no", fl_skysub="no", fl_fluxcal="no", \
		#	fl_inter="no", fl_vardq="yes", mdffile=mdf, mdfdir="./", \
		#	bias=masterBias, logfile=logFile, verbose="no")
		
		rawFlat = "S20080405S0070"
		#v1895_cal.reduceIm(rawFlat, slits="both", rawpath="flat/", \
		#	fl_gscrrej="no", fl_wavtran="no", fl_skysub="no", fl_fluxcal="no", \
		#	fl_inter="no", fl_vardq="yes", mdffile=mdf, mdfdir="./", \
		#	bias=masterBias, logfile=logFile, verbose="no")
		
		
		# determine wavelength solution for each central wavelength
		logFile = "v1895_waveCal.log"

		arcs = ["S20080405S0109", "S20080405S0110"]
		extFlats = ["ergS20080405S0066", "ergS20080405S0070"]
		for i in range(2):
			#_ = v1895_cal.reduceIm(arcs[i], slits="both", rawpath="arc/", \
			#	fl_bias="no", fl_gscrrej="no", fl_wavtran="no", fl_skysub="no", \
			#	fl_fluxcal="no", fl_inter="no", mdffile=mdf, mdfdir="./", \
			#	reference=extFlats[i], recenter="no", trace="no", \
			#	logfile=logFile, verbose="no")
			continue
		
		arcList = "@arcFiles.txt"
		#v1895_cal.calWave("erg" + arcList, coordlist="gmos$data/GCALcuar.dat", \
		#	fl_inter="yes", threshold=25., nlost=10, ntarget=15, \
		#	logfile=logFile, verbose="no")
		
		
		# reduce the lamp flat (incl. removal of scattered light and QE correction)
		logFile = "v1895_subScat.log"
		#v1895_cal.reduceFlats(self.flatList, subscat="yes", rawpath="flat/", \
		#	slits="both", fl_gscrrej="no", fl_wavtran="no", fl_skysub="no", \
		#	fl_fluxcal="no", fl_inter="no", fl_vardq="yes", mdffile=mdf, \
		#	mdfdir="./", bias=masterBias, logfile=logFile, verbose="no")
		
		# TODO: replace loop with call to reduceFlats
		arc = "S20080405S0109"
		imList = "flatFiles_473.txt"
		logFile = "v1895_corrQE.log"
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

		# TODO: replace loop with call to reduceFlats
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

	#logFile = "v1895_reduceCals.log"
	v1895_cal = ReduceIfuCals(biasList="biasFiles.txt", \
		flatList="flatFiles.txt")

