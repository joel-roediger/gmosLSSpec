from pyraf import iraf
from iraf import gemini, gmos
from iraf import images, imutil

class reduceSpectrum():
	def __init__(self, rawIm, rawPath):
		self.rawIm = rawIm
		self.rawPath = rawPath
		self.options = {}
		# print "Working on image " + self.rawIm

	# function to view image within ds9
	def viewIm(self, image, frameNo=1):
		iraf.gdisplay(image, frame=frameNo)

	# function to reduce the raw spectroscopic image and view the result
	def reduce(self, image, pref="gs", bias="no", gscrrej="no", crspec="no", \
		qecorr="no", flat="no", fixpix="yes", vardq="no", biasPath="", \
		flatPath="", refimPath="", imPath="", bpmPath="gmos$data/chipgaps.dat", \
		logPath="", verbose="yes", frameNo=1):

		# TODO: stash relevant task options in attributes [using dict of dicts?]
		# self.options['gsreduce'] = {}
		# self.options['gsreduce'][...] = ...

		# next delete previous copies of the output (should they exist)
		iraf.imdel("g" + image)
		iraf.imdel(pref + image)

		# run gsreduce and view output
		iraf.gsreduce(image, outpref=pref, fl_bias=bias, fl_gscrrej=gscrrej, \
			fl_crspec=crspec, fl_qecorr=qecorr, fl_flat=flat, fl_fixpix=fixpix, \
			fl_vardq=vardq, bias=biasPath, flatim=flatPath, qe_refim=refimPath, \
			rawpath=imPath, bpm=bpmPath, logfile=logPath, verbose=verbose)
		
		image = pref + image
		self.viewIm(image, frameNo=frameNo)

		return

	# function to apply the wavelength calibration found from arc exposures
	def transform(self, image, lamRef, pref="t", vardq="no", logPath="", \
		verbose="yes", frameNo=2):

		# TODO: stash relevant task options in attributes

		# next delete previous copy of the output (should it exist)
		iraf.imdel(pref + image)

		# run gstransform and view output
		iraf.gstransform(image, outpref=pref, wavtraname=lamRef, \
			fl_vardq=vardq, logfile=logPath, verbose=verbose)

		image = pref + image
		self.viewIm(image, frameNo=frameNo)

		return

	# function to subtract the sky lines and continuum from an image
	def skysub(self, image, pref="s", vardq="no", sample="*", func="chebyshev", \
		order=1, low=2.5, high=2.5, inter="no", logPath="", verbose="yes", \
		frameNo=3):

		# next delete previous copy of the output (should it exit)
		iraf.imdel(pref + image)

		# run gsskysub and view output
		iraf.gsskysub(image, outpref=pref, fl_vardq=vardq, long_sample=sample, \
			function=func, order=order, low_rej=low, high_rej=high, \
			fl_inter=inter, logfile=logPath, verbose=verbose)

		image = pref + image
		self.viewIm(image, frameNo=frameNo)
		
		return


if __name__ == "__main__":
	# default initializations
	imPath = "data/"
	image = "N20121004S0120.fits"

	log = "reduceTarget_n7078_S0120.log"

	# create instance of the reduceSpectrum class
	n7078_spec = reduceSpectrum(image, imPath)

	# reduce the raw spectroscopic image and view output
	pref = "gs"
	biasIm = "data/bias_N20120928-1010.fits"
	flatIm = "data/flat_850_N20121004_v2.fits"
	n7078_spec.reduce(image, pref=pref, bias="yes", crspec="yes", flat="yes", \
		fixpix="no", vardq="yes", biasPath=biasIm, flatPath=flatIm, \
		imPath=imPath, logPath=log, verbose="no")
	image = pref + image

	# apply the wavelength solution to the reduced image and view output
	pref = "t"
	lamRef = "gsN20121004S0855"
	n7078_spec.transform(image, lamRef, pref=pref, vardq="yes", logPath=log, \
		verbose="no")
	image = pref + image

	# remove sky lines and continuum from the reduced image and view output
	pref = "s1"
	sample = "20:500,1800:2280"
	n7078_spec.skysub(image, pref=pref, vardq="yes", sample=sample, order=5, \
		high=1., inter="yes", logPath=log, verbose="no")
	image = pref + image

