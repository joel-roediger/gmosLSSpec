from pyraf import iraf
from iraf import gemini, gmos
from iraf import images, imutil
from iraf import noao, onedspec

class reduceSpectrum():
	def __init__(self, rawIm, rawPath):
		self.rawIm = rawIm
		self.rawPath = rawPath
		self.options = {}

	# function to view image within ds9
	def viewIm(self, image, frameNo=1):
		iraf.gdisplay(image, frame=frameNo)

		return

	# function to stash settings used to run any given task
	def stashSettings(self, task, **kwargs):
		self.options[task] = {}
		for kwarg in kwargs:
			self.options[task][kwarg] = kwargs[kwarg]

		return

	# function to reduce the raw spectroscopic image
	def reduce(self, image, pref="gs", bias="no", gscrrej="no", crspec="no", \
	 	qecorr="no", flat="no", fixpix="yes", vardq="no", biasPath="", \
		flatPath="", refimPath="", imPath="", bpmPath="gmos$data/chipgaps.dat", \
		logPath="", verbose="yes", frameNo=1):

		print "RUNNING GSREDUCE"

		# stash task settings
		self.stashSettings('gsreduce', outpref=pref, fl_bias=bias, \
			bias=biasPath, fl_flat=flat, flatim=flatPath, fl_qecorr=qecorr, \
			qe_refim=refimPath, fl_gscrrej=gscrrej, fl_crspec=crspec, \
			fl_fixpix=fixpix, fl_vardq=vardq, bpm=bpmPath, logfile=logPath)

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

		print "RUNNING GSTRANSFORM"

		# stash task settings
		self.stashSettings('gstransform', outpref=pref, wavtraname=lamRef, \
			fl_vardq=vardq, logfile=logPath)

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

		print "RUNNING GSSKYSUB"

		# stash task settings
		self.stashSettings('gsskysub', outpref=pref, fl_vardq=vardq, \
			long_sample=sample, function=func, order=order, low_rej=low, \
			high_rej=high, fl_inter=inter, logfile=logPath)

		# next delete previous copy of the output (should it exit)
		iraf.imdel(pref + image)

		# run gsskysub and view output
		iraf.gsskysub(image, outpref=pref, fl_vardq=vardq, long_sample=sample, \
			function=func, order=order, low_rej=low, high_rej=high, \
			fl_inter=inter, logfile=logPath, verbose=verbose)

		image = pref + image
		self.viewIm(image, frameNo=frameNo)
		
		return

	# function to extract the spectrum from the reduced, transformed, 
	# and sky-subtracted image
	def extract(self, image, pref="e", apRef="", width=1., inter="no", \
		center="yes", trace="yes", weights="none", vardq="no", logPath="", \
		verbose="yes", view=False):

		print "RUNNING GSEXTRACT"

		# stash task settings
		self.stashSettings('gsextract', outpref=pref, refimages=apRef, \
			apwidth=width, fl_inter=inter, recenter=center, trace=trace, \
			weights=weights, fl_vardq=vardq, logfile=logPath)

		# next delete previous copy of the output (should it exit)
		iraf.imdel(pref + image)

		# run gsextract and view output (if directed)
		iraf.gsextract(image, outpref=pref, refimages=apRef, apwidth=width, \
			fl_inter=inter, recenter=center, trace=trace, weights=weights, \
			fl_vardq=vardq, logfile=logPath, verbose=verbose)		

		# TODO: change 'view' to 'viewSpec' and move splot call to separate fcn
		if view:
			image = pref + image + "[sci]"
			iraf.splot(image)

		return


if __name__ == "__main__":
	# default initializations
	imPath = "data/"
	image = "N20121004S0120.fits"

	log = "reduceTarget_n7078_S0120.log"

	# create instance of the reduceSpectrum class
	n7078_850_q10 = reduceSpectrum(image, imPath)

	# reduce the raw spectroscopic image and view output
	pref = "gs"
	biasIm = "data/bias_N20120928-1010.fits"
	flatIm = "data/flat_850_N20121004_v2.fits"
	n7078_850_q10.reduce(image, pref=pref, bias="yes", crspec="yes", \
		flat="yes", fixpix="no", vardq="yes", biasPath=biasIm, flatPath=flatIm, \
		imPath=imPath, logPath=log, verbose="no")
	image = pref + image

	# apply the wavelength solution to the reduced image and view output
	pref = "t"
	lamRef = "gsN20121004S0855"
	n7078_850_q10.transform(image, lamRef, pref=pref, vardq="yes", logPath=log, \
		verbose="no")
	image = pref + image

	# remove sky lines and continuum from the reduced image and view output
	pref = "s1"
	sample = "20:500,1800:2280"
	n7078_850_q10.skysub(image, pref=pref, vardq="yes", sample=sample, order=5, \
		high=1., inter="no", logPath=log, verbose="no")
	image = pref + image

	# extract cluster spectrum
	# start by tracing position of a bright star as a function of wavelength
	weights = "variance"
	n7078_850_q10.extract(image, inter="yes", weights=weights, logPath=log)

	# next set up the first cluster aperture based on the trace of the star 
	pref = "e1"
	apRef = image
	n7078_850_q10.extract(image, pref=pref, apRef=apRef, width=10., \
		inter="yes", center="no", trace="no", weights=weights, vardq="yes", \
		logPath=log, verbose="no", view=True)

	# sanity check on stashed keywords
	print ''
	for key in n7078_850_q10.options.keys():
		print key, n7078_850_q10.options[key]

