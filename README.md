python-code
===========

Everything so far relies on a specific folder structure, which is:

EXPERIMENT : series of recording sessions that makes up a study
SESSION : a single recording session

basedir : this is the parent folder for the EXPERIMENT, in it should be the following folders
	/data : this is where the b*.mat files live (the MATLAB fileconversion files)
	/experimentfiles
	 	*SESSION*.txt : penetration coordinates
		*SESSION*.jpg : scopephoto with penetrations marked on it (you make this during the experiment)
		ScopePhoto1.jpg : the scopephoto
	/fileconversion : this is where the hdf files are stored after conversion from *.mat. Each file is one recording site, filenames are like this *prefix*###.h5. Prefixes are RF for rf, VOC for vocalization, etc.
	
	cfs.txt : contains the characteristic frequency indices (not actual cfs, but which number stimulus was played. this file is created by RF.add_bf_man(), which just plots the RF as an (no. attenuations x no. frequencies) image)
	


fileconversion : converts .mat files to python-friendly .h5 files

