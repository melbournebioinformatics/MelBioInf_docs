# Data Formats and Pre-processing

Mass Spectrometry data analysis is plagued by an overabundance of file formats.  The good news is that the Mass Spec community, including many instrument vendors have developed a standard file format for raw data, [mzML](http://www.psidev.info/mzml_1_0_0%20).  The bad news is that many of the old formats are still in widespread use, and most instruments don't produce it natively.  The reference implementation of the mzML standard is a software suite called [ProteoWizard](http://proteowizard.sourceforge.net/). ProteoWizard includes a very handy tool called `msconvert` that is capable of converting raw data from most instruments into `mzML` or into one of many other formats.  In addition to format conversion, `msconvert` can also perform a wide variety of noise filtering and peak-picking functions to prepare data for analysis.  A typical pre-processing involves;

1. Conversion from instrument .raw to mzML
2. Peak picking on both MS1 and MS2 data using vendor-native peak picking routines (built in to msconvert)
3. Denoising of MS2 data either by thresholding or by keeping only the largest peaks withing a moving window
4. Convert spectrum identifiers into a standardized format


To convert files from raw instrument native formats to `mzML` a windows PC is required.  If you need to do this, be sure to download [ProteoWizard](http://proteowizard.sourceforge.net/) with _vendor reader support_ .  This package comes with `MSConvertGUI` which allows conversion of raw files using a graphical interface.  Once files are in `mzML` or `mgf` format they can be converted to various other formats using the `msconvert3` tool in Galaxy.
