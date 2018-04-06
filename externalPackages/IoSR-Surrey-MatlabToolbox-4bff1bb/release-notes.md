# v2.8 - 17th June 2017

- Added kernelDensity function.
- boxPlot: Added ‘violin’ methods and properties for making violin plots.
- boxPlot: Added various violin-related options. Also modified the scatter offset to use kernel density rather than histogram.
- boxPlot: Added themeColor property, allowing the various colours for the theme to be changed via a single property.
- calcSnr: Prevent rank deficient warning when output is zeros.

# v2.7 - 28th March 2017

- Added true functional boxplot.
- Added alternative perceptual centroid function.
- Added caching of various dependent properties to iosr.bss.mixture class.
- Added identifiers to all warnings and errors.
- Updated loudness weighting calculations for A and C to better match IEC 61672. Added normalisation ISO 226 curve at 1kHz to bring it in to line with other weighting functions.
- Improved extrapolation of ISO 226 function. Added params output to return reference values.
- Updated SOFA API version in installer.
- Bug fixes and documentation improvements.

# v2.6 - 6th March 2017

Added a number of properties and methods to the iosr.bss.mixture class. Added whiskers property to iosr.statistics.functionalSpreadPlot. Other minor bug fixes and documentation improvements.

# v2.5 - 20th February 2017

Generalised iosr.bss.mixture class to support arbitrary spatial configurations, as well as none. Added decomposition properties to allow direct calculation of ideal masks, and application of masks. Currently, only STFT is supported; will add gammatone in due course. Added a number of of measures related to the mask and to signal overlap.

A few bug fixes.

# v2.4.1 - 18th January 2017

Added optional first argument to boxPlot and functionalSpreadPlot to specify axes in which to plot.

# v2.4 - 15th December 2016

Added functionalSpreadPlot for making functional box plots and related plots. Also added statsPlot class and moved some boxPlot methods to the class, as they are shared with functionalSpreadPlot.

# v2.3.2 - 13th October 2016

* Added check for OCTAVE in irStats, and a basic check for a unique peak.
* Corrected mistake in boxPlot documentation (pulled from TGabor).
* Moved important bits of documentation into chXcorr.m (the only separately-documented portion of the code).

# v2.3.1 - 25th July 2016

* Minor tweak to magnitude calculation in matchEQ.
* Redefined TIR in mixture class as target w.r.t. sum of interfering sources. Updated documentation.

# v2.3 - 19th July 2016

* Corrected code to restore current directory when installation is complete.
* Added function to generate BSS mixtures by combining sources in various ways.
* Added property to iosr.bss.mixture to return interferer filenames as char array. Also corrected bug setting properties when interferer comprises multiple sources.

# v2.2.3 - 12th July 2016

Corrected boxPlot bug whereby x-separator line would disappear when setting y-axis limits to inf.

# v2.2.2 - 10th July 2016

Corrected calls to other toolbox functions.

# v2.2.1 - 8th July 2016

Fixed erroneous default 'method' in boxPlot.

# v2.2 - 7th July 2016

Added install function that downloads dependencies and sets path. Fixed bug in boxPlot where some set functions check the old value rather than the new value.

# v2.1 - 22nd June 2016

Added match EQ function.

# v2 - 6th June 2016

Restructured toolbox into a Matlab package in order to appropriately restrict namespace. Consolidated some file/function names into a consistent format.
