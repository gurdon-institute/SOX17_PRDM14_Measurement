# Segments nuclei in C1 and calls positive or negative for each nucleus in other channels using inclusive empirically determined thresholds.
# Channel identities are parsed from image titles by taking substring matching "\w+-\d{3}", where the leading word is the gene name and the 
# trailing digits are the fluorophore wavelength.
#
# - by Richard Butler, Gurdon Institute Imaging Facility

import re, math

from ij import IJ, WindowManager, Prefs, ImagePlus
from ij.plugin import Duplicator
from ij.process import ImageStatistics, ImageProcessor
from ij.measure import Calibration, Measurements, ResultsTable
from ij.gui import Roi, ShapeRoi, Overlay

from java.awt import Color
from java.awt.geom import Ellipse2D

from org.jfree.chart import JFreeChart, ChartFactory, ChartFrame
from org.jfree.chart.plot import PlotOrientation
from org.jfree.data import Range
from org.jfree.data.xy import XYSeries, XYSeriesCollection


####### SETTINGS ########
MINA = 30.0		#minimum nucleus area (µm²)
MAXA = 300.0	#maximum nucleus area (µm²)

thresholds = { "GFP":35, "P14":55, "SOX17":45, "AP2y":20, "SOX2":45, "Venus":float("inf"), "TIR1":float("inf") }
#########################


def getMask(imp, chan):
	dup = Duplicator()
	mask = dup.run(imp, chan, chan, 1,1, 1,1)
	sigma = 0.2
	IJ.run(mask, "Gaussian Blur...", "sigma="+str(sigma)+" scaled")
	if chan==1:
		method = "Otsu"
	else:
		method = "MaxEntropy"
	IJ.setAutoThreshold(mask, method+" dark")
	Prefs.blackBackground = True
	IJ.run(mask, "Convert to Mask", "")
	IJ.run(mask, "Close-", "")
	IJ.run(mask, "Watershed", "")
	return mask

def plot2D(points, Ca, Cb):
	maxIntensity = 255.0
	dataset = XYSeriesCollection()

	seriesNN = XYSeries(channels[Ca+1]+" -ve "+channels[Cb+1]+" -ve")
	seriesPP = XYSeries(channels[Ca+1]+" +ve "+channels[Cb+1]+" +ve")
	seriesNP = XYSeries(channels[Ca+1]+" -ve "+channels[Cb+1]+" +ve")
	seriesPN = XYSeries(channels[Ca+1]+" +ve "+channels[Cb+1]+" -ve")
	for p in points:
		posA = channels[Ca+1] in thresholds and p[Ca]>thresholds[ channels[Ca+1] ]
		posB = channels[Cb+1] in thresholds and p[Cb]>thresholds[ channels[Cb+1] ]
		if posA and posB:
			seriesPP.add(p[Cb], p[Ca])
		elif posA:
			seriesPN.add(p[Cb], p[Ca])
		elif posB:
			seriesNP.add(p[Cb], p[Ca])
		else:
			seriesNN.add(p[Cb], p[Ca])
	dataset.addSeries(seriesNN)
	dataset.addSeries(seriesPN)
	dataset.addSeries(seriesNP)
	dataset.addSeries(seriesPP)
	
	chart = ChartFactory.createScatterPlot( title+" - "+channels[Cb+1]+" vs "+channels[Ca+1], channels[Cb+1], channels[Ca+1], dataset, PlotOrientation.VERTICAL, False,True,False )
	plot = chart.getPlot()
	plot.getDomainAxis().setRange(Range(0.00, maxIntensity), True, False)
	plot.getRangeAxis().setRange(Range(0.00, maxIntensity), True, False)
	renderer = chart.getPlot().getRenderer()
	
	renderer.setSeriesPaint(0, Color(64,64,64)) #NN
	renderer.setSeriesPaint(1, Color(0,255,0)) #PN
	renderer.setSeriesPaint(2, Color(0,0,255)) #NP
	renderer.setSeriesPaint(3, Color(0,255,255)) #PP

	shape = Ellipse2D.Float(-1,-1,3,3)
	renderer.setSeriesShape(0, shape )
	renderer.setSeriesShape(1, shape )
	renderer.setSeriesShape(2, shape )
	renderer.setSeriesShape(3, shape )
	
	frame = ChartFrame(title+" - "+channels[Cb+1]+" vs "+channels[Ca+1], chart)
	frame.setSize(800, 800)
	frame.setLocationRelativeTo(None)
	frame.setVisible(True)

def analyse():
	rt = ResultsTable()
	ol = Overlay()

	masks = [getMask(imp, c) for c in range(1, C+1)]
	DAPImask = masks[0]
	IJ.run(DAPImask, "Create Selection", "")
	DAPIRoi = DAPImask.getRoi()
	rois = ShapeRoi(DAPIRoi).getRois()
	for c,mask in enumerate(masks):
		if c==0: continue
		IJ.run(mask, "Create Selection", "")
		signalRoi = mask.getRoi()
		signalRoi.setPosition(c+1)
		colour = Color.BLUE
		if c==1: colour = Color.CYAN
		elif c==2: colour = Color.GREEN
		elif c==3: colour = Color.RED
		signalRoi.setStrokeColor(colour)
		ol.add(signalRoi)
	
	ip = [imp.getStack().getProcessor(c+1) for c in range(0, C)]
	points = []
	row = 0
	npos = [0 for i in range(C)]
	pos13 = 0
	pos12 = 0
	for roi in rois:
		ip[0].setRoi(roi)
		area = ImageStatistics.getStatistics(ip[0], Measurements.AREA, cal).area
		if area>=MINA and area<=MAXA:
			roi.setStrokeColor(Color.YELLOW)
			bounds = roi.getBounds()
			rt.setValue("X", row, (bounds.x+bounds.width/2)*cal.pixelWidth)
			rt.setValue("Y", row, (bounds.y+bounds.height/2)*cal.pixelHeight)
			rt.setValue("Area", row, area)
			stats = [i for i in range(0,C)]
			call = [False for i in range(C)]
			for c in range(1, C):
				ip[c].setRoi(roi)
				stats[c] = ImageStatistics.getStatistics(ip[c], Measurements.MEAN, cal)
				rt.setValue("C"+str(c+1)+" "+channels[c]+" Mean", row, stats[c].mean)
	
				masks[c].setRoi(roi)
				proportion = masks[c].getStatistics().mean/255
				rt.setValue("C"+str(c+1)+" "+channels[c]+" Proportion", row, proportion)

				if channels[c] in thresholds:
					thresh = thresholds[ channels[c] ]
					pos = stats[c].mean >= thresh
					if pos:
						npos[c] += 1
						call[c] = True
					rt.setValue("C"+str(c+1)+" "+channels[c]+" Positive?", row, str(pos))
			if len(call)>3 and call[1] and call[3]:
				pos13 += 1
			if len(call)>2 and call[1] and call[2]:
				pos12 += 1
				
			points.append( [s.mean for s in stats[1:]] )
			row += 1
			ol.add(roi)
	imp.setOverlay(ol)
	title = imp.getTitle()
	rt.show(title+" Results")

	if len(channels)>3:
		plot2D(points, 0, 2)
	if len(channels)>2:
		plot2D(points, 0, 1)
	results = ResultsTable.getResultsTable()
	cr = results.getCounter()
	results.setValue("Image", cr, title)
	results.setValue("total nuclei", cr, len(points))
	for c in range(1, C):
		results.setValue(channels[c]+" +ve", cr, npos[c])
	if len(channels) > 3:
		results.setValue(channels[1]+" and "+channels[3]+" +ve", cr, pos13)
	if len(channels) > 2:
		results.setValue(channels[1]+" and "+channels[2]+" +ve", cr, pos12)
	results.show("Results")

imp = WindowManager.getCurrentImage()
imp.killRoi()
C = imp.getNChannels()
cal = imp.getCalibration()
unit = cal.getUnit()
if re.match( "[Mm]icro.*", unit  ) is not None:
	unit = u"\u00b5"+"m"
Aunit = unit+u"\u00b2"

title = imp.getTitle()
chans = re.findall("\w+-\d{3}", title)
channels = ["" for c in range(0,C)]
channels[0] = "DAPI"
for chan in chans:
	split = chan.split("-")
	if split[1] == "488" and C>1:
		channels[1] = split[0]
	elif split[1] == "568" and C>2:
		channels[2] = split[0]
	elif split[1] == "647" and C>3:
		channels[3] = split[0]

analyse()

