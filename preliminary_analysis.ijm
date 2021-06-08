/************* parameters *************/
path = "/Users/gisele.miranda/Desktop/neuritness 17 nov 2020/Missed ROIs LAMIQ/";
maxProjPath = "/Users/gisele.miranda/Desktop/neuritness 17 nov 2020/Missed ROIs LAMIQ_MAX_projection/";
nBins = 70; // number of histogram bins
increment = 500;
channel = "FITC";
/**************************************/

bufferSummary = "image;area;mean;std;min;max;median\n";
bufferHist = "image;";

startBin = 0;
for(binCount=0; binCount<nBins; binCount++) {
	if(binCount == (nBins-1)) {
		bufferHist = bufferHist + "[" + startBin + ", " + (startBin+increment) + "[\n";
	} else {
		bufferHist = bufferHist + "[" + startBin + ", " + (startBin+increment) + "[;";
	}
	startBin = startBin + increment;
}

function euclideanDist(x0,y0,x1,y1) {
	return Math.sqrt(Math.pow((x1-x0), 2) + Math.pow(y1-y0, 2));
}

/************* analyze each file of the selected folder *************/
SNIs = getFileList(path);
for(cont=0; cont<SNIs.length; cont++) {
	folder = SNIs[cont];
	images = getFileList(path+folder);
	
	if(startsWith(folder, "SNI") || startsWith(folder, "Neg SNI")) { 
		channelImg = false;
		// open all channels and search for the query channel
		for(cont2=0; cont2<images.length; cont2++) {
			image = images[cont2];
			
			if(matches(image, ".*" + channel + ".*") && !startsWith(image, "Ne_")) { // check if file is the query channel
				channelImg = true;
				run("Bio-Formats Importer", "open=[" + path + folder + image + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
				rename("inputImg");

				print(path+folder+image);

				// read original FITC channel again and perform segmentation with the distance transform
				wid = getWidth();
				hei = getHeight();

				// generate epithelium contour + apical and basal layers to crop the epithelium region
				newImage("Temp", "8-bit black", wid, hei, 1);
				print(maxProjPath + replace(folder, "/", "") + "_AB.zip");
				roiManager("Open", maxProjPath + replace(folder, "/", "") + "_AB.zip");
				roiManager("Show All");
							
				roiManager("Select", 0);
				Roi.getCoordinates(x1, y1);
				run("Create Mask");								
				selectWindow("Temp");
				roiManager("Select", 1);				
				Roi.getCoordinates(x2, y2);
				run("Create Mask");

				// check the distance betweeen the ending points of each line and connect the points with minimum distance
				dA = euclideanDist(x1[0],y1[0],x2[0],y2[0]);
				dB = euclideanDist(x1[0],y1[0],x2[(x2.length-1)],y2[(y2.length-1)]);
				
				selectWindow("Mask");
				if(dA < dB) {
					makeLine(x1[0], y1[0], x2[0], y2[0]);
					run("Create Mask");
					makeLine(x1[(x1.length-1)], y1[(y1.length-1)], x2[(x2.length-1)], y2[(y2.length-1)]);
					run("Create Mask");
				} else {
					makeLine(x1[0], y1[0], x2[(x2.length-1)], y2[(y2.length-1)]);
					run("Create Mask");
					makeLine(x1[(x1.length-1)], y1[(y1.length-1)], x2[0], y2[0]);
					run("Create Mask");
				}
				run("Fill Holes");
				run("Select None");

				selectWindow("ROI Manager"); 
				run("Close");
				close("Temp");

				selectWindow("Mask");
				rename("Mask-epithelium");
				run("Create Selection");

				selectWindow("inputImg");
				run("Restore Selection");
				run("Clear Outside");

				run("Set Measurements...", "area mean standard median min redirect=None decimal=3");
				run("Measure");
				
				bufferSummary = bufferSummary + replace(folder, "/", "") + ";" + getResultString("Area", 0) + ";" + getResultString("Mean", 0) + ";" + getResultString("StdDev", 0) + ";" + getResultString("Min", 0) + ";" + getResultString("Max", 0) + ";" + getResultString("Median", 0) + "\n";
				selectWindow("Results"); 
				run("Close");

				selectWindow("inputImg");
				getHistogram(values, counts, nBins, 0, (nBins*increment));
				line = ";";
				for(i=0; i<nBins; i++) {
					if(i == (nBins-1)) {
						line = line + counts[i] + "\n";
					}
					else line = line + counts[i] + ";";
				}
				bufferHist = bufferHist + replace(folder, "/", "") + line;

				run("Close All");
			}
		}
		
		if(!channelImg) {
			bufferSummary = bufferSummary + replace(folder, "/", "") + "\n";
			bufferHist = bufferHist + replace(folder, "/", "") + "\n";
		}
	}
}

// create file that contains the histograms of red and green channel
summaryPath = path+"preliminary_statistics_" + channel + "_channel.csv";
summaryPathHist = path+"preliminary_statistics_histogram_" + channel + "_channel.csv";

if(File.exists(summaryPath))
	File.delete(summaryPath);

summaryFile = File.open(summaryPath);
print(summaryFile, bufferSummary);
File.close(summaryFile);

if(File.exists(summaryPathHist))
	File.delete(summaryPathHist);

summaryFileHist = File.open(summaryPathHist);
print(summaryFileHist, bufferHist);
File.close(summaryFileHist);
