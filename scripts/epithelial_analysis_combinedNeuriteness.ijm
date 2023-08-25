
// ===============================
// AUTHOR : Gisele Miranda 
// IVESTIGATORS : Annelie Tjernlund, Mathias Franzén Boger, Gabriela Edfeldt, Kristina Broliden
// CREATE DATE : 2023 - 05 - 10
// PURPOSE : Structural analysis of the cervical epithelial tissue
// NOTES : Required plugins - MorphoLibJ (https://imagej.net/MorphoLibJ)
// ===============================


/********************** parameters ***********************/

// path to the input directory containing the SNIs
#@ File (label = "SNI directory", style = "directory") path
#@ File (label = "MAX directory", style = "directory") maxProjPath

useScale = true;
pixPerMic = 3.07693;
area_connected_threshold = 5000;
/*********************************************************/

// get parameters per channel
exp_id_cy3  = getParametersPerChannel("cy3"); 
exp_id_cy5  = getParametersPerChannel("Cy5"); 
exp_id_fitc = getParametersPerChannel("FITC");

dir = getFileList(path);
bufferMeasures = "";
overlayPath = path + File.separator + "combinedOverlay" + File.separator;
if(!File.exists(overlayPath)) File.makeDirectory(overlayPath);

for(countSNI=0; countSNI<dir.length; countSNI++) { // for each SNI
	print("countSNI: " + countSNI);
	sni = dir[countSNI];
	images = getFileList(path+sni);
	
	// check if input folder is a SNI. It should start with "SNI" or "Neg SNI" (for negative control)
	if(startsWith(sni, "SNI") || startsWith(sni, "Neg SNI")) {
		
		// cread output folder for combined neuriteness
		outpath = path + File.separator + sni + "combinedOutput" + File.separator;
		if(!File.exists(outpath)) File.makeDirectory(outpath);
		
		sni = replace(sni, "/", "");
		print(sni);

		// check if experiment folders exist
		exp_output_cy3 = path + File.separator + sni + File.separator + "output_" + exp_id_cy3 + "/";
		exp_output_cy5 = path + File.separator + sni + File.separator + "output_" + exp_id_cy5 + "/";
		exp_output_fitc = path + File.separator + sni + File.separator + "output_" + exp_id_fitc + "/";
		
		print("exp_output_cy3: " + exp_output_cy3);
		print("exp_output_cy5: " + exp_output_cy5);
		print("exp_output_fitc: " + exp_output_fitc);
		
		bufferMeasures = bufferMeasures + sni + ";";
		
		// check if all output folders exist
		if(File.exists(exp_output_cy3) && File.exists(exp_output_cy5) && File.exists(exp_output_fitc)) {
			
			// get combined neuriteness
			files = newArray("", "", "");
			neur_cy3_path  = exp_output_cy3 + "thresholded_neuriteness_segmented_region_cy3.tif";
			if(File.exists(neur_cy3_path)) files[0] = neur_cy3_path;
			neur_cy5_path  = exp_output_cy5 + "thresholded_neuriteness_segmented_region_cy5.tif";
			if(File.exists(neur_cy5_path)) files[1] = neur_cy5_path;
			neur_fitc_path = exp_output_fitc + "thresholded_neuriteness_segmented_region_FITC.tif";
			if(File.exists(neur_fitc_path)) files[2] = neur_fitc_path;
			
			combID_neur = getCombinedMasks(files, outpath, "combined_neuriteness.tif");
			
			// get combined segmentation 
			files = newArray("", "", "");
			seg_cy3_path  = exp_output_cy3 + "segmented_region_cy3.tif";
			if(File.exists(seg_cy3_path)) files[0] = seg_cy3_path;
			seg_cy5_path  = exp_output_cy5 + "segmented_region_cy5.tif";
			if(File.exists(seg_cy5_path)) files[1] = seg_cy5_path;
			seg_fitc_path = exp_output_fitc + "segmented_region_FITC.tif";
			if(File.exists(seg_fitc_path)) files[2] = seg_fitc_path;
			
			combID_seg = getCombinedMasks(files, outpath, "combined_segmentation.tif");
			
			if(combID_neur == 1 && combID_seg == 1) { // then combined masks exists
				
				// get areas of combined neuriteness and combined segmentation
				comb_neur_area = getCombinedArea("combined_neuriteness.tif", outpath);
				comb_seg_area = getCombinedArea("combined_segmentation.tif", outpath);
				
				// get total area of the epithelium
				tot_area_epithelium = getCombinedArea("mask_epithelium.tif", exp_output_cy3);
				
				// open all masks to get upper, lower and parabasal layers
				open(outpath + "combined_neuriteness.tif");
				rename("combined_neuriteness");
				open(outpath + "combined_segmentation.tif");
				rename("combined_segmentation");
				
				open(exp_output_cy3 + "mask_epithelium.tif"); // open mask epithelium
				rename("mask_epithelium");
				
				open(exp_output_cy3 + "mask_nuclei.tif"); // open nuclei mask
				rename("mask_nuclei");
				
				run("Duplicate...", "temp_ROIs");
				rename("temp_ROIs");

				// get apical and basal lines
				roiManager("Open", maxProjPath + File.separator + sni + "_AB.zip");
				roiManager("Show All");
				createLayer(0, exp_output_cy3);
				selectWindow("temp_ROIs");
				createLayer(1, exp_output_cy3);
				selectWindow("ROI Manager"); 
				run("Close");
				close("temp_ROIs");
				
				// get area of parabasal layer by adding the mask_channel to the mask_nuclei
				imageCalculator("Add create", "combined_segmentation","mask_nuclei");
				selectWindow("Result of combined_segmentation");
				rename("parabasal_marker");
				run("Fill Holes");

				selectWindow("parabasal_marker");
				if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
				else run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
				run("Set Measurements...", "area mean standard min redirect=None decimal=4");
				run("Analyze Particles...", "clear include summarize");
				Table.rename("Summary", "Results");
				tot_area_parabasal = getResult("Total Area", 0);
				selectWindow("Results"); 
				run("Close");
				
				// get area of upper layer (between parabasal and superficial layer)
				imageCalculator("Subtract create", "mask_epithelium","parabasal_marker");
				selectWindow("Result of mask_epithelium");
				run("Analyze Particles...", "pixel show=Nothing display clear include add");
				big_ind = findBiggestROI();
				roiManager("select", big_ind);
				run("Create Mask");
				rename("upper");
				selectWindow("Results"); 
				run("Close");
				selectWindow("ROI Manager"); 
				run("Close");
				close("Result of mask_epithelium");
				selectWindow("upper");
				
				if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
				else run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
				run("Analyze Particles...", "clear summarize");
				Table.rename("Summary", "Results");
				tot_area_upper = getResult("Total Area", 0);
				selectWindow("Results"); 
				run("Close");
				
				// get area of lower layer (between upper + segmented marker and basal layer)
				selectWindow("combined_segmentation");
				imageCalculator("Add create", "upper","combined_segmentation");
				run("Fill Holes");
	
				imageCalculator("Subtract create", "mask_epithelium","Result of upper");
				rename("lower");
				if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
				else run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
				run("Analyze Particles...", "clear include summarize");
				Table.rename("Summary", "Results");
				tot_area_lower = getResult("Total Area", 0);
				selectWindow("Results"); 
				run("Close");

				selectWindow("parabasal_marker");
				saveAs("Tif", outpath + "parabasal_layer.tif");
				selectWindow("upper");
				saveAs("Tif", outpath + "upper_layer.tif");
				selectWindow("lower");
				saveAs("Tif", outpath + "lower_layer.tif");
				
				// distance transform of segmented channel marker to apical layer
				distA = calculateDistance("combined_segmentation", "apical", false);
				close("EDT-mask_channel");
				selectWindow("apical");
				run("Remove Overlay");
				
				// distance transform of segmented channel marker to basal layer
				distB = calculateDistance("combined_segmentation", "basal", false);
				close("EDT-mask_channel");
				selectWindow("basal");
				run("Remove Overlay");
				run("Close All");
				
				// get intact net measures
				intact_net_measures = 0;
				flooded_area = 0;
				
				intact_net_measures = getIntactNet_ConnectedComponent("", sni, outpath, overlayPath);
				//this functions also returns the area, however we're not using it here now
				getIntactNet(sni, outpath, overlayPath);
				// getFloodedArea depends on the results generated on function getIntactNet_ConnectedComponent
				flooded_area = getFloodedArea(outpath, exp_output_cy3, overlayPath, maxProjPath + File.separator + sni + "_AB.zip", sni);
				if (flooded_area == -1) flooded_area = tot_area_epithelium; // then flooded_area will be equal the area of the epithelium
				
				intact_net_measures = split(intact_net_measures," ");
				area_intact_net_cc_holes = parseFloat(intact_net_measures[0]);
				area_intact_net_cc = parseFloat(intact_net_measures[1]);
				
				// copy measures to buffer
				bufferMeasures = bufferMeasures + tot_area_epithelium + ";" + comb_seg_area + ";" + comb_neur_area + ";" + 
								 (comb_seg_area/tot_area_epithelium) + ";" + (comb_neur_area/tot_area_epithelium) + ";" + 
								 (comb_neur_area/comb_seg_area) + ";" + tot_area_parabasal + ";" + tot_area_upper + ";" +
								 tot_area_lower + ";" + distA + ";" + distB + ";" + area_intact_net_cc + ";" + 
								 area_intact_net_cc_holes + ";" + (area_intact_net_cc/comb_neur_area) + ";" + 
								 flooded_area + ";" + (flooded_area/tot_area_epithelium) + "\n";
			} else {
				print("not possible to generate combined masks");
			}
		} else {
			print("missing output folders per channel");
		}
	}
}

summaryPath = path + "/summary_combined.csv";

if(File.exists(summaryPath))
	File.delete(summaryPath);

summaryFile = File.open(summaryPath);
print(summaryFile, "sni;area epithelium;area combined segmentation;area combined neuriteness;relative area - area combined segmentation per total epithelium;realtive area - area combined neuriteness per total epithelium;realtive area - area combined neuriteness per area combined segmentation; area - marker of interest + parabasal layer;area - upper layer;area - lower layer;average height between marker of interest and apical layer;std height between marker of interest and apical layer;min height between marker of interest and apical layer;max height between marker of interest and apical layer;average height between marker of interest and basal layer;std height between marker of interest and basal layer;min height between marker of interest and basal layer;max height between marker of interest and basal layer;area - intact net;area - intact net with holes;relative area - intact net per neuriteness total area;flooded area;relative area - flooded area per total epithelium\n");
print(summaryFile, bufferMeasures);
File.close(summaryFile);

/****************** functions ****************************/

function getParametersPerChannel(ch) { 
	title = "Input parameters for " + ch;
	
	Dialog.create("Input parameters for " + ch + " channel");
	Dialog.addString("Threshold method:", "Otsu");
	Dialog.addString("Correction factor:", "1");
	Dialog.addString("Distance threshold:", "7");
	Dialog.addString("Area threshold:", "70000");
	Dialog.addString("Radius dilation:", "5");
	Dialog.addString("Experiment ID:", "id");
	Dialog.show();
	
	channel_config = Dialog.getString() + "_correc-factor-" + Dialog.getString() + "_dist-" + Dialog.getString() + "_area-" + Dialog.getString() + "_radDilation-" + Dialog.getString() + "_" + Dialog.getString();
	return channel_config;
}

function getCombinedMasks(files, out_path, mask_name) {
	comb = 0;
	
	if(matches(files[0], "") && matches(files[1], "") && matches(files[2], "")) { // no neuriteness found
		comb = 0;
		print("all empty");
	}
	
	if(!matches(files[0], "") && !matches(files[1], "") && matches(files[2], "")) { // combine cy3 and cy5
		open(files[0]);
		rename("cy3");
		open(files[1]);
		rename("cy5");
		
		mergeMasks("cy3","cy5");
		comb = 1;
		print("cy3 and cy5");
	}
	
	if(!matches(files[0], "") && matches(files[1], "") && !matches(files[2], "")) { // combine cy3 and fitc
		open(files[0]);
		rename("cy3");
		open(files[2]);
		rename("fitc");
		
		mergeMasks("cy3","fitc");
		comb = 1;
		print("cy3 and fitc");
	}
	
	if(matches(files[0], "") && !matches(files[1], "") && !matches(files[2], "")) { // combine cy5 and fitc
		open(files[1]);
		rename("cy5");
		open(files[2]);
		rename("fitc");
		
		mergeMasks("cy5","fitc");
		comb = 1;
		print("cy5 and fitc");
	}
	
	if(!matches(files[0], "") && matches(files[1], "") && matches(files[2], "")) { // return cy3
		open(files[0]);
		rename("merged");
		comb = 1;
		
		print("cy3");
	}
	
	if(matches(files[0], "") && !matches(files[1], "") && matches(files[2], "")) { // return cy5
		open(files[1]);
		rename("merged");
		comb = 1;
		
		print("cy5");
	}
	
	if(matches(files[0], "") && matches(files[1], "") && !matches(files[2], "")) { // return fitc
		open(files[2]);
		rename("merged");
		comb = 1;
		
		print("fitc");
	}
	
	if(!matches(files[0], "") && !matches(files[1], "") && !matches(files[2], "")) { // combine all
		print("all");
		open(files[0]);
		rename("cy3");
		open(files[1]);
		rename("cy5");
		open(files[2]);
		rename("fitc");
		
		mergeMasks("cy3", "cy5");
		rename("merged_cy3_cy5");
		mergeMasks("merged_cy3_cy5", "fitc");
		
		comb = 1;
	}
	
	if(comb >= 0) saveAs("Tif", out_path + mask_name);
	run("Close All");
	return comb;
}

function mergeMasks(img_1, img_2) {
	imageCalculator("Add create", img_1, img_2);
	rename("merged");
}

function createLayer(ind, out_path) {
	roiManager("Select", ind);
	rName = Roi.getName();
	run("Create Mask");
	
	open(out_path + "mask_epithelium.tif");
	rename("epithelium");
	run("Create Selection");
	selectWindow("Mask");
	run("Restore Selection");
	run("Clear Outside");
	run("Select None");

	run("Set Measurements...", "area mean standard min redirect=None");
	run("Set Scale...", "distance=1 known=1 unit=unit");
	run("Analyze Particles...", "size=500-Infinity show=Masks");
	setOption("BlackBackground", true);
	run("Convert to Mask");

	if(rName == "A") rename("apical");
	else rename("basal");

	close("Mask");
	close("epithelium");
}

function calculateDistance(refLayer, targetLayer, areaValue) {
	selectWindow(refLayer);
	run("Invert");
	run("Exact Euclidean Distance Transform (3D)");
	rename("EDT-"+refLayer);
	
	selectWindow(targetLayer);
	run("Set Measurements...", "area mean standard min redirect=EDT-"+refLayer+" decimal=4");

	run("Analyze Particles...", "clear add");
	selectWindow("EDT-"+refLayer);
	run("Conversions...", "weighted");
	setOption("ScaleConversions", false);
	run("16-bit");
	run("From ROI Manager");

	// get number of pixels in reference layer
	selectWindow("EDT-"+refLayer);
	roiManager("Select", 0);
	roiManager("Measure");
	area = getResult("Area", 0);
	selectWindow("Results"); 
	run("Close");
	
	// get the other measures: average, std, min and max
	if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
	selectWindow("EDT-"+refLayer);
	roiManager("Select", 0);
	roiManager("Measure");

	if(useScale) {
		mean = getResult("Mean", 0)/pixPerMic;
		std = getResult("StdDev", 0)/pixPerMic;
		min = getResult("Min", 0)/pixPerMic;
		max = getResult("Max", 0)/pixPerMic;
	} else {
		mean = getResult("Mean", 0);
		std = getResult("StdDev", 0);
		min = getResult("Min", 0);
		max = getResult("Max", 0);
	}

	selectWindow(refLayer);
	run("Invert");

	selectWindow("Results"); 
	run("Close");
	selectWindow("ROI Manager"); 
	run("Close"); 

	buffer = "";
	if(areaValue)
		buffer = "" + area + ";" + mean + ";" + std + ";" + min + ";" + max;
	else
		buffer = "" + mean + ";" + std + ";" + min + ";" + max;

	return buffer;
} 

// returns the area of the entire neuriteness network
function getCombinedArea(img, out_path) {
	
	open(out_path + img);
	rename("ref_mask");
	
	if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
	run("Set Measurements...", "area mean standard min redirect=None decimal=4");
	run("Analyze Particles...", "clear summarize");
	Table.rename("Summary", "Results");
	tot_area = getResult("Total Area", 0);

	selectWindow("Results"); 
	run("Close");
	
	run("Close All");
	
	return(tot_area);
}

function getIntactNet_ConnectedComponent(channel, sni, out_path, overlay_path) {
	// open combined neuriteness image
	open(out_path + "combined_neuriteness.tif");
	if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
		
	rename("inputImg");
	run("Duplicate...", "markerImg");
	rename("markerImg");
	
	// filter connected components by size and create filtered mask
	run("Set Measurements...", "area mean standard min redirect=markerImg decimal=4");
	run("Analyze Particles...", "size="+area_connected_threshold+"-Infinity show=Masks clear include add"); // include holes
	count = roiManager("count");
	selectWindow("ROI Manager"); 
	run("Close");
	
	if(count != 0) {
		run("Convert to Mask");
		run("Create Selection");
		selectWindow("markerImg");
		run("Restore Selection");
		run("Clear Outside");
		run("Select None");
	
		selectWindow("Mask of markerImg");
		//if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
		run("Measure");
		area_with_holes = getResult("Area", 0);
		selectWindow("Results"); 
		run("Close");

		// measure MFI of filtered connected components and save LUT image (glasbey) - not taking into account holes in the binary mask
		selectWindow("markerImg");
		run("Duplicate...", "markerImg_intactNet");
		rename("markerImg_intactNet");
		
		run("Set Measurements...", "area mean standard min redirect=inputImg decimal=4");
		run("Analyze Particles...", "size="+area_connected_threshold+"-Infinity show=[Count Masks] clear");
		run("glasbey");
		saveAs("PNG", out_path + "intact_net_connected.png");
		saveAs("PNG", overlay_path + "intact_net_connected_" + sni);
		
		selectWindow("markerImg");
		run("Create Selection");
		selectWindow("inputImg");
		run("Restore Selection");
		run("Measure");
		
		if(nResults > 0) {
			tot_area = getResult("Area", 0);
			// concatenate both area measures
			m = d2s(area_with_holes,4) + " " + d2s(tot_area,4);
		} 
	} else m = "0 0";
	
	run("Close All");
	return m;
}

function getIntactNet(sni, out_path, overlay_path) {

	open(out_path + "combined_neuriteness.tif");
		
	run("Median...", "radius=4");
	rename("inputImg");
	run("Duplicate...", "markerImg");
	rename("markerImg");
	run("Invert");
	
	run("Marker-controlled Watershed", "input=inputImg marker=markerImg binary calculate use");
	run("16-bit");
	setThreshold(2, 65535);
	run("Convert to Mask");
	saveAs("Tif", out_path + "watershed.tif");
	saveAs("Tif", overlay_path + "watershed_" + sni);
	rename("watershed");

	run("Morphological Filters", "operation=Dilation element=Square radius=6");
	imageCalculator("AND create", "inputImg","watershed-Dilation");
	if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
	run("Set Measurements...", "area mean standard min redirect=None decimal=4");
	selectWindow("Result of inputImg");
	run("Analyze Particles...", "clear summarize"); // get area of the intact net
	Table.rename("Summary", "Results");
	tot_area = getResult("Total Area", 0);

	selectWindow("Results"); 
	run("Close");
	
	saveAs("Tif", out_path + "intact_net_watershed.tif");
	saveAs("Tif", overlay_path + "watershed_" + sni);

	run("Close All");
	return(tot_area);
}

function getFloodedArea(out_path, out_path_channel, overlay_path, maxProj_path, sni) {
	open(out_path_channel + "mask_epithelium.tif");
	
	wid = getWidth();
	hei = getHeight();
	newImage("temp", "8-bit black", wid, hei, 1);
	
	roiManager("Open", maxProj_path);
	roiManager("Show All");
	
	roiManager("Select", 0);
	run("Create Mask");
	rename("apical");
	
	selectWindow("temp");
	roiManager("Select", 1);
	run("Create Mask");
	rename("basal");
	
	selectWindow("temp");
	close();
	selectWindow("ROI Manager"); 
	run("Close");

	if(File.exists(out_path + "intact_net_connected.png")) {
		open(out_path + "intact_net_connected.png");
		setThreshold(1, 65535, "raw");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		rename("neuriteness");
			
		selectWindow("mask_epithelium.tif");
		run("Duplicate...", "epithelium_outline");
		rename("epithelium_outline");
		run("Outline");
		
		// create marker for whatershed
		selectWindow("apical");
		run("Morphological Filters", "operation=Dilation element=Square radius=20");
		rename("marker");
		imageCalculator("AND create", "mask_epithelium.tif","marker");
		rename("marker_processed");
		run("Erode");
		
		imageCalculator("Add create", "epithelium_outline","neuriteness");
		rename("mask_watershed");
		run("Invert");
		run("Morphological Reconstruction", "marker=marker_processed mask=mask_watershed type=[By Dilation] connectivity=4");
		
		saveAs("Tif", overlay_path + "flooded_region_" + sni);
		saveAs("Tif", out_path + "flooded_region.tif");
		rename("marker_processed-rec");
		
		if(useScale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
		run("Set Measurements...", "area mean standard min redirect=None decimal=4");
		run("Analyze Particles...", "clear summarize");
		Table.rename("Summary", "Results");
		area_flooded = getResult("Total Area", 0);
		selectWindow("Results"); 
		run("Close");
	} else area_flooded = -1;
	
	run("Close All");
	return(area_flooded);
}

function findBiggestROI() {
	biggestInd = 0; // index of the first biggest ROI
	biggestROI = 0;
	
	for (i=0; i<nResults; i++) {
		area_i = getResult("Area", i);
		if(area_i > biggestROI) {
			biggestROI = area_i;
			biggestInd = i;
		}
	}

	return biggestInd;
}

