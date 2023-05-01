///*
//* TODO documentation for this script
//*/

// Get directories from the user via the GUI
#@ File (label = "Input directory", style = "directory") raw
#@ File (label = "Output directory", style = "directory") out
#@ File (label = "Path to modeldef", style = "file") modeldef
#@ String (label = "Caffemodel name") caffemodel

// call to function
processFolder();

// function to do segmentation on <count> images in the test set (-1 for all)
// note: make sure the caffemodel file is in home directory or get a cryptic error
// note: seemingly need to connect to the RSA keyfile once manual before this script will work
function processFolder() {
	// make a sorted list of raw images
	raw_list = getFileList(raw);
	count = raw_list.length;//TODO1 replace this with raw_list.length for a full segmentation
	if (count == -1){
		count = raw_list.length;
	}
	for (i = 0; i < count; i++) {
		if(endsWith(raw_list[i], ".png")){
		// open raw image
		    open(raw + "/" + raw_list[i]);
		    call('de.unifreiburg.unet.SegmentationJob.processHyperStack', 'modelFilename=' + modeldef + ',Memory (MB):=4000,weightsFilename=' + caffemodel + ',gpuId=GPU 0, useRemoteHost=false, processFolder=process-folder/, average=rotate, keepOriginal=true, outputScores=false,outputSoftmaxScores=false');
		    saveAs("PNG", out + File.separator + raw_list[i]);//changed the extra .png problem
		    run("Close All");//new
		}
	 }
	 run("Close All");
}
