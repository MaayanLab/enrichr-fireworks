var textures = new Textures()


//changed /graph/ to number corresponding to correct category, the int will change with toggle
var sd = new ScatterData({
	resultid: result_id,
	currentgraph: currentgraph,
	currenttest: currenttest,
	currentlayout: currentlayout,
	url: 'result/'+currenttest+'/'+currentgraph+'/'+currentlayout+'/'+ result_id,
	textures: textures,
})

if(currentgraph==0)
	sdvgraph='Diseases and Drugs';
if(currentgraph==1)
	sdvgraph='Transcription';
if(currentgraph==2)
	sdvgraph='Cell Type';
if(currentgraph==3)
	sdvgraph='Pathway';
if(currenttest.localeCompare('othertest')==0)
	sdvtest='Chi Square';
else 
	sdvtest='Fishers Test';
if(currentlayout.localeCompare('cy')==0)
	sdvlayout='Allegro Edge Repulsive';
else
	sdvlayout='Basic Edge Repulsive';

var sdv = new Scatter3dView({
	model: sd,
	textures: textures,
	// pointSize: 0.1, 
	pointSize: 12,
	is3d: false,
	colorKey: 'score',
	shapeKey: 'library',
	labelKey:['geneset','library','score'],
	testtype:currenttest,
	layouttype: currentlayout,
	networkKey: sdvgraph,
	graphtype: currentgraph,
	layoutKey: sdvlayout,
	testkey: sdvtest

})

var overlay = new Overlay({scatterPlot: sdv})


var legend = new Legend({scatterPlot: sdv, h: window.innerHeight-200})

var topnscores = new Scores({scatterPlot:sdv})

var controler = new Controler({scatterPlot: sdv, h: window.innerHeight/2, w: 200})

var sigSimSearch = new SigSimSearchForm({scatterPlot: sdv, container: "#controls1", result_id: result_id})

//controler.render();
//var resultModalBtn = new ResultModalBtn({scatterPlot: sdv, container: document.body, result_id: result_id})

//var resultModal = new ResultModal({scatterPlot: sdv});

