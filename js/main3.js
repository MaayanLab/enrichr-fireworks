var textures = new Textures()


var sd = new ScatterData({

	url: 'graph/0/cy',
	textures: textures,
})
//container is set so that mouse hover in View is accurate when header is in place
var container=document.getElementById('body')
var width= container.clientWidth;
var height=container.clientHeight;

var sdv = new Scatter3dView({
	container: container,
	WIDTH: width,
	HEIGHT: height,
	model: sd,
	textures: textures,
	// pointSize: 0.1, 
	pointSize: 12,
	is3d: false,
	colorKey: 'library',
	shapeKey: 'library',
	labelKey: ['geneset','library'],
})
var legend = new Legend({scatterPlot: sdv, h: window.innerHeight, container: container})

var controler = new Controler({scatterPlot: sdv, h: window.innerHeight/2, w: 200, container: container})



var sigSimSearch = new SigSimSearchForm({scatterPlot: sdv, container: "#controls1"})

var overlay = new Overlay({scatterPlot: sdv})

