var textures = new Textures()


var sd = new ScatterData({

	url: 'graph/0/cy',
	textures: textures,
})


var sdv = new Scatter3dView({
	model: sd,
	textures: textures,
	// pointSize: 0.1, 
	pointSize: 12,
	is3d: false,
	colorKey: 'library',
	shapeKey: 'library',
	labelKey: ['geneset','library'],
})
var legend = new Legend({scatterPlot: sdv, h: window.innerHeight})

var controler = new Controler({scatterPlot: sdv, h: window.innerHeight-800, w: 200})

//var search = new SearchSelectize({scatterPlot: sdv, container: "#controls"})


var sigSimSearch = new SigSimSearchForm({scatterPlot: sdv, container: "#controls1"})

var overlay = new Overlay({scatterPlot: sdv})

