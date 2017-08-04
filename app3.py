# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 16:27:34 2017

@author: Charlotte
"""
#python -c 'import foo; print foo.hello()'
import os, sys
import json
import StringIO
import numpy as np
np.random.seed(10)
#import requests
import pandas as pd

#import re
#from werkzeug.routing import Rule, RequestRedirect

from flask import Flask, request, redirect, render_template, send_from_directory, abort, Response
enter_point='/enrichr-fireworks' 
allcyjs=np.array(["Disease(noknn&nothresh25).gml.cyjs","TF(noknn&nothresh25)weightedthresh99.gml.cyjs","celltypejul31(noknn&nothresh25).gml.cyjs","Ontologyjul31(knn&thresh25).gml.cyjs","Pathwayjul31(noknn&nothresh25).gml.cyjs"])
allcyjs2=np.array(["Diseaserandom.json","TFpositionsnoknn99(0.2).txt","CellTypepositions.txt","Ontology(.001).txt","Path(real.001).txt"])
allmetadata=np.array(["Disease(noknn&nothresh25)names.txt","TF(noknn&nothresh25)weightedthresh99names.txt","celltypejul31(noknn&nothresh25)names.txt","Ontologyjul31(knn&thresh25)names.txt","Pathwayjul31(noknn&nothresh25)names.txt"])
app2 = Flask(__name__, static_url_path=enter_point, static_folder=os.getcwd())

app2.config['SEND_FILE_MAX_AGE_DEFAULT'] = 6


from orm3 import *

app2.debug=True


@app2.before_first_request
def load_globals2():
    global metadatalist, graph_df_list, graph_df_list_2, genelistnames, genesetlist, Res 
    Res=pd.DataFrame()
    cyjslist=[]
    cyjslist2=[]
    metadatalist=[]
    graph_df_list=[]
    graph_df_list_2=[]
    for i in range(allcyjs.size):
        cyjslist.append(allcyjs[i])
        cyjslist2.append(allcyjs2[i])
        metadatalist.append(pd.read_csv(allmetadata[i]))
    
    for i in range(len(cyjslist)):
        graph_df_list.append(load_graph(cyjslist[i],metadatalist[i]))
        graph_df_list_2.append(load_graph(cyjslist2[i],metadatalist[i]))
# taken from https://stackoverflow.com/questions/17714571/creating-a-dictionary-from-a-txt-file-using-python   
    genelistnames={1:'TF(noknn&nothresh25)weightedthresh99genes.txt',2:'celltypejul31(noknn&nothresh25)genes.txt',3:'Ontologyjul31(knn&thresh25)genes.txt',0:'Disease(noknn&nothresh25)genes.txt',4:'Pathwayjul31(noknn&nothresh25)genes.txt'}
    #genelistnames={0:'Diseasesgeneset.txt',1:'TFgeneset.txt',2:'CellTypegeneset.txt'}
    genesetlist={}
    #each key in genesetlist corresponds to a list(of all genesets) of lists(of genes)
    for (k,v) in genelistnames.items():
        genesetlist[k]=[]
        with open(v,'r') as f:
            for line in f:
                spl=line.split()
                li=[]
                for i in spl[1:]:
                    li.append(i)
                genesetlist[k].append(li)
    return 

@app2.route(enter_point+'/')
def index_page():
    return render_template('index.html',
        script='main3',
        enter_point=enter_point,
        result_id='hello')  

@app2.route(enter_point + '/graph/<string:graph>/<string:layout>', methods=['GET'])
def load_graph_layout_coords(graph,layout):
	if request.method == 'GET':
		index=int(graph)
		if layout=='cy':
			print graph_df_list[index].shape
			return graph_df_list[index].reset_index().to_json(orient='records')
		print graph_df_list_2[index].shape
		return graph_df_list_2[index].reset_index().to_json(orient='records')       

@app2.route(enter_point + '/graph/<string:layout>/<string:graph>', methods=['GET'])
def load_graph_layout_coords_own(layout,graph):
	if request.method == 'GET':
		index=int(graph)
		print graph_df_list_2[index].shape
		return graph_df_list_2[index].reset_index().to_json(orient='records')    

@app2.route(enter_point + '/sig_ids', methods=['GET'])
def get_all_sig_ids():
	if request.method == 'GET':
		cyjs_filename =CYJS #os.environ['CYJS']
		json_data = json.load(open('notebooks/%s' % cyjs_filename, 'rb'))
		json_data = json_data['elements']['nodes']
		sig_ids = [rec['data']['name'] for rec in json_data]
		return json.dumps({'sig_ids': sig_ids, 'n_sig_ids': len(sig_ids)})    
#
@app2.route(enter_point+'/<path:filename>')
def send_file(filename):
    return send_from_directory(app2.static_folder,filename)



@app2.route(enter_point + '/search', methods=['POST'])
def post_to_sigine():
	#Endpoint handling signature similarity search, POST the up/down genes 
	#to the RURL and redirect to the result page.
	if request.method == 'POST':
		up_genes = request.form.get('upGenes', '').split()
		gene_sets = GeneSets(up_genes)
		fisherresult = gene_sets.enrich(genesetlist)
		otherresult = gene_sets.enrichother(genesetlist)
		rid = gene_sets.save()
		print rid
		return redirect(enter_point + '/result/' + rid, code=302)
    
#@app2.route(enter_point+'/searchother', methods=['POST'])
#def post_to_sigineother():
#    if request.method=='POST':
#        up_genes=request.form.get('upGenes','').split()
#        gene_sets=GeneSets(up_genes)
#        result=gene_sets.enrichother(genesetlist)
#        #rid=gene_sets.saveanother()  #saves into same db as other
#        rid=gene_sets.save()
#        return redirect(enter_point+'/result/'+rid, code=302)




    
   #this gets called by the result javascript method
@app2.route(enter_point + '/result/<string:testtype>/<string:graph>/<string:layout>/<string:result_id>', methods=['GET'])
def result(testtype,graph, layout, result_id):
	"""
	Retrieve a simiarity search result using id and combine it
	with graph layout.
	"""
	index=int(graph)
	if layout=='cy':
		graph_df=graph_df_list[index]  
	else:
		graph_df=graph_df_list_2[index]
	# retrieve enrichment results from db
	result_obj = EnrichmentResult(result_id,testtype,index)
	# bind enrichment result to the network layout
	graph_df_res = result_obj.bind_to_graph(graph_df)

	return graph_df_res.reset_index().to_json(orient='records')


@app2.route(enter_point + '/result/<string:result_id>', methods=['GET'])
def result_page(result_id):
	#The result page.
	
	return render_template('index.html', 
		script='result', 
		enter_point=enter_point,
		result_id=result_id)

@app2.route(enter_point+'/inputgenes/<string:result_id>',methods=['GET'])
def gene_page(result_id):
    result_obj=EnrichmentResult(result_id,'fishertest',0)
    result_genes=result_obj.get_genes()
    return result_genes.to_json(orient='values')    
#    

@app2.route(enter_point+'/topn/<string:testtype>/<string:graph>/<string:result_id>',methods=['GET'])
def topn_page(testtype,graph,result_id):
    graph=int(graph)
    topn=EnrichmentResult(result_id,testtype,graph)
    topn=topn.get_topn(graph_df_list[graph]) 
    return topn.to_json(orient='records')

@app2.route(enter_point + '/result/download/<string:result_id>', methods=['POST','GET'])
def result_download(result_id):
	#To download the results to a csv file.
	s=StringIO.StringIO()
	filename=('enrichment_search_result_(%s).csv' % result_id)
	result_df={}
	frames=[]
	    
	for i in range(len(graph_df_list)):
		result_obj = EnrichmentResult(result_id,'fishertest',i)
    	# Prepare a DataFrame for the result
        #can we implement this for all

		graph_df=graph_df_list[i]
		scores = result_obj.result['score']
		result_df[i] = pd.DataFrame({'p_value': scores, 
	    	'geneset': graph_df['geneset'],
	    	'library': graph_df['library'],
	    	})\
			.sort_values('p_value', ascending=True)
	# Write into memory
		frames.append(result_df[i])
	# Prepare response
	final=pd.concat(frames)
	final.to_csv(s)	
	resp = Response(s.getvalue(), mimetype='text/csv')
	resp.headers['Content-Disposition'] = 'attachment; filename=%s' \
		% filename
	return resp



if __name__ == '__main__':
	app2.run(host='0.0.0.0',port=5000, threaded=True)
#host='0.0.0.0',