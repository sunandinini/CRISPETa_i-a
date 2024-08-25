#!/usr/bin/env python
import re, sys, os
from Bio.Seq import Seq
import subprocess

import pandas as pd
import csv, argparse, sys
import pickle
import PyPDF2

global model

import argparse, tempfile, os, re, math, sys, resource, time, datetime
from subprocess import call
from os import listdir
from PyPDF2 import PdfFileMerger


#Biopithon
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from operator import itemgetter

#graphics
import numpy as np
try:
	import chart_studio
	import plotly
	
	_plotly=1
except ImportError:
	_plotly=0
if _plotly==1:
	import chart_studio.plotly as py
	import plotly.graph_objs as go
	from plotly.graph_objs import Scatter, Layout

# read input and output files from command line
def plots(file_name, results_regions):
	#if not os.path.exists(file_name+'_graphics'):
	#	os.makedirs(file_name+'_graphics')
	if not os.path.isfile(file_name):
		return 
	f = open(file_name,'r')

	pairs = 0
	regions = 0
	switch = 0
	n_pairs = []
	score1 = []
	score2 = []
	pscore = []
	dist = []
	dist1 = []
	dist2 = []
	diff_dist = []
	picking_round = []
	previous_ID = ''
	for line in f:
		if line.startswith('target_id'):
			continue
		row = line.strip().split('\t')
		current_ID = row[0].split(':')[0]
		
		if current_ID == previous_ID:
			pairs += 1
		else:
			if switch == 1:
				n_pairs.append(int(pairs))
			regions += 1
			pairs = 1
		
		switch = 1
		score1.append(float(row[11]))
		score2.append(float(row[25]))
			
		pscore.append(float(row[11]) + float(row[25]))
		dist1.append(int(row[12]))
		dist2.append(int(row[26]))
		if (row[3] < row[16]):
			diff_dist.append(abs(int(row[3])-int(row[16])))
		else:
			diff_dist.append(abs(int(row[2])-int(row[17])))
		previous_ID = current_ID
		
	dist = dist1+dist2
	n_pairs.append(int(pairs))
	iscore=score1+score2
	dir = './plots/'
	
	fig = def_fig(n_pairs+[0]*results_regions[0], 'Number of pairs designed for each target region', '# pairs', 'Frequency')
	html1 = plotly.offline.plot(fig, auto_open=False, show_link=False, link_text='',output_type='div')
	pdf1 = plotly.io.write_image(fig, dir+'pdf1.pdf')
	

	fig = def_fig(iscore, 'Individual sgRNA scores', 'scores', 'Frequency') ##Individual sgRNA scores
	html2 = plotly.offline.plot(fig, auto_open=False, show_link=False, link_text='',output_type='div')
	pdf2 = plotly.io.write_image(fig, dir+'pdf2.pdf')


	fig = def_fig(pscore, 'Paired sgRNA scores', 'scores', 'Frequency') ##Paired sgRNA scores
	html3 = plotly.offline.plot(fig, auto_open=False, show_link=False, link_text='',output_type='div')
	pdf3 = plotly.io.write_image(fig, dir+'pdf3.pdf')


	fig = def_fig(dist, 'Distance from TSS', 'distance (nt)', 'Frequency') ##Distance of pairs from TSS
	html4 = plotly.offline.plot(fig, auto_open=False, show_link=False, link_text='',output_type='div')
	pdf4 = plotly.io.write_image(fig, dir+'pdf4.pdf')
	
	fig = def_fig(diff_dist, 'Distance between pairs', 'distance (nt)', 'Frequency') ##Distance between pairs
	html5 = plotly.offline.plot(fig, auto_open=False, show_link=False, link_text='',output_type='div')
	pdf5 = plotly.io.write_image(fig, dir+'pdf5.pdf')
	
	chartpie = def_chartpie(n_pairs+[0]*results_regions[0])
	html6 = plotly.offline.plot(chartpie, auto_open=False, show_link=False, link_text='',output_type='div')
	pdf6 = plotly.io.write_image(chartpie, dir+'pdf6.pdf')
	


	
	#Create html file with all images using 4 htmml codes (1 code from each hist)
	html_string = '\n'.join([html2,html3,html4,html1,html5,html6])
	fo = open(dir+'images.html','w+')
	fo.write('\n'.join([html1,html6,html2,html3,html4,html5]))
	fo.close()
	
	pdfs = [dir+'pdf1.pdf', dir+'pdf6.pdf', dir+'pdf2.pdf', dir+'pdf3.pdf', dir+'pdf4.pdf', dir+'pdf5.pdf']

	merger = PdfFileMerger()

	for pdf in pdfs:
		merger.append(pdf)

	merger.write(dir+'_compiled.pdf')
	merger.close()
	

	#Converts to pdf
	#try:
	#	import pdfkit
	#	_pdfkit=1
	#except ImportError:
	#	_pdfkit=0
	#	print "\n\tpdfkit module not installed. Install it to obtain pdf with graphics."
	
	#if _pdfkit==1:
	#	config1 = pdfkit.configuration(wkhtmltopdf="/usr/local/bin/wkhtmltopdf")
	#	pdfkit.from_string(html_string,file_name+'images.pdf',options={'quiet':''}, configuration=config1)


def def_chartpie(values):

	zeros = values.count(0)
	tens = values.count(10)
	rest = len(values)-zeros-tens
	fig = {
	'data': [{	'labels': ['n=0', '0<n<10', 'n>=10'],
				'values': [zeros,rest,tens],
				'type': 'pie'}],
	'layout': {	'title': 'Number of pairs designed for each target region'}
	}

	return fig

def def_fig(data, title, xname, yname):

	data = [go.Histogram(x=data)]
	layout = def_layout(title, xname, yname)
	fig = go.Figure(data=data, layout=layout)

	return fig


def def_layout(title, xname, yname, xfont='Courier New, monospace', yfont='Courier New, monospace', xsize=18, ysize=18, xcolor='#7f7f7f', ycolor='#7f7f7f'):

	layout = go.Layout(	title=title,
					xaxis=dict(	title=xname,
								titlefont=dict(
									family=xfont,
									size=xsize,
									color=xcolor)
								),
					yaxis=dict(	title=yname,
								titlefont=dict(
									family=yfont,
									size=ysize,
									color=ycolor)
								))
	return layout
