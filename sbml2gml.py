#!/usr/bin/python

############################################################################
#
#	Converter SBML to GML
#
#	Created by Paulo Burke at Federal University of Sao Paulo (2015)
#
#	The code is free to anyone use and edit but without any warranty.
#
#	Last update: 02/03/2015
#
############################################################################

import networkx as nx
from libsbml import *
import sys, getopt


def check(value, message=None):
	"""If 'value' is None, prints an error message constructed using
	'message' and then exits with status code 1. If 'value' is an integer,
	it assumes it is a libSBML return status code. If the code value is
	LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
	prints an error message constructed using 'message' along with text from
	libSBML explaining the meaning of the code, and exits with status code 1.
	"""
	if value == None:
		print('LibSBML returned a null value trying to ' + str(message) + '.')
		print('Exiting.')
		sys.exit(1)
	elif type(value) is int:
		if value == LIBSBML_OPERATION_SUCCESS:
			if message:
				print message + '...'
		else:
			print('Error encountered at ' + str(message) + '.')
			print('LibSBML returned error code ' + str(value) + ': "' \
			+ OperationReturnValue_toString(value).strip() + '"')
			print('Exiting.')
			sys.exit(1)
	else:
		if message:
			print message



def getSpeciesList(l):
	l = l.getListOfSpecies()
	sp = []
	for i in range(l.size()):
		sp.append(l[i].getId())
	return sp


def getReactionsList(l):
	l = l.getListOfReactions()
	r = []
	for i in range(l.size()):
		r.append(l[i].getId())
	return r

def main(argv):
#Options treatment ----------------------------------------------

	inputfile = ''
	outputfile = ''
	usestring = 'Usage: ./nxSBML2GML.py.py -i <inputfile(SBML)> -o <outputfile>'

	try:
		opts, args = getopt.getopt(argv,"hi:o:",["input=","output="])
	except getopt.GetoptError:
		print(usestring)
	for opt, arg in opts:
		if opt == '-h':
			print usestring
			sys.exit()
		elif opt in ("-i","--input"):
			inputfile = arg
		elif opt in ("-o","--output"):
			outputfile = arg
	if inputfile=='' or outputfile=='':
		print('Error - Check if all the parameters were set')
		print(usestring)
		sys.exit()

#Reading SBML ----------------------------------------------------

	print('Reading SBML...')
	try:
		reader = SBMLReader()
		document = reader.readSBMLFromFile(inputfile)
		check(document)
	except ValueError:
		print('Error - Could not read SBML Document object')
		sys.exit(1)

	SBMLNetwork = document.getModel()
	check(SBMLNetwork)

	species = getSpeciesList(SBMLNetwork)
	reactions = getReactionsList(SBMLNetwork)

	g = nx.DiGraph()
	molDict = {}
	count = 0
	for sp in species:
		s = SBMLNetwork.getSpecies(sp)
		if s!=None:
			annotation = s.getAnnotationString().replace("<annotation>","").replace("</annotation>","")
			name = s.getName()
			cell_location = s.getCompartment()
			g.add_node(count,name=sp,bdName=name,annotation=annotation, type='molecule',cellLocation=cell_location)
			molDict[sp] = count
			count=count+1

	reactDict = {}
	for react in reactions:
		r = SBMLNetwork.getReaction(react)
		if r!=None:
			annotation = r.getAnnotationString().replace("<annotation>","").replace("</annotation>","")
			name = r.getName()
			notes = r.getNotes()
			if notes != None:
				process = XMLNode.convertXMLNodeToString(notes.getChild(0)).replace('<p xmlns="http://www.w3.org/1999/xhtml">',"").replace('</p>',"").replace('<>',"")
			g.add_node(count,name=react,bdName=name,annotation=annotation, type='reaction',process=process)
			count=count+1
			reactDict[react] = count
			for m in range(r.getNumReactants()):
				mol = r.getReactant(m)
				if mol!=None:
					st = mol.getStoichiometry()
					if st == 0:
						st=-1
						print('weight zero:',mol.getSpecies(),react)
					g.add_edge(molDict[mol.getSpecies()],reactDict[react],weight=st,type='reactant')

			for m in range(r.getNumProducts()):
				mol = r.getProduct(m)
				if mol!=None:
					st = mol.getStoichiometry()
					if st != 0:
						g.add_edge(reactDict[react],molDict[mol.getSpecies()],weight=st,type='product')

			for m in range(r.getNumModifiers()):
				mol = r.getModifier(m)
				if mol!=None:
					g.add_edge(molDict[mol.getSpecies()],reactDict[react],weight=1,type='modifier')

#-------------------------------------------------------------------------------------------------------


	deg = nx.degree(g.to_undirected())
	for d in deg:
		if d[1] == 0:
			g.remove_node(d[0])

	for n in g.nodes(data=True):
		if not len(n[1].keys()):
			g.remove_node(n[0])

	nx.write_gml(g,outputfile)


if __name__ == '__main__':
	main(sys.argv[1:])








