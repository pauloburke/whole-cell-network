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

from igraph import *
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
	usestring = 'Usage: ./igraphSBML.py -i <inputfile(SBML)> -o <outputfile>'

	try:
		opts, args = getopt.getopt(argv,"hi:o:",["input=","output="])
	except getopt.GetoptError:
		print usestring
	for opt, arg in opts:
		if opt == '-h':
			print usestring
			sys.exit()
		elif opt in ("-i","--input"):
			inputfile = arg
		elif opt in ("-o","--output"):
			outputfile = arg
	if inputfile=='' or outputfile=='':
		print 'Error - Check if all the parameters were set'
		print usestring
		sys.exit()

#Reading SBML ----------------------------------------------------

	print 'Reading SBML...'
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

	g = Graph()
	g = g.as_directed()

	g.es["weight"] = 1.0
	g.es["type"] = 'None'

	for sp in species:
		s = SBMLNetwork.getSpecies(sp)
		if s!=None:
			annotation = s.getAnnotationString().replace("<annotation>","").replace("</annotation>","")
			name = s.getName()
			cell_location = s.getCompartment()
			g.add_vertex(sp)
			g.vs(len(g.vs)-1)['bdName'] = name
			g.vs(len(g.vs)-1)['annotation'] = annotation
			g.vs(len(g.vs)-1)['type'] = 0
			g.vs(len(g.vs)-1)['cellLocation'] = cell_location

	for react in reactions:
		r = SBMLNetwork.getReaction(react)
		if r!=None:
			annotation = r.getAnnotationString().replace("<annotation>","").replace("</annotation>","")
			name = r.getName()
			notes = r.getNotes()
			if notes != None:
				process = XMLNode.convertXMLNodeToString(notes.getChild(0)).replace('<p xmlns="http://www.w3.org/1999/xhtml">',"").replace('</p>',"").replace('<>',"")
			g.add_vertex(react)
			g.vs(len(g.vs)-1)['bdName'] = name
			g.vs(len(g.vs)-1)['annotation'] = annotation
			g.vs(len(g.vs)-1)['type'] = 1
			g.vs(len(g.vs)-1)['process'] = process
			for m in range(r.getNumReactants()):
				mol = r.getReactant(m)
				if mol!=None:
					m = g.vs.find(name=mol.getSpecies())
					g.add_edge(m,len(g.vs)-1)
					st = mol.getStoichiometry()
					if st == 0:
						st = -1
					g.es(len(g.es)-1)['weight'] = st
					g.es(len(g.es)-1)['type'] = 'reactant'
			for m in range(r.getNumProducts()):
				mol = r.getProduct(m)
				if mol!=None:
					m = g.vs.find(name=mol.getSpecies())
					g.add_edge(len(g.vs)-1,m)
					g.es(len(g.es)-1)['weight'] = mol.getStoichiometry()
					g.es(len(g.es)-1)['type'] = 'product'
			for m in range(r.getNumModifiers()):
				mol = r.getModifier(m)
				if mol!=None:
					m = g.vs.find(name=mol.getSpecies())
					g.add_edge(m,len(g.vs)-1)
					g.es(len(g.es)-1)['weight'] = 1
					g.es(len(g.es)-1)['type'] = 'modifier'


	g.vs['Degree'] = g.degree()

	Graph.write(g,outputfile)
	




if __name__ == '__main__':
	main(sys.argv[1:])








