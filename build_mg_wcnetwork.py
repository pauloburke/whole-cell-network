#!/usr/bin/python

############################################################################
#
#	Whole-Cell Network
#
#	Created by Paulo Burke at Federal University of Sao Paulo (2015)
#
#	The code is free for anyone to use and edit but without any warranty.
#
############################################################################


import sys, getopt, subprocess, os
import MySQLdb as mdb
from libsbml import *


#Define key proteins WIDs:

RNAPOLIMERASE='RNA_POLYMERASE_HOLOENZYME'
NUSA = 'MG_141_MONOMER'
NUSB = 'MG_027_MONOMER'
NUSG = 'MG_054_MONOMER'
GREA = 'MG_282_MONOMER'
S10 = 'MG_150_MONOMER'
METPEPTDASE = 'MG_172_MONOMER'
DEFORMYLASE = 'MG_106_DIMER'
RIBOSOME = 'RIBOSOME_70S'
RIBOSOMES_SUB = ['RIBOSOME_30S','RIBOSOME_50S']
RIBONUCLEASE = 'MG_104_MONOMER'
DEAD_AMP = 'MG_425_DIMER'
TRNAHYDROLASE = 'MG_083_MONOMER'
RNASEIII = 'MG_367_DIMER'
RNASEJ = 'MG_139_DIMER'
RNASEP = 'MG_0003_465'
RSGA = 'MG_110_MONOMER'
RIBOSOME_30S = 'RIBOSOME_30S'
RIBOSOME_50S = 'RIBOSOME_50S'
IF = ['MG_173_MONOMER','MG_142_MONOMER','MG_196_MONOMER']
EFTU = 'MG_451_DIMER'
EFTS = 'MG_433_DIMER'
EF4 = 'MG_138_MONOMER'
EFG = 'MG_089_DIMER'
EFP = 'MG_026_MONOMER'
RRF = 'MG_435_MONOMER'
RF = 'MG_258_MONOMER'
LEPA = 'MG_138_MONOMER'
TRIGGER = 'MG_238_MONOMER'
RNA_HELICASE = 'MG_308_MONOMER'
GRPE = 'MG_201_DIMER'
DNAJ = 'MG_019_DIMER'
RNADECAY_RNASES = [RNASEIII,RNASEJ]
TRNA_FMET = 'TRNA_FMET'


#Ribosome 30s binding factors
RB30SBF = ['MG_143_MONOMER']
#Ribosome 50s binding factors
RB50SBF = ['MG_384_MONOMER','MG_442_MONOMER','MG_387_MONOMER','MG_335_MONOMER','MG_329_MONOMER']

#Protein Decay
PEPTDASE = ['MG_046_DIMER','MG_324_MONOMER','MG_391_HEXAMER','MG_183_MONOMER','MG_020_MONOMER']
CLP_PROTEASE = 'MG_355_HEXAMER'
FTSH_PROTEASE = 'MG_457_HEXAMER'
LA_PROTEASE = 'MG_239_HEXAMER'
SSRA = 'MG_0004'
SSRA_BP = 'MG_059_MONOMER'

ATP = 'ATP'
TTP = 'TTP'
CTP = 'CTP'
GTP = 'GTP'
URA = 'URA'
PPI = 'PPI'
PI = 'PI'
AMP = 'AMP'
TMP = 'TMP'
CMP = 'CMP'
GMP = 'GMP'
UMP = 'UMP'
ADP = 'ADP'
TDP = 'TDP'
CDP = 'CDP'
GDP = 'GDP'
UDP = 'UDP'
H = 'H'
H2O = 'H2O'
MET = 'MET'
FOR = 'FOR'
DNAA = 'MG_469_7MER_ATP'
DNAB = 'MG_349_HEXAMER'
DNAB2 = 'MG_364_MONOMER'
DATP = 'DATP'
DTTP = 'DTTP'
DCTP = 'DCTP'
DGTP = 'DGTP'
REPLISOME = 'REPLISOME'
SSBINDPROT = 'MG_091_OCTAMER'
GYRASE = 'DNA_GYRASE'
TOPOISOMERASE = 'MG_203_204_TETRAMER'
TOPOISOMERASE1 = 'MG_122_MONOMER'
CHROMO_SEGREGATION_PROTEIN = 'MG_213_214_298_6MER_ADP'
REPLICATION_CR = 'REPLICATION_CHROMOSOME'
DUPLICATED_CR = 'DUPLICATED_CHROMOSOME'
DAMN_CR = 'DAMAGED_CHROMOSOME'
FTSZ_RING = "MG_224_9MER_GDP"

CR = 'CHROMOSOME'

MEMBRANE = 'MEMBRANE'
TERMINAL_ORGANELLE = 'TERMINAL_ORGANELLE'

MEMBRANE_ATTACH_REACT = ['PG160',ATP,H2O]
MEMBRANE_ATTACH_PROD = ['SNGLYP',H,PI,ADP]
MEMBRANE_ATTACH_MOD = ['MG_086_MONOMER','MG_055_170_277_464_476_20MER','MG_210_MONOMER','MG_210_MONOMER','MG_072_DIMER','MG_0001_048']

MEMBRANE_SECR_REACT = [GTP,H2O]
MEMBRANE_SECR_PROD = [H,PI,GDP]
MEMBRANE_SECR_MOD = ['MG_0001_048','MG_297_MONOMER','MG_055_170_277_464_476_20MER']

AMINOACYLATION_TYPE = '<annotation>aminoacylation</annotation>'
RNA_TYPE = '<annotation>RNA</annotation>'
PROTEIN_TYPE = '<annotation>Protein Monomer</annotation>'
TRANSCRIPTION_TYPE = '<annotation>Transcription Reaction</annotation>'
TRANSCRIPTION2_TYPE = '<annotation>preRNA Processing Reaction</annotation>'
TRANSLATION_TYPE = '<annotation>Translation Reaction</annotation>'
TRNAHYDRO_PREFIX = 'PeptidylTRNAHs'
PROTPOLY_PREFIX = 'ProteinPolymerization'
RNAPOLY_PREFIX = 'RnaPolymerization'
RNACLEAVAGE_PREFIX = 'RNACleavage'
RIBOASSEMBLY_PREFIX = 'GTPHs_RibosomeAssembly'
METCLEAV_PREFIX = 'NTerminalMethionineClevage'
DEFORMYLATION_PREFIX = 'Deformylation'
MODIFICATION_PREFIX = 'MODIFICATION'
PROTDECAY_PREFIX = 'ProteinDegradation'
DNABINDREACT_PREFIX = 'DNABoundProteinRelease'
DNAPRIMERPOLY_PREFIX = 'DNAPrimerPolymerization'
CHROMOSEGREGATION_PREFIX = 'GTPHs_ChromosomeSegregation'
FTSZPOLY_PREFIX = 'FtsZ_polymerization'
DNAAPOLY_PREFIX = 'DnaAATP_polmerization'
DNAAPOLY_PREFIX2 = 'DNAPolymerization'
DNADamage_PREFIX = 'DNADamage'
organelProcess = 'Process_TerminalOrganelleAssembly'

complexToDelete = ['DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE',
					'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX','DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE','DNA_POLYMERASE_HOLOENZYME','MG_001_DIMER',
					'MG_469_2MER_1ATP_ADP','MG_469_3MER_2ATP_ADP','MG_469_4MER_3ATP_ADP','MG_469_5MER_4ATP_ADP','MG_469_6MER_5ATP_ADP','MG_469_7MER_6ATP_ADP','RIBOSOME_30S_IF3','RIBOSOME_70S']
					
reactionsToDelete = ['ProLipoproteinDiacylglycerylTransferase','ATPHs_ProteinTranslocation_SecA','ATPHs_ProteinTranslocation_SRP','ATPHs_DNASupercoiling_topoisomeraseIV','ATPHs_DNASupercoiling_gyrase',
						'ATPHs_DNASupercoiling_topoisomeraseI','NTerminalMethionineClevage']

AMINOACIDS = {
'ALA' : 'A',
'ARG' : 'R',
'ASN' : 'N',
'ASP' : 'D',
'CYS' : 'C',
'GLU' : 'E',
'GLN' : 'Q',
'GLY' : 'G',
'HIS' : 'H',
'ILE' : 'I',
'LEU' : 'L',
'LYS' : 'K',
'MET' : 'M',
'PHE' : 'F',
'PRO' : 'P',
'SER' : 'S',
'THR' : 'T',
'TRP' : 'W',
'TYR' : 'Y',
'VAL' : 'V',
'SEC' : 'U',
'PYL' : 'O'
}

AMINOACIDS_INVERSE = dict((v,k) for k,v in AMINOACIDS.iteritems())

deletedReactions = reactionsToDelete
deletedMolecules = complexToDelete

 
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
 

#-------------------------------------------------------------------------------------

def mysqlQuery(c,query):
	try:
		c.execute(query)
		result = c.fetchall()
		return result
	except mdb.Error, e:
		print "MySQL error %d: %s" % (e.args[0],e.args[1])
		sys.exit(1)

#M node creator ---------------------------------------------------------------------

def createMoleculeVertice(net,wid,name,tp,compartment='c'):
	s = net.createSpecies()
	check(s)
	check(s.setId(wid.replace("-","_")))
	check(s.setName(name))
	check(s.setCompartment(compartment))
	check(s.setConstant(False))
	check(s.setInitialAmount(0))
	check(s.setSubstanceUnits('mole'))
	check(s.setBoundaryCondition(False))
	check(s.setHasOnlySubstanceUnits(False))	
	check(s.appendAnnotation(tp))

#Create Reactant ---------------------------------------------------------------------

def createReactant(r,Id,q):
	sp = r.getReactant(Id)
	if sp == None:
		species = r.createReactant()
		check(species)
		check(species.setSpecies(Id.replace("-","_")))
		check(species.setConstant(True))
		check(species.setStoichiometry(q))
	else:
		check(sp.setStoichiometry(sp.getStoichiometry()+q))

#Create Product ---------------------------------------------------------------------

def createProduct(r,Id,q):
	sp = r.getProduct(Id)
	if sp == None:
		species = r.createProduct()
		check(species)
		check(species.setSpecies(Id.replace("-","_")))
		check(species.setConstant(True))
		check(species.setStoichiometry(q))
	else:
		check(sp.setStoichiometry(sp.getStoichiometry()+q))

#Create Modifier ---------------------------------------------------------------------

def createModifier(r,Id):
	sp = r.getModifier(Id)
	if sp == None:
		modifier = r.createModifier()
		check(modifier)
		check(modifier.setSpecies(Id))

#Species List ----------------------------------------------------------------------

def getSpeciesList(l):
	l = l.getListOfSpecies()
	sp = []
	for i in range(l.size()):
		sp.append(l[i].getId())
	return sp


#MAIN ---------------------------------------------------------------------------------
 
def main(argv):

	sqlhost = 'localhost'
	sqluser = ''
	sqlpasswd = ''
	dbname = ''
	outputfile = ''
	usestring = 'Usage: ./build_mg_wcnetwork.py -u <mysqluser> -p <mysqlpswd> -d <dbname> -o <outputfile> [optional: -l <host>, -w (write deleted reactions and molecules to file)]'
	writeDeletedReactions = 0


#Options treatment --------------------------------------------------------------------

	try:
		opts, args = getopt.getopt(argv,"hwu:p:d:o:l:",["user=","passwd=","database=","output=","host"])
	except getopt.GetoptError:
		print usestring
	for opt, arg in opts:
		if opt == '-h':
			print usestring
			sys.exit()
		elif opt == '-w':
			writeDeletedReactions = 1
		elif opt in ("-u","--user"):
			sqluser = arg
		elif opt in ("-p","--passwd"):
			sqlpasswd = arg
		elif opt in ("-d","--database"):
			dbname = arg
		elif opt in ("-o","--output"):
			outputfile = arg
		elif opt in ("-l","--host"):
			sqlhost = arg
	if sqluser=='' or dbname=='' or outputfile=='':
		print 'Error - Check if all the parameters were set'
		print usestring
		sys.exit()

	print 'Whole-Cell Network'
	print 'Created by Paulo Burke at Federal University of Sao Paulo (2015)'
	print 'Initializing...'
	print 'Connecting WholeCellKB database named: ' + dbname + ' at ' +sqlhost+ '...'



#EMBOSS transeq check ----------------------------------

	try:
		with open(os.devnull, 'wb') as devnull:
			process = subprocess.Popen('transeq -h',shell=True,stdin=subprocess.PIPE,stdout=devnull, stderr=subprocess.STDOUT)
	except OSError as e:
		if e.errno == os.errno.ENOENT:
			print 'EMBOSS transeq command not found. Check if is it installed.'
			sys.exit(1)
		else:
			print 'Can\'t run EMBOSS transeq command for not known reason. Check if is it installed.'
			sys.exit(1)

#MySQL connection --------------------------------------
	try:
		con = mdb.connect(host=sqlhost, user=sqluser, passwd=sqlpasswd,db=dbname)
		c = con.cursor(mdb.cursors.DictCursor)
		print "Database connected."
	except mdb.Error, e:
		print "MySQL error %d: %s" % (e.args[0],e.args[1])
		sys.exit(1)



#Creating SBML -----------------------------------------
	try:	
		f = open(outputfile,'w')
	except ValueError:
		print 'Error - Can\'t create a file.'
		sys.exit(1)

	print 'Starting SBML module...'
	try:
		document = SBMLDocument(3, 1)
	except ValueError:
		print('Error - Could not create SBMLDocumention object')
		sys.exit(1)

	network = document.createModel()
	check(network, 'Creating model')


#Reading and creating compartments ---------------------

	counter = 0;
	print 'Reading compartments...'
	compartments = mysqlQuery(c,'select public_entry.wid, public_entry.name from public_compartment inner join public_entry on public_entry.id=public_compartment.parent_ptr_species_component_id')
	for comp in compartments:
		compart = network.createCompartment()
		check(compart)
		check(compart.setId(comp['wid'].replace("-","_")))
		check(compart.setName(comp['name']))
		check(compart.setConstant(True))
		check(compart.setSize(1))
		check(compart.setSpatialDimensions(3))
		check(compart.setUnits('litre'))
		counter = counter+1
	print str(counter)+' compartments created.'


#Retrieving network vertices from database -------------

	print 'Starting molecules reading...'

	queries = [
			['Transcriptional Unit','select public_entry.wid, public_entry.name from public_transcriptionunit inner join public_entry on public_entry.id = public_transcriptionunit.parent_ptr_molecule_id'],
			['Protein Monomer','select public_entry.wid, public_entry.name, public_proteinmonomer.localization_id from public_proteinmonomer inner join public_entry on public_entry.id = public_proteinmonomer.parent_ptr_protein_id'],
			['Protein Complex','select public_entry.wid, public_entry.name from public_proteincomplex inner join public_entry on public_entry.id = public_proteincomplex.parent_ptr_protein_id'],
			['Metabolite','select public_entry.wid, public_entry.name from public_metabolite inner join public_entry on public_entry.id = public_metabolite.parent_ptr_molecule_id']
		]
	for query in queries:
		print 'Reading '+query[0]+'...'
		result = mysqlQuery(c,query[1])
		counter = 0;
		for row in result:
			if query[0] == 'Protein Monomer':
				loc = 'c'
				loc = mysqlQuery(c,"select wid from public_entry where id = "+str(row['localization_id']))
				loc = loc[0]['wid']
				createMoleculeVertice(network,row['wid'].replace("-","_"),row['name'],query[0],loc)
			else:				
				createMoleculeVertice(network,row['wid'].replace("-","_"),row['name'],query[0])
			counter = counter + 1
		print str(counter)+' '+query[0]+' created.'
		
	createMoleculeVertice(network,TERMINAL_ORGANELLE,'Terminal Organelle','Cellular Structure')
	createMoleculeVertice(network,MEMBRANE,'Membrane','Cellular Structure')


#Retrieving network reaction -----------------

	speciesList = getSpeciesList(network)
	modificationReactions = []
	organelAssemblyReactions = []
	print 'Starting general reactions reading...'
	reactions = mysqlQuery(c,'select public_entry.wid, public_entry.name, public_reaction.modification_id, public_reaction.enzyme_id, public_reaction.processes_id, public_reaction.parent_ptr_species_component_id, public_reaction.direction from public_reaction inner join public_entry on public_entry.id = public_reaction.parent_ptr_species_component_id')
	counter = 0
	for row in reactions:
		if row['wid'].replace("-","_") in reactionsToDelete:
			continue
		if row['modification_id'] == None:
			counter = counter+1
			r = network.createReaction()
			check(r)
			check(r.setId(row['wid'].replace("-","_")))
			check(r.setName(row['name']))
			check(r.setReversible(False))
			check(r.setFast(False))
			rev_r = None
			if row['direction'] == 'r':
				rev_r = network.createReaction()
				check(rev_r)
				check(rev_r.setId(row['wid'].replace("-","_").replace("-","_")+"_reverse"))
				check(rev_r.setName(row['name']+" reverse"))
				check(rev_r.setReversible(False))
				check(rev_r.setFast(False))
			pathway = mysqlQuery(c,'select public_entry.wid from public_reaction_pathways inner join public_entry on public_entry.id = public_reaction_pathways.pathway_id where public_reaction_pathways.reaction_id = '+str(row['parent_ptr_species_component_id']))
			if len(pathway):
				pathway = pathway[0]['wid']
				check(r.appendAnnotation(pathway))
				if row['direction'] == 'r':
					check(rev_r.appendAnnotation(pathway))
			if row['processes_id'] != None:
				process = mysqlQuery(c,'select public_entry.wid from public_entry where public_entry.id='+str(row['processes_id']))
				if len(process):
					process = process[0]['wid']
					check(r.setNotes(process.strip('Process_'),True))
					if row['direction'] == 'r':
						check(rev_r.setNotes(process.strip('Process_'),True))
			else:
				check(r.setNotes('Other',True))
				if row['direction'] == 'r':
					check(rev_r.setNotes('Other',True))
			reactants = mysqlQuery(c,'select public_entry.wid,public_entry.name, public_reactionstoichiometryparticipant.coefficient, public_reactionstoichiometryparticipant.compartment_id from public_reaction_stoichiometry inner join public_reactionstoichiometryparticipant on public_reaction_stoichiometry.reactionstoichiometryparticipant_id=public_reactionstoichiometryparticipant.id inner join public_entry on public_entry.id=public_reactionstoichiometryparticipant.molecule_id where reaction_id = '+str(row['parent_ptr_species_component_id']))
			if len(reactants)>0:
				reactantCount = 0
				productCount = 0
				for molecule in reactants:
					compartment = mysqlQuery(c,'select wid from public_entry where id ='+str(molecule['compartment_id']))[0]['wid']
					if compartment == 'e':
						molecule['wid'] = molecule['wid']+'_extracell'
						molecule['name'] = molecule['name']+' Extracellular'
						createModifier(r,MEMBRANE)
						check(r.setNotes('Transmembrane Transport',True))
						if row['direction'] == 'r':
							check(rev_r.setNotes('Transmembrane Transport',True))
					if molecule['wid'].replace("-","_") not in speciesList:
						createMoleculeVertice(network,molecule['wid'].replace("-","_"),molecule['name'],"Metabolite","e")
						speciesList.append(molecule['wid'].replace("-","_"))
					if molecule['coefficient']<0:
						createReactant(r,molecule['wid'],molecule['coefficient']*-1)
						if row['direction'] == 'r':
							createProduct(rev_r,molecule['wid'],molecule['coefficient']*-1)
						reactantCount = reactantCount+1
					else:
						createProduct(r,molecule['wid'],molecule['coefficient'])
						if row['direction'] == 'r':
							createReactant(rev_r,molecule['wid'],molecule['coefficient'])
						productCount = productCount+1
				if reactantCount==0:
					print 'Warning: reaction '+row['wid'].replace("-","_")+' has no reactant!'
				if productCount==0:
					print 'Warning: reaction '+row['wid'].replace("-","_")+' has no product!'
			if row['enzyme_id'] != None:
				enzymes = mysqlQuery(c,'select public_entry.wid from public_enzymeparticipant inner join public_entry on public_entry.id=public_enzymeparticipant.protein_id where public_enzymeparticipant.id ='+str(row['enzyme_id']))
				for enz in enzymes:
					createModifier(r,enz['wid'])
					if row['direction'] == 'r':
						createModifier(rev_r,enz['wid'])
		else:
			#Modification Reaction
			process = mysqlQuery(c,'select wid from public_entry where id = '+str(row['processes_id']))
			process = process[0]['wid']
			if process == organelProcess:
				organelAssemblyReactions.append(row['wid'])
			else:
				modificationReactions.append(row['wid'])
	print str(counter)+' general reactions created.'


#Assign types -----------------------------------------------------------------

	print 'Annotating molecule types...'
	counter = 0
	
	#species
	for i in range(network.getNumSpecies()):
		sp = network.getSpecies(long(i))
		if sp != None:
			spwid = sp.getId()
			typeId = mysqlQuery(c,"select public_speciescomponent_type.type_id from public_entry inner join public_speciescomponent_type on public_speciescomponent_type.speciescomponent_id=public_entry.id where public_entry.wid='"+spwid+"'")
			for tId in typeId:
				typeWid = mysqlQuery(c,'select public_entry.wid from public_entry where public_entry.id='+str(tId['type_id']))
				sp.appendAnnotation(typeWid[0]['wid'].replace("-","_"))
				counter = counter+1

	#reactions
	for i in range(network.getNumReactions()):
		r = network.getReaction(long(i))
		if r != None:
			if r.getAnnotationString().replace("<annotation>","").replace("</annotation>","") != '':
				rwid = r.getId()
				typeId = mysqlQuery(c,"select public_speciescomponent_type.type_id from public_entry inner join public_speciescomponent_type on public_speciescomponent_type.speciescomponent_id=public_entry.id where public_entry.wid='"+rwid+"'")
				for tId in typeId:
					typeWid = mysqlQuery(c,'select public_entry.wid from public_entry where public_entry.id='+str(tId['type_id']))
					r.appendAnnotation(typeWid[0]['wid'].replace("-","_"))
					counter = counter+1
	print str(counter)+' types annotated.'



#Remove ribossomal assembly reactions -----------------------------------------
	
	print 'Removing ribossomal assembly old reactions...'
	counter = 0
	maxi = network.getNumReactions()
	i=0
	while i<maxi:
		if i<network.getNumReactions():
			r = network.getReaction(long(i))
			if r != None:
				if RIBOASSEMBLY_PREFIX in r.getId():
					deletedReactions.append(r.getId())
					network.removeReaction(r.getId())
					counter = counter+1
					i = i-1
		i = i+1
	print str(counter)+' ribossomal assembly old reactions removed.'


#DNA Replication --------------------------------------------------------------

	genome = mysqlQuery(c,'select sequence from public_chromosome limit 1')
	genome = genome[0]['sequence']

	print 'Creating DNA replication and cytokinesis reactions...'
	counter = 0

	#create Chromosome molecule
	createMoleculeVertice(network,REPLICATION_CR,'Chromosome ready to replication','Cellular Structure')

	#Replication iniciation reaction
	counter += 1
	r = network.createReaction()
	check(r)
	check(r.setId('REPLICATION_INITIATION'))
	check(r.setName('Replication Initiation'))
	check(r.setReversible(False))
	check(r.setFast(False))
	check(r.appendAnnotation('Replication Reaction'))
	check(r.setNotes('Replication',True))
	#add Transcriptional Units	
	tUnits = mysqlQuery(c,'select public_entry.wid from public_transcriptionunit inner join public_entry on public_entry.id = public_transcriptionunit.parent_ptr_molecule_id')
	for tu in tUnits:
		createReactant(r,tu['wid'].replace("-","_"),1)
	#add DnaA
	createModifier(r,DNAA)
	createModifier(r,DNAB)
	createModifier(r,DNAB2)
	
	#Add chromosome
	createProduct(r,REPLICATION_CR,1)

	#create Replication Reaction
	counter += 1
	r = network.createReaction()
	check(r)
	check(r.setId('REPLICATION_ELONGATION'))
	check(r.setName('Replication Elongation'))
	check(r.setReversible(False))
	check(r.setFast(False))
	check(r.appendAnnotation('Replication Reaction'))
	check(r.setNotes('Replication',True))

	#Add chromosome
	createReactant(r,REPLICATION_CR,1)
	#Add supercoiling proteins
	createModifier(r,GYRASE)
	createModifier(r,TOPOISOMERASE)
	createModifier(r,TOPOISOMERASE1)

	#Add replisome
	createModifier(r,REPLISOME)
	#Add dNTP
	ACounter = 0
	TCounter = 0
	CCounter = 0
	GCounter = 0
	for n in genome:
		if n == 'A':
			ACounter = ACounter+1
		if n == 'T':
			TCounter = TCounter+1
		if n == 'C':
			CCounter = CCounter+1
		if n == 'G':
			GCounter = GCounter+1
	createReactant(r,DATP,2*TCounter)
	createReactant(r,DTTP,2*ACounter)
	createReactant(r,DCTP,2*GCounter)
	createReactant(r,DGTP,2*CCounter)
	createReactant(r,ATP,TCounter)
	createReactant(r,URA,ACounter)
	createReactant(r,CTP,GCounter)
	createReactant(r,GTP,CCounter)
	#Add ss Binding protein
	createModifier(r,SSBINDPROT)

	createProduct(r,PPI,(ACounter+TCounter+CCounter+GCounter)*3)
	createProduct(r,AMP,TCounter)
	createProduct(r,UMP,ACounter)
	createProduct(r,CMP,GCounter)
	createProduct(r,GMP,CCounter)

	#DNA repair
	#Reading repair reactions
	repReacts = mysqlQuery(c,"select public_reaction.parent_ptr_species_component_id from public_reaction inner join public_entry on public_entry.id=public_reaction.processes_id where public_entry.wid='Process_DNARepair'")
	for rep in repReacts:
		reaction = mysqlQuery(c,'select public_entry.wid from public_reaction inner join public_entry on public_entry.id=public_reaction.parent_ptr_species_component_id where public_reaction.parent_ptr_species_component_id='+str(rep['parent_ptr_species_component_id']))
		reaction = reaction[0]['wid'].replace("-","_")
		repR = network.getReaction(reaction)
		if repR!=None:
			modN = repR.getNumModifiers()
			for i in range(modN):
				mod = repR.getModifier(long(i))
				if mod!=None:
					createModifier(r,mod.getSpecies())
			network.removeReaction(reaction)

	#Create new DNAs
	createMoleculeVertice(network,DUPLICATED_CR,'Duplicated Chromosome','Cellular Structure')
	createProduct(r,DUPLICATED_CR,1)



	#Create Segregation reaction
	counter += 1
	r = network.createReaction()
	check(r)
	check(r.setId('CELLULAR_DIVISION'))
	check(r.setName('Cellular Division'))
	check(r.setReversible(False))
	check(r.setFast(False))
	check(r.appendAnnotation('Replication Reaction'))
	check(r.setNotes('Replication',True))

	#Add chromosome
	createReactant(r,DUPLICATED_CR,1)

	segReacts = mysqlQuery(c,"select public_reaction.parent_ptr_species_component_id from public_reaction inner join public_entry on public_entry.id=public_reaction.processes_id where public_entry.wid='Process_ChromosomeSegregation'")
	for seg in segReacts:
		reaction = mysqlQuery(c,'select public_entry.wid from public_reaction inner join public_entry on public_entry.id=public_reaction.parent_ptr_species_component_id where public_reaction.parent_ptr_species_component_id='+str(seg['parent_ptr_species_component_id']))
		reaction = reaction[0]['wid'].replace("-","_")
		segR = network.getReaction(reaction)
		if segR !=None:
			modN = segR.getNumModifiers()
			for i in range(modN):
				mod = segR.getModifier(long(i))
				if mod!=None:
					createModifier(r,mod.getSpecies())
					createReactant(r,GTP,1)
					createProduct(r,GDP,1)
					createProduct(r,PI,1)
			network.removeReaction(reaction)
	createModifier(r,CHROMO_SEGREGATION_PROTEIN)
	createModifier(r,MEMBRANE)
	createModifier(r,TERMINAL_ORGANELLE)
	createModifier(r,FTSZ_RING)

		
	print counter,'DNA replication and cytokinesis reactions created.'


#Removing DNA Polymerization Old Reactions ------------------------------------

	print 'Removing DNA Replication old reactions...'
	counter = 0
	maxi = network.getNumReactions()
	i=0
	while i<maxi:
		if i<network.getNumReactions():
			r = network.getReaction(long(i))
			if r != None:
				if DNABINDREACT_PREFIX in r.getId() or DNAPRIMERPOLY_PREFIX in r.getId() or DNAAPOLY_PREFIX2 in r.getId() or CHROMOSEGREGATION_PREFIX in r.getId() or FTSZPOLY_PREFIX in r.getId() or DNAAPOLY_PREFIX in r.getId() or DNADamage_PREFIX in r.getId():
					deletedReactions.append(r.getId())
					network.removeReaction(r.getId())
					counter = counter+1
					i = i-1
		i = i+1
	reacts = mysqlQuery(c,"select public_reaction.parent_ptr_species_component_id from public_reaction inner join public_entry on public_entry.id=public_reaction.processes_id where public_entry.wid='Process_Replication' or public_entry.wid='Process_ChromosomeSegregation'")
	for reac in reacts:
		reaction = mysqlQuery(c,'select public_entry.wid from public_reaction inner join public_entry on public_entry.id=public_reaction.parent_ptr_species_component_id where public_reaction.parent_ptr_species_component_id='+str(rep['parent_ptr_species_component_id']))
		reaction = reaction[0]['wid'].replace("-","_")
		segR = network.getReaction(reaction)
		if segR !=None:
			deletedReactions.append(reaction)
			network.removeReaction(reaction)
			counter = counter+1
	print str(counter)+' DNA Replication old reactions removed.'
		


#Transcription reactions ------------------------------------------------------
	polycistrons = {}
	counter = 0
	decayCounter = 0
	processingCounter = 0
	print 'Creating DNA transcription reactions and RNA decay reactions...'

	modificationMolecules = mysqlQuery(c,'select public_entry.wid from public_modificationreaction inner join public_entry on public_entry.id=public_modificationreaction.molecule_id')
	modificationMolecules = [mol['wid'] for mol in modificationMolecules]
	modificationMolecules = [i for i in set(modificationMolecules) if not ('_MONOMER' in i)]
	#print modificationMolecules

	transUnits = mysqlQuery(c,'select public_entry.wid, public_entry.name, public_transcriptionunit.parent_ptr_molecule_id from public_transcriptionunit inner join public_entry on public_entry.id = public_transcriptionunit.parent_ptr_molecule_id')
	for tr in transUnits:
		counter = counter+1
		r = network.createReaction()
		check(r)
		check(r.setId(tr['wid'].replace("-","_")+'_TRANSCRIPTION'))
		check(r.setName(tr['name']+' Transcription Reaction'))
		check(r.setReversible(False))
		check(r.setFast(False))
		check(r.appendAnnotation('Transcription Reaction'))
		check(r.setNotes('Transcription',True))
		#ADD RNA Polimerase
		createModifier(r,RNAPOLIMERASE)
		#ADD NusA, B, G, GreA and S10
		createModifier(r,NUSA)
		createModifier(r,NUSB)
		createModifier(r,NUSG)
		createModifier(r,S10)
		createModifier(r,GREA)
		createModifier(r,TOPOISOMERASE1)
		#ADD transcription unit
		createModifier(r,tr['wid'])
		#add transcription factors
		tFactors = mysqlQuery(c,'select public_entry.wid from public_transcriptionalregulation inner join public_entry on public_entry.id=public_transcriptionalregulation.transcription_factor_id where public_transcriptionalregulation.transcription_unit_id='+str(tr['parent_ptr_molecule_id']))	
		for tf in tFactors:
			createModifier(r,tf['wid'])
		#add genes (mRNAs)
		genes = mysqlQuery(c,'select public_entry.wid,public_gene.coordinate,public_gene.length,public_gene.direction from public_transcriptionunit_genes inner join public_entry on public_entry.id=public_transcriptionunit_genes.gene_id inner join public_gene on public_gene.parent_ptr_molecule_id=public_transcriptionunit_genes.gene_id where public_transcriptionunit_genes.transcriptionunit_id='+str(tr['parent_ptr_molecule_id']))

		if len(genes)>0:
			mRNA = ''
			if len(genes)>1:
				namesGenes = [g['wid'] for g in genes]
				processingFlag = len(set(namesGenes).intersection(modificationMolecules))
				if processingFlag:
					processingCounter = processingCounter+1
					if (tr['wid'].replace("-","_")+'_PRE_RNA') not in speciesList:
						createMoleculeVertice(network,tr['wid'].replace("-","_")+'_PRE_RNA',tr['wid'].replace("-","_")+' Pre RNA',"preRNA")
						speciesList.append(tr['wid'].replace("-","_")+'_PRE_RNA')
					createProduct(r,tr['wid']+'_PRE_RNA',1)
					#Create Pre RNA processing reaction
					pr = network.createReaction()
					check(pr)
					check(pr.setId(tr['wid'].replace("-","_")+'_RNA_PROCESSING'))
					check(pr.setName(tr['name']+' RNA Processing'))
					check(pr.setReversible(False))
					check(pr.setFast(False))
					check(pr.appendAnnotation('preRNA Processing Reaction'))
					check(pr.setNotes('Transcription',True))
					#Add RNASEIII
					createModifier(pr,RNASEIII)
					#ADD H2O	
					createReactant(pr,H2O,len(genes))
					#ADD PRE RNA
					createReactant(pr,tr['wid']+'_PRE_RNA',1)
					mRNA = tr['wid'].replace("-","_")+'_PRE_RNA'
				else:
					mRNA = tr['wid'].replace("-","_")+'_POLYMRNA'
					if mRNA not in speciesList:
						createMoleculeVertice(network,mRNA,tr['wid'].replace("-","_")+' Poly mRNA',"Poly mRNA")
						speciesList.append(mRNA)
			else:
				mRNA = genes[0]['wid']
				if mRNA not in speciesList:
					if '_' in mRNA:
						createMoleculeVertice(network,mRNA,mRNA,'mRNA')
					else:
						createMoleculeVertice(network,mRNA,mRNA,'tRNA')		
					speciesList.append(mRNA)

			#Create mRNA
			createProduct(r,mRNA,1)
			if not processingFlag:
				polycistrons[mRNA] = []	
			ACounter=0
			UCounter=0
			CCounter=0
			GCounter=0
			for gene in genes:
				if processingFlag:
					if gene['wid'] not in speciesList:
						if '_' in gene['wid']:
							createMoleculeVertice(network,gene['wid'],gene['wid'],'mRNA')
						elif 'MGrrnA' in gene['wid']:
							createMoleculeVertice(network,gene['wid'],gene['wid'],'rRNA')
						else:
							createMoleculeVertice(network,gene['wid'],gene['wid'],'tRNA')
						speciesList.append(gene['wid'])
					createProduct(pr,gene['wid'],1)
					polycistrons[gene['wid']] = []
					polycistrons[gene['wid']].append(gene['wid'])
					
				else:
					if not any(gene['wid'] in mm for mm in modificationMolecules):	
						polycistrons[mRNA].append(gene['wid'])
				ini = gene['coordinate']-1
				direction = gene['direction']
				end = 0;
				if direction == 'f':
					end = ini+gene['length']
					for i in range(ini,end):
						if genome[i] == "A":
							UCounter = UCounter+1
						elif genome[i] == "T":
							ACounter = ACounter+1
						elif genome[i] == "C":
							GCounter = GCounter+1
						elif genome[i] == "G":
							CCounter = CCounter+1
				elif direction == 'r':
					end = ini-gene['length']+1
					for i in range(end,ini+1):
						if genome[i] == "A":
							UCounter = ACounter+1
						elif genome[i] == "T":
							ACounter = TCounter+1
						elif genome[i] == "C":
							GCounter = CCounter+1
						elif genome[i] == "G":
							CCounter = GCounter+1
				
			#Creating RNA decay reaction
			for rnase in RNADECAY_RNASES:
				decayCounter = decayCounter+1
				dr = network.createReaction()
				check(dr)
				check(dr.setId(mRNA+'_'+rnase+'_DECAY'))
				check(dr.setName(mRNA+' '+rnase+' Decay Reaction'))
				check(dr.setReversible(False))
				check(dr.setFast(False))
				check(dr.appendAnnotation('RNA Decay Reaction'))
				check(dr.setNotes('RNA Decay',True))
				#ADD Ribonuclease
				createModifier(dr,RIBONUCLEASE)
				#ADD RNASE
				createModifier(dr,rnase)
				#ADD RNA Helicase
				createModifier(dr,DEAD_AMP)
				#ADD RNA
				createReactant(dr,mRNA,1)
				#Add AMPS
				createProduct(dr,AMP,ACounter)
				#Add MPURACILS
				createProduct(dr,UMP,UCounter)
				#Add GMPS
				createProduct(dr,GMP,GCounter)
				#Add CMPS
				createProduct(dr,CMP,CCounter)
				#Add H
				createProduct(dr,H,ACounter+UCounter+CCounter+GCounter)
				# ADD DEAD ATPs
				createReactant(dr,ATP,ACounter+UCounter+CCounter+GCounter)
				createProduct(dr,ADP,ACounter+UCounter+CCounter+GCounter)
				createProduct(dr,PI,ACounter+UCounter+CCounter+GCounter)
				#ADD H2O
				createReactant(dr,H2O,(ACounter+UCounter+CCounter+GCounter)/2)

							
			#Add ATPS
			createReactant(r,ATP,ACounter)
			#Add URACILS
			createReactant(r,URA,UCounter)
			#Add GTPS
			createReactant(r,GTP,GCounter)
			#Add CTPS
			createReactant(r,CTP,CCounter)
			#Add PPis
			createProduct(r,PPI,ACounter+UCounter+CCounter+GCounter)			
			
		else:
			print 'Warning: reaction '+tr['wid'].replace("-","_")+'_TRANSCRIPTION'+' has no product!'
	print str(counter)+' transcription reactions created.'
	print str(processingCounter)+' pre-RNA processing reactions created.'
	print str(decayCounter)+' RNA decay reactions created.'
		

#Removing RNA polymerization old reactions ------------------------

	print 'Removing RNA Polymerization old reactions...'
	counter = 0
	maxi = network.getNumReactions()
	i=0
	while i<maxi:
		if i<network.getNumReactions():
			r = network.getReaction(long(i))
			if r != None:
				if RNAPOLY_PREFIX in r.getId():
					deletedReactions.append(r.getId())
					network.removeReaction(r.getId())
					counter = counter+1
					i = i-1
		i = i+1
	print str(counter)+' RNA polymerization old reactions removed.'

#Removing RNA cleavage old reactions ------------------------

	print 'Removing RNA Cleavage old reactions...'
	counter = 0
	maxi = network.getNumReactions()
	i=0
	while i<maxi:
		if i<network.getNumReactions():
			r = network.getReaction(long(i))
			if r != None:
				if RNACLEAVAGE_PREFIX in r.getId():
					deletedReactions.append(r.getId())
					network.removeReaction(r.getId())
					counter = counter+1
					i = i-1
		i = i+1
	print str(counter)+' RNA cleavage old reactions removed.'


#EF Reactions ---------------------------------------------------------------

	print 'Creating elongation factor activation and deactivation reactions...'
	EFs = [EFTU,EFG,EFP,RRF]
	counter = 0

	for EF in EFs:
		#EF Activation
		counter += 1
		r = network.createReaction()
		check(r)
		check(r.setId(EF+'_ACTIVATION'))
		check(r.setName(EF+' Activation'))
		check(r.setReversible(False))
		check(r.setFast(False))
		check(r.appendAnnotation('Elongation Factor Modification'))
		check(r.setNotes('Protein Modification',True))
		createReactant(r,EF,1)
		createReactant(r,GTP,1)
		if EF+'_'+GTP not in speciesList:
			createMoleculeVertice(network,EF+'_'+GTP,EF+' activated with '+GTP,'Elongation Factor')
			speciesList.append(EF+'_'+GTP)
		createProduct(r,EF+'_'+GTP,1)

		#EF GDP Dissociation
		counter += 1
		r = network.createReaction()
		check(r)
		check(r.setId(EF+'_DISSOCIATION'))
		check(r.setName(EF+' GDP Dissociation'))
		check(r.setReversible(False))
		check(r.setFast(False))
		check(r.appendAnnotation('Elongation Factor Modification'))
		check(r.setNotes('Transcription',True))
		createProduct(r,EF,1)
		createProduct(r,GDP,1)
		createModifier(r,EFTS)
		if EF+'_'+GDP not in speciesList:
			createMoleculeVertice(network,EF+'_'+GDP,EF+' inactive with '+GDP,'Elongation Factor')
			speciesList.append(EF+'_'+GDP)
		createReactant(r,EF+'_'+GDP,1)

	print counter,'elongation factor activation and deactivation reactions created.'

#Modification Reactions -------------------------------------------------------

	tRNAdict = {}

	print 'Reading RNA modification reactions...'
	aminoacylationCounter = 0
	tRNACounter = 0
	rnaHydrolysisCounter = 0
	maturationCounter = 0
	for reaction in modificationReactions:
		react = mysqlQuery(c,"select public_reaction.modification_id, public_reaction.enzyme_id, public_reaction.parent_ptr_species_component_id from public_entry inner join public_reaction on public_reaction.parent_ptr_species_component_id=public_entry.id where public_entry.wid='"+reaction+"'")
		react = react[0]
		modification = mysqlQuery(c,'select public_entry.wid from public_modificationreaction inner join public_entry on public_entry.id=public_modificationreaction.molecule_id where public_modificationreaction.id='+str(react['modification_id']))
		modification = modification[0]['wid']
		reactants = mysqlQuery(c,'select public_entry.wid,public_entry.name, public_reactionstoichiometryparticipant.coefficient from public_reaction_stoichiometry inner join public_reactionstoichiometryparticipant on public_reaction_stoichiometry.reactionstoichiometryparticipant_id=public_reactionstoichiometryparticipant.id inner join public_entry on public_entry.id=public_reactionstoichiometryparticipant.molecule_id where reaction_id = '+str(react['parent_ptr_species_component_id']))
		#Aminoacylation Reactions creation
		if 'Aminoacylation' in reaction:
			aminoacylationCounter = aminoacylationCounter+1
			aa = ''
			if len(reactants)>0:
				if modification == 'MG488':
					aa = 'FMET'
				else:
					for molecule in reactants:
						if AMINOACIDS.has_key(molecule['wid']):
							aa = molecule['wid']
							break
				
				if 'TRNA_'+aa not in speciesList:
					createMoleculeVertice(network,'TRNA_'+aa,'tRNA with '+aa,'Aminoacylated tRNA')
					speciesList.append('TRNA_'+aa)
					tRNACounter = tRNACounter+1
				if not tRNAdict.has_key(aa):
					tRNAdict[aa] = []
				tRNAdict[aa].append(modification)
				r = network.createReaction()
				check(r)
				check(r.setId(modification+'_AMINOACYLATION'))
				check(r.setName(modification+' Aminoacylation Reaction'))
				check(r.setReversible(False))
				check(r.setFast(False))
				check(r.appendAnnotation('Aminoacylation Reaction'))
				check(r.setNotes('Aminoacylation',True))
				#ADD tRNA
				createReactant(r,modification,1)
				#ADD AA
				createReactant(r,aa,1)
				#ADD ATP
				createReactant(r,ATP,1)
				#ADD tRNA with AA
				createProduct(r,'TRNA_'+aa,1)
				#ADD PPI
				createProduct(r,PPI,1)
				#ADD AMP
				createProduct(r,AMP,1)
				#ADD Enzymes
				if react['enzyme_id'] != None:
					enzymes = mysqlQuery(c,'select public_entry.wid from public_enzymeparticipant inner join public_entry on public_entry.id=public_enzymeparticipant.protein_id where public_enzymeparticipant.id ='+str(react['enzyme_id']))
					for enz in enzymes:
						createModifier(r,enz['wid'])
				#Aminoacylated tRNA EFTU addition
				if aa != 'FMET':
					checkReaction = network.getReaction('TRNA_'+aa+'_TU_ADDITION')
					if checkReaction == None:
						r = network.createReaction()
						check(r)
						check(r.setId('TRNA_'+aa+'_TU_ADDITION'))
						check(r.setName('TRNA_'+aa+' EFTU Addition'))
						check(r.setReversible(False))
						check(r.setFast(False))
						check(r.appendAnnotation('EFTU TRNA Complex Formation'))
						check(r.setNotes('Aminoacylation',True))
						createReactant(r,'TRNA_'+aa,1)
						createReactant(r,EFTU+'_'+GTP,2)
						if 'TRNA_'+aa+'_TU' not in speciesList:
							createMoleculeVertice(network,'TRNA_'+aa+'_TU','tRNA with '+aa+' and EFTU','Aminoacylated tRNA')
							speciesList.append('TRNA_'+aa+'_TU')
						createProduct(r,'TRNA_'+aa+'_TU',1)
				
		#RNAs MATURATION
		elif network.getSpecies(modification).getAnnotationString() == RNA_TYPE:
			maturationCounter = maturationCounter+1
			bases = [AMP,UMP,CMP,GMP]
			r = network.getReaction('PRE_'+modification+'_MATURATION')
			if r == None:
				r = network.createReaction()
				check(r)
				check(r.setId('PRE_'+modification+'_MATURATION'))
				check(r.setName('PRE_'+modification+' Maturation Reaction'))
				check(r.setReversible(False))
				check(r.setFast(False))
				check(r.appendAnnotation('RNA Maturation Reaction'))
				check(r.setNotes('Transcription',True))
				if len(reactants)>0:
					for molecule in reactants:
						match = 0
						for base in bases:
							if base in molecule['wid']:
								match = 1
								break
						if not match:
							if molecule['coefficient']<0:
								createReactant(r,molecule['wid'],molecule['coefficient']*-1)
							else:
								createProduct(r,molecule['wid'],molecule['coefficient'])
					#ADD RIBONUCLEASE P
					createModifier(r,RNASEP)
					createReactant(r,H2O,1)
					createProduct(r,H,1)
					#Create product and pre product				
					createProduct(r,modification,1)
					if 'PRE_'+modification not in speciesList:
						createMoleculeVertice(network,'PRE_'+modification,'PRE_'+modification,'preRNA')
						speciesList.append('PRE_'+modification)
					createReactant(r,'PRE_'+modification,1)
					#Correct transcription reaction
					reactionFound = 0
					for i in range(network.getNumReactions()):
						ptrnar = network.getReaction(long(i))
						if ptrnar != None:
							if ptrnar.getAnnotationString() == TRANSCRIPTION_TYPE or ptrnar.getAnnotationString() == TRANSCRIPTION2_TYPE:
								for j in range(ptrnar.getNumProducts()):
									gene = ptrnar.getProduct(long(j))
									if gene != None:
										if gene.getSpecies() == modification:
											ptrnar.removeProduct(modification)
											createProduct(ptrnar,'PRE_'+modification,1)
			else:
				if len(reactants)>0:
					for molecule in reactants:
						match = 0
						for base in bases:
							if base in molecule['wid']:
								match = 1
								break
						if not match:
							if molecule['coefficient']<0:
								createReactant(r,molecule['wid'],molecule['coefficient']*-1)
							else:
								createProduct(r,molecule['wid'],molecule['coefficient'])
			#ADD Enzymes
			if react['enzyme_id'] != None:
				enzymes = mysqlQuery(c,'select public_entry.wid from public_enzymeparticipant inner join public_entry on public_entry.id=public_enzymeparticipant.protein_id where public_enzymeparticipant.id ='+str(react['enzyme_id']))
				for enz in enzymes:
					if r.getModifier(enz['wid'])==None:
						createModifier(r,enz['wid'])
			

			
	#tRNA Hydrolysis reaction
	for aa in tRNAdict:
		rnaHydrolysisCounter = rnaHydrolysisCounter+1
		hr = network.createReaction()
		check(hr)
		check(hr.setId('TRNA_'+aa+'_HYDROLYSIS'))
		check(hr.setName('TRNA_'+aa+' Hydrolysis Reaction'))
		check(hr.setReversible(False))
		check(hr.setFast(False))
		check(hr.appendAnnotation('tRNA Hydrolysis Reaction'))
		check(hr.setNotes('Aminoacylation',True))
		#ADD Ribonuclease
		createModifier(hr,TRNAHYDROLASE)
		#ADD TRNA with AA
		createReactant(hr,'TRNA_'+aa,1)
		#ADD H2O
		createReactant(hr,H2O,1)
		#Add AA
		createProduct(hr,aa,1)
		#ADD tRNA
		for trna in tRNAdict[aa]:		
			createProduct(hr,trna,1.0/len(tRNAdict[aa]))


	print str(rnaHydrolysisCounter)+' Aminoacylation reactions created.'
	print str(tRNACounter)+' tRNA with Aminoacid created.'
	print str(aminoacylationCounter)+' tRNA Hydrolysis reactions created.'
	print str(maturationCounter)+' RNA Maturation reactions created.'



#Removing aminoacylated tRNAs hydrolysis old reactions ------------------------

	print 'Removing aminoacylated tRNAs hydrolysis old reactions...'
	counter = 0
	maxi = network.getNumReactions()
	i=0
	while i<maxi:
		if i<network.getNumReactions():
			r = network.getReaction(long(i))
			if r != None:
				if TRNAHYDRO_PREFIX in r.getId():
					deletedReactions.append(r.getId())
					network.removeReaction(r.getId())
					counter = counter+1
					i = i-1
		i = i+1
	print str(counter)+' aminoacylated tRNAs hydrolysis old reactions removed.'



#Translation Reactions ----------------------------------------------------------

	counter = 0
	ICcounter = 0
	speciesCounter = 0
	decayCounter = 0
	print 'Reading Protein synthesis reactions...'

	

	aaCounter = dict((v,k) for k,v in AMINOACIDS.iteritems())
	for key in aaCounter:
		aaCounter[key] = 0

	for mRNA in polycistrons:
		proteinCount=0
		if len(polycistrons[mRNA]):
			gene = mysqlQuery(c,'select public_entry.id,public_gene.coordinate,public_gene.length,public_gene.direction from public_entry inner join public_gene on public_gene.parent_ptr_molecule_id=public_entry.id where public_entry.wid="'+polycistrons[mRNA][0]+'"')
			gene = gene[0]
			proteinCount = mysqlQuery(c,'select count(*) as c from public_proteinmonomer where public_proteinmonomer.gene_id='+str(gene['id']))
			proteinCount = int(proteinCount[0]['c'])
		if proteinCount:
			N = len(polycistrons[mRNA])
			#Initiation Complex Formation
			ICcounter +=1
			ir = network.createReaction()
			check(ir)
			check(ir.setId(mRNA+'_TRANSLATION_IC_FORMATION'))
			check(ir.setName(mRNA+' Translation Initiation Complex Formation'))
			check(ir.setReversible(False))
			check(ir.setFast(False))
			check(ir.appendAnnotation('Translation Initiation Reaction'))
			check(ir.setNotes('Translation',True))

			createModifier(ir,RIBOSOME_30S)
			createModifier(ir,IF[0])
			createModifier(ir,IF[1])
			createModifier(ir,IF[2])
			createModifier(ir,RNA_HELICASE)
			createReactant(ir,EFP+'_'+GTP,N)
			createReactant(ir,TRNA_FMET,N)
			createReactant(ir,GTP,N) #From IF2
			#Add mRNA
			createModifier(ir,mRNA.replace("-","_"))
			if mRNA.replace("-","_")+'_TRANSLATION_IC' not in speciesList:
				createMoleculeVertice(network,mRNA.replace("-","_")+'_TRANSLATION_IC',mRNA.replace("-","_")+' Translation Initiation Complex','Translation Initiation Complex')
				speciesList.append(mRNA.replace("-","_")+'_TRANSLATION_IC')
			createProduct(ir,mRNA.replace("-","_")+'_TRANSLATION_IC',1)



			#Translation Reaction
			r = network.createReaction()
			check(r)
			check(r.setId(mRNA.replace("-","_")+'_TRANSLATION'))
			check(r.setName(mRNA+' Translation Reaction'))
			check(r.setReversible(False))
			check(r.setFast(False))
			check(r.appendAnnotation('Translation Reaction'))
			check(r.setNotes('Translation',True))
			#ADD Iniciation complex
			createReactant(r,mRNA.replace("-","_")+'_TRANSLATION_IC',1)
			
			createModifier(r,RIBOSOME_30S)
			createModifier(r,RIBOSOME_50S)
			createModifier(r,METPEPTDASE)
			createModifier(r,DEFORMYLASE)
			#ADD mRNA
			createModifier(r,mRNA.replace("-","_"))
			#ADD EF products
			createModifier(r,LEPA)
			
			for geneName in polycistrons[mRNA]:
				gene = mysqlQuery(c,'select public_entry.id,public_gene.coordinate,public_gene.length,public_gene.direction from public_entry inner join public_gene on public_gene.parent_ptr_molecule_id=public_entry.id where public_entry.wid="'+geneName+'"')
				gene = gene[0]
				protein = mysqlQuery(c,'select public_entry.wid,public_entry.name,public_proteinmonomer.parent_ptr_protein_id from public_proteinmonomer inner join public_entry on  public_entry.id=public_proteinmonomer.parent_ptr_protein_id where public_proteinmonomer.gene_id='+str(gene['id']))

				protein = protein[0]
				counter = counter+1		

				cell_local = network.getSpecies(protein['wid']).getCompartment()

				#Membrane protein
				if cell_local == 'm':
					for m in MEMBRANE_ATTACH_REACT:
						createReactant(r,m,1)
					for m in MEMBRANE_ATTACH_PROD:
						createProduct(r,m,1)
					for m in MEMBRANE_ATTACH_MOD:
						createModifier(r,m)
				
				#Secreted
				if cell_local == 'e':
					for m in MEMBRANE_SECR_REACT:
						createReactant(r,m,1)
					for m in MEMBRANE_SECR_PROD:
						createProduct(r,m,1)
					for m in MEMBRANE_SECR_MOD:
						createModifier(r,m)

				#ADD Release factors
				createModifier(r,RF)
				createReactant(r,RRF+'_'+GTP,1)
				createProduct(r,RRF+'_'+GDP,1)
				
				#Add EFP
				createProduct(r,EFP+'_'+GDP,N)
				createProduct(r,GDP,N)
				createProduct(r,PI,N)
				createProduct(r,EFP+'_'+GDP,1)
				createProduct(r,GDP,2)
				createProduct(r,PI,2)			
				
				#ADD formylase and methionine aminopeptidase sub-products

				createProduct(r,MET,1)
				createReactant(r,H2O,2)
				createProduct(r,FOR,1)
				createProduct(r,H,1)
		
				#Add prosthetic groups
				pGroups = mysqlQuery(c,'select public_entry.wid, public_prostheticgroupparticipant.coefficient from public_protein_prosthetic_groups inner join public_prostheticgroupparticipant on public_prostheticgroupparticipant.id=public_protein_prosthetic_groups.prostheticgroupparticipant_id inner join public_entry on public_entry.id=public_prostheticgroupparticipant.metabolite_id where public_protein_prosthetic_groups.protein_id='+str(protein['parent_ptr_protein_id']))
				for prost in pGroups:
					if prost['coefficient'] != None:
						createReactant(r,prost['wid'],prost['coefficient'])
				#Add chaperones
				chaperones = mysqlQuery(c,'select public_entry.wid from public_protein_chaperones inner join public_entry on public_entry.id=public_protein_chaperones.to_protein_id where public_protein_chaperones.from_protein_id='+str(protein['parent_ptr_protein_id']))
				for chap in chaperones:
					if chap != None:
						createModifier(r,chap['wid'])
				if len(chaperones)>0 and cell_local != 'e':
					createModifier(r,TRIGGER)
				if len(chaperones)>0:
					createModifier(r,GRPE)
					createModifier(r,DNAJ)
				#Add protein monomer
				createProduct(r,protein['wid'],1)
				#Add tRNAs and AAs
				if gene != None:
					#translating using EMBOSS transeq
					ini = gene['coordinate']-1
					direction = gene['direction']
					end = 0;
					if direction == 'f':
						end = ini+gene['length']-1
					elif direction == 'r':
						end = ini-gene['length']+1
						aux = ini
						ini = end
						end = aux
					rnaSeq = genome[ini:end]
					if direction == 'r':
						rnaSeq = rnaSeq[::-1]
					try:
						tmpSeq = open('tmpSeq.fasta.tmp','w')
					except:
						print 'Can\'t create a tmp file.'
						sys.exit(1)
					tmpSeq.write(rnaSeq)
					tmpSeq.close()
					command='transeq -table 4 -trim tmpSeq.fasta.tmp tmpProt.fasta.tmp 1>&2'
					with open(os.devnull, 'wb') as devnull:
						process = subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,stdout=devnull, stderr=subprocess.STDOUT)
					process.stdin.close()
					process.wait()
					try:
						tmpProt = open('tmpProt.fasta.tmp','r')
					except:
						print 'Can\'t create a tmp file.'
						sys.exit(1)
					tmpProt.readline()
					protSeq = tmpProt.read()		
					tmpProt.close()
					os.system('rm *.tmp')
					protSeq = ''.join(protSeq.split('\n'))
					
					#Adding tRNAs with AAs as reactants
					for key in aaCounter:
						aaCounter[key] = 0
					for aa in protSeq:
						if AMINOACIDS_INVERSE.has_key(aa):
							aaCounter[aa] = aaCounter[aa]+1			
					waterSum = 0
					for key in aaCounter:
						if aaCounter[key] > 0:
							#create tRNA with aa.
							if 'TRNA_'+AMINOACIDS_INVERSE[key]+'_TU' not in speciesList:
								createMoleculeVertice(network,'TRNA_'+AMINOACIDS_INVERSE[key]+'_TU','tRNA with '+AMINOACIDS_INVERSE[key]+' with EFTU','Aminoacylated tRNA')
								speciesList.append('TRNA_'+AMINOACIDS_INVERSE[key]+'_TU')
								speciesCounter = speciesCounter+1
							#add tRNA with aa as reactant
							createReactant(r,'TRNA_'+AMINOACIDS_INVERSE[key]+'_TU',aaCounter[key])
							waterSum = waterSum+aaCounter[key]
							#add tRNA without aa as product
							if tRNAdict.has_key(AMINOACIDS_INVERSE[key]):
								for i in range(len(tRNAdict[AMINOACIDS_INVERSE[key]])):
									createProduct(r,tRNAdict[AMINOACIDS_INVERSE[key]][i],float(aaCounter[key])/len(tRNAdict[AMINOACIDS_INVERSE[key]]))

					#Adding waters as products
					createProduct(r,H2O,waterSum)

					protSize = len(protSeq)

					#Add EFTUs and EFGs
					createProduct(r,EFTU+'_'+GDP,protSize*2)
					createReactant(r,EFG+'_'+GTP,protSize)
					createProduct(r,EFG+'_'+GDP,protSize)
					createProduct(r,PI,2*protSize+1) #EFTUs, EFGs and RRF
				

					#Create degradation reactions

					#FTSH
					decayCounter = decayCounter+1
					dr = network.createReaction()
					check(dr)
					check(dr.setId(protein['wid'].replace("-","_")+'_FTSH_DEGRADATION'))
					check(dr.setName(protein['name']+' FtsH Degradation Reaction'))
					check(dr.setReversible(False))
					check(dr.setFast(False))
					check(dr.appendAnnotation('Protein Degradation Reaction'))
					check(dr.setNotes('Protein Decay',True))
					createReactant(dr,protein['wid'],1)
					createReactant(dr,ATP,int((protSize/15)*8.3))
					createReactant(dr,H2O,int((protSize/15)*2))
					createModifier(dr,FTSH_PROTEASE)
					for pept in PEPTDASE:
						createModifier(dr,pept)
					createModifier(dr,SSRA)
					createModifier(dr,SSRA_BP)
					createProduct(dr,ADP,int((protSize/15)*8.3))
					createProduct(dr,H,int(protSize/15))
					createProduct(dr,PI,int((protSize/15)*8.3))
					for aa in aaCounter:
						if aaCounter[aa]>0:
							createProduct(dr,AMINOACIDS_INVERSE[aa],aaCounter[aa])
					for prost in pGroups:
						if prost['coefficient'] != None:
							createProduct(dr,prost['wid'],prost['coefficient'])

					#Lon
					decayCounter = decayCounter+1
					dr = network.createReaction()
					check(dr)
					check(dr.setId(protein['wid'].replace("-","_")+'_LON_DEGRADATION'))
					check(dr.setName(protein['name']+' Lon Degradation Reaction'))
					check(dr.setReversible(False))
					check(dr.setFast(False))
					check(dr.appendAnnotation('Protein Degradation Reaction'))
					check(dr.setNotes('Protein Decay',True))
					createReactant(dr,protein['wid'],1)
					createReactant(dr,ATP,int((protSize/20)*6))
					createReactant(dr,H2O,int((protSize/20)*2))
					createModifier(dr,LA_PROTEASE)
					createModifier(dr,CLP_PROTEASE)
					for pept in PEPTDASE:
						createModifier(dr,pept)
					createProduct(dr,ADP,int((protSize/20)*6))
					createProduct(dr,H,int(protSize/20))
					createProduct(dr,PI,int((protSize/20)*6))
					for aa in aaCounter:
						if aaCounter[aa]>0:
							createProduct(dr,AMINOACIDS_INVERSE[aa],aaCounter[aa])
					for prost in pGroups:
						if prost['coefficient'] != None:
							createProduct(dr,prost['wid'],prost['coefficient'])
			

	print ICcounter,'translation initiation complex formation reactions created.'
	print str(counter)+' translation reactions created.'
	print str(decayCounter)+' protein decay reactions created.'
	print str(speciesCounter)+' new molecules (tRNA with aminoacid) created.'



#Complex formation reactions --------------------------------------------------

	counter = 0
	decayCounter = 0
	inconpleteComplexDegrad = {}
	print 'Creating complex formaMG_224_9MER_GDPtion reactions...'
	complexes = mysqlQuery(c,'select public_proteincomplex.parent_ptr_protein_id, public_entry.wid, public_entry.name from public_proteincomplex inner join public_entry on public_entry.id=public_proteincomplex.parent_ptr_protein_id')
	for comp in complexes:
		if comp['wid'].replace("-","_") in complexToDelete:
			continue
		counter = counter+1
		r = network.createReaction()
		check(r)
		check(r.setId(comp['wid'].replace("-","_")+'_FORMATION'))
		check(r.setName(comp['name']+' Formation Reaction'))
		check(r.setReversible(False))
		check(r.setFast(False))
		check(r.appendAnnotation('Complex Formation Reaction'))
		check(r.setNotes('Complex Formation',True))
		reactants = mysqlQuery(c,'select public_entry.wid,public_entry.name, public_entry.model_type, public_proteincomplexbiosythesisparticipant.coefficient from public_proteincomplex_biosynthesis inner join public_proteincomplexbiosythesisparticipant on public_proteincomplexbiosythesisparticipant.id=public_proteincomplex_biosynthesis.proteincomplexbiosythesisparticipant_id inner join public_entry on public_entry.id=public_proteincomplexbiosythesisparticipant.molecule_id where proteincomplex_id='+str(comp['parent_ptr_protein_id']))
		if len(reactants)>0:
			reactantCount = 0
			productCount = 0
			for molecule in reactants:			
				if molecule['wid'].replace("-","_") not in speciesList:
					createMoleculeVertice(network,molecule['wid'].replace("-","_"),molecule['name'],"Protein Complex")
					speciesList.append(molecule['wid'].replace("-","_"))
				if molecule['coefficient']<0:
					createReactant(r,molecule['wid'],molecule['coefficient']*-1)
					reactantCount = reactantCount+1
				else:
					createProduct(r,molecule['wid'],molecule['coefficient'])
					productCount = productCount+1
			#Add RSGA in RIBOSOME_30S
			if comp['wid']== RIBOSOME_30S:
				createModifier(r,RSGA)
				createReactant(r,GTP,1+len(RB30SBF))
				createReactant(r,H2O,2+len(RB30SBF))
				createProduct(r,PI,1+len(RB30SBF))
				createProduct(r,H,2+len(RB30SBF))
				createProduct(r,GDP,1+len(RB30SBF))
				for sub in RB30SBF:
					createModifier(r,sub)
			#Ribosome 50S
			if comp['wid']== RIBOSOME_50S:
				createReactant(r,GTP,len(RB50SBF))
				createReactant(r,H2O,len(RB50SBF))
				createProduct(r,PI,len(RB50SBF))
				createProduct(r,H,len(RB50SBF))
				createProduct(r,GDP,len(RB50SBF))
				for sub in RB50SBF:
					createModifier(r,sub)

			#Create degradation reactions

			#Lon
			decayCounter = decayCounter+1
			dr = network.createReaction()
			check(dr)
			check(dr.setId(comp['wid'].replace("-","_")+'_LON_DEGRADATION'))
			check(dr.setName(comp['name'].replace("-","_")+' Lon Degradation Reaction'))
			check(dr.setReversible(False))
			check(dr.setFast(False))
			check(dr.appendAnnotation('Protein Degradation Reaction'))
			check(dr.setNotes('Protein Decay',True))
			createReactant(dr,comp['wid'],1)
			createModifier(dr,LA_PROTEASE)
			createModifier(dr,CLP_PROTEASE)
			if comp['wid']== RIBOSOME_30S or comp['wid']== RIBOSOME_50S:
				createModifier(dr,RNASEIII)
				createModifier(dr,RNASEJ)
			for molecule in reactants:
				if molecule['coefficient']<0:
					molecule['wid'] = molecule['wid'].replace("-","_")
					if molecule['model_type'] == 'ProteinComplex':
						degradReact = network.getReaction(molecule['wid']+'_LON_DEGRADATION')
						if degradReact!=None and not molecule['wid'] in inconpleteComplexDegrad:
							degradN = degradReact.getNumProducts()
							for i in range(degradN):
								degradProd = degradReact.getProduct(long(i))
								createProduct(dr,degradProd.getSpecies(),degradProd.getStoichiometry())
							degradN = degradReact.getNumReactants()
							for i in range(degradN):
								degradR = degradReact.getReactant(long(i))
								if degradR.getSpecies() != molecule:
									createReactant(dr,degradR.getSpecies(),degradR.getStoichiometry())
						else:
							if not inconpleteComplexDegrad.has_key(comp['wid'].replace("-","_")):
								inconpleteComplexDegrad[comp['wid'].replace("-","_")] = []
							inconpleteComplexDegrad[comp['wid'].replace("-","_")].append(molecule['wid'])
					elif molecule['model_type'] == 'ProteinMonomer':
						degradReact = network.getReaction(molecule['wid']+'_LON_DEGRADATION')
						if degradReact != None:
							degradN = degradReact.getNumProducts()
							for i in range(degradN):
								degradProd = degradReact.getProduct(long(i))
								createProduct(dr,degradProd.getSpecies(),degradProd.getStoichiometry())
							degradN = degradReact.getNumReactants()
							for i in range(degradN):
								degradR = degradReact.getReactant(long(i))
								if degradR.getSpecies() != molecule:
									createReactant(dr,degradR.getSpecies(),degradR.getStoichiometry())
						else:
							print "Degradation reaction not found for "+molecule['wid']
	
	#Filling the complex degradation reactions
	while(len(inconpleteComplexDegrad)):
		keysToDelete = []
		for key in inconpleteComplexDegrad:
			complexDegradReact = network.getReaction(key+'_LON_DEGRADATION')
			for molecule in inconpleteComplexDegrad[key]:
				degradReact = network.getReaction(molecule+'_LON_DEGRADATION')
				if degradReact != None and not molecule in inconpleteComplexDegrad.keys():
					degradN = degradReact.getNumProducts()
					for i in range(degradN):
						degradProd = degradReact.getProduct(long(i))
						createProduct(complexDegradReact,degradProd.getSpecies(),degradProd.getStoichiometry())
					inconpleteComplexDegrad[key].remove(molecule)
			if not len(inconpleteComplexDegrad[key]):
				keysToDelete.append(key)
		for key in keysToDelete:
			del inconpleteComplexDegrad[key]
					
					
	#Terminal organel Assembly
	counter += 1
	r = network.createReaction()
	check(r)
	check(r.setId('TERMINAL_ORGANELLE_ASSEMBLY'))
	check(r.setName('Terminal Organelle Assembly'))
	check(r.setReversible(False))
	check(r.setFast(False))
	check(r.appendAnnotation('Complex Formation Reaction'))
	check(r.setNotes('Complex Formation',True))
	
	createModifier(r,MEMBRANE)
	for reaction in organelAssemblyReactions:
		react = mysqlQuery(c,"select public_reaction.modification_id, public_reaction.enzyme_id, public_reaction.parent_ptr_species_component_id from public_entry inner join public_reaction on public_reaction.parent_ptr_species_component_id=public_entry.id where public_entry.wid='"+reaction+"'")
		react = react[0]
		modification = mysqlQuery(c,'select public_entry.wid from public_modificationreaction inner join public_entry on public_entry.id=public_modificationreaction.molecule_id where public_modificationreaction.id='+str(react['modification_id']))
		modification = modification[0]['wid']
		
		createReactant(r,modification,1)
		if react['enzyme_id'] != None:
			enzymes = mysqlQuery(c,'select public_entry.wid from public_enzymeparticipant inner join public_entry on public_entry.id=public_enzymeparticipant.protein_id where public_enzymeparticipant.id ='+str(react['enzyme_id']))
			for enz in enzymes:
				if r.getModifier(enz['wid'])==None:
					createModifier(r,enz['wid'])
	createProduct(r,TERMINAL_ORGANELLE,1)
	

	print str(counter)+' complex formation reactions created.'
	print str(decayCounter)+' complex degradation reactions created.'


#Protein Modification --------------------------------------------------------

	print 'Reading protein modification reactions...'
	proteinModCounter = 0
	for reaction in modificationReactions:
		react = mysqlQuery(c,"select public_reaction.modification_id, public_reaction.enzyme_id, public_reaction.parent_ptr_species_component_id from public_entry inner join public_reaction on public_reaction.parent_ptr_species_component_id=public_entry.id where public_entry.wid='"+reaction+"'")
		react = react[0]
		modification = mysqlQuery(c,'select public_entry.wid from public_modificationreaction inner join public_entry on public_entry.id=public_modificationreaction.molecule_id where public_modificationreaction.id='+str(react['modification_id']))
		modification = modification[0]['wid']
		reactants = mysqlQuery(c,'select public_entry.wid,public_entry.name, public_reactionstoichiometryparticipant.coefficient from public_reaction_stoichiometry inner join public_reactionstoichiometryparticipant on public_reaction_stoichiometry.reactionstoichiometryparticipant_id=public_reactionstoichiometryparticipant.id inner join public_entry on public_entry.id=public_reactionstoichiometryparticipant.molecule_id where reaction_id = '+str(react['parent_ptr_species_component_id']))

		#Protein Maturation
		if network.getSpecies(modification).getAnnotationString() == PROTEIN_TYPE:
			proteinModCounter = proteinModCounter+1
			r = network.getReaction(modification+'_MODIFICATION')
			if r == None:
				r = network.createReaction()
				check(r)
				check(r.setId(modification+'_MODIFICATION'))
				check(r.setName(modification+' Modification Reaction'))
				check(r.setReversible(False))
				check(r.setFast(False))
				check(r.appendAnnotation('Protein Modification Reaction'))
				check(r.setNotes('Protein Modification',True))
				if len(reactants)>0:
					for molecule in reactants:
						match = 0
						for aa in AMINOACIDS:
							if aa in molecule['wid']:
								match = 1
								break
						if not match:
							if molecule['coefficient']<0:
								createReactant(r,molecule['wid'],molecule['coefficient']*-1)
							else:
								createProduct(r,molecule['wid'],molecule['coefficient'])
					#Create product and pre product				
					createProduct(r,modification,1)
					if modification+'_INACTIVE' not in speciesList:
						createMoleculeVertice(network,modification+'_INACTIVE',modification+' Inactive','Inactive Protein')
						speciesList.append(modification+'_INACTIVE')
					createReactant(r,modification+'_INACTIVE',1)
					#Correct transcription reaction
					
					reactionFound = 0
					for i in range(network.getNumReactions()):
						ptrnar = network.getReaction(long(i))
						if ptrnar != None:
							if ptrnar.getAnnotationString() == TRANSLATION_TYPE:
								for j in range(ptrnar.getNumProducts()):
									gene = ptrnar.getProduct(long(j))
									if gene != None:
										if gene.getSpecies() == modification:
											ptrnar.removeProduct(modification)
											createProduct(ptrnar,modification+'_INACTIVE',1)
			else:
				if len(reactants)>0:
					for molecule in reactants:
						match = 0
						for aa in AMINOACIDS:
							if aa in molecule['wid']:
								match = 1
								break
						if not match:
							if molecule['coefficient']<0:
								createReactant(r,molecule['wid'],molecule['coefficient']*-1)
							else:
								createProduct(r,molecule['wid'],molecule['coefficient'])
			#ADD Enzymes
			if react['enzyme_id'] != None:
				enzymes = mysqlQuery(c,'select public_entry.wid from public_enzymeparticipant inner join public_entry on public_entry.id=public_enzymeparticipant.protein_id where public_enzymeparticipant.id ='+str(react['enzyme_id']))
				for enz in enzymes:
					if r.getModifier(enz['wid'])==None:
						createModifier(r,enz['wid'])

	print str(proteinModCounter)+' Protein Modification reactions created.'


#Removing protein polymerization old reactions -------------------------------

	print 'Removing Protein Polymerization and Protein Processing old reactions...'
	counter = 0
	maxi = network.getNumReactions()
	i=0
	while i<maxi:
		if i<network.getNumReactions():
			r = network.getReaction(long(i))
			if r != None:
				if PROTPOLY_PREFIX in r.getId() or DEFORMYLATION_PREFIX in r.getId() or METCLEAV_PREFIX in r.getId():
					deletedReactions.append(r.getId())
					network.removeReaction(r.getId())
					counter = counter+1
					i = i-1
		i = i+1
	print str(counter)+' protein polymerization old reactions removed.'



#Removing protein degradation old reactions -------------------------------

	print 'Removing Protein Degradation old reactions...'
	counter = 0
	maxi = network.getNumReactions()
	i=0
	while i<maxi:
		if i<network.getNumReactions():
			r = network.getReaction(long(i))
			if r != None:
				if PROTDECAY_PREFIX in r.getId():
					deletedReactions.append(r.getId())
					network.removeReaction(r.getId())
					counter = counter+1
					i = i-1
		i = i+1
	print str(counter)+' protein degradation old reactions removed.'




#Look for reactions without reactants or reagents -----------------------------
	counter = 0
	print 'Looking for reactions inconsistency...'
	maxi = network.getNumReactions()
	i=0
	while i<maxi:
		if i<network.getNumReactions():
			r = network.getReaction(long(i))
			if r!=None:
				np = r.getNumProducts()
				nr = r.getNumReactants()
				if np==0 or nr==0:
					counter = counter+1
					print 'Warning: ',
					deletedReactions.append(r.getId())
					check(network.removeReaction(i),'removing reaction '+str(r.getId())+'...')
					i=i-1
		i = i+1
	print str(counter)+' reactions removed.'


#Writing discarded reactions and molecules to file ----------------------------

	if writeDeletedReactions:
		print 'Writing discarded reactions and molecules to file.'
		delF = open('deleted_reactions.txt','w')
		for r in deletedReactions:
			delF.write(r+'\n')
		delF.close()
		delF = open('deleted_molecules.txt','w')
		for r in deletedMolecules:
			delF.write(r+'\n')
		delF.close()


#Writing SMBL to file ---------------------------------------------------------

	print "Writing network to file " +outputfile+ "..."
	try:
		f.write(writeSBMLToString(document))
	except ValueError:
		print "Error - Can\'t write SBML document to the file " +outputfile+"."
		sys.exit(1)

#Closing database connection --------------------------------------------------

	print "Closing database connection..."
	if con:
		con.close();

	print 'The whole-cell network is done.'
 
if __name__ == '__main__':
	main(sys.argv[1:])
	
	
