import json
import molmass
import networkx as nx
from collections import Counter
import re

def react_element(mol,loc,sto=1):
    return {'molecule':mol,'location':loc,'stoichiometry':sto}

def complimentary_base(b,target='DNA'):
    if b not in ['A','T','C','G']:
        print("WARNING: Base must be 'A','T','C',or 'G'. Got: "+dbase)
        return False
    if b == 'A':
        if target=='RNA':
            return 'U'
        return 'T'
    if b == 'T':
        return 'A'
    if b == 'C':
        return 'G'
    if b == 'G':
        return 'C'
    
def complimentary_sequence(seq,target='DNA'):
    new_seq = ''
    for b in seq:
        new_seq+=complimentary_base(b,target)
    return new_seq

class Jsonable(object):
    def reprJSON(self):
        d = dict()
        for a, v in self.__dict__.items():
            if (hasattr(v, "reprJSON")):
                d[a] = v.reprJSON()
            else:
                d[a] = v
        return d
    
    def __repr__(self):
        return str(self.reprJSON())
    
    def __str__(self):
        return str(self.reprJSON())
    
class ComplexEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj,'reprJSON'):
            return obj.reprJSON()
        else:
            return json.JSONEncoder.default(self, obj)

class Cell(Jsonable):
    def __init__(self):
        self.organism = None
        self.locations = []
        self.molecules = {}
        self.reactions = {}
        self.chromosome_regions_length = {}
        self.chromosome_regions = {}
        self.dna_regions = []
        self.genes = {}
        self.transcription_units = {}
        self.chromosomes = {}
        self.translation_process = None
        self.transcription_process = None
        self.replication_process = None
        self.protein_degradation_process = None
        self.rna_degradation_process = None
        self.created_molecules = {}
        self.unknown_function_rna_sequences = []
        
    
    def solve_inconcistencies(self,handler):
        mol_inconcistencies = {}
        missing_mols = []
        missing_reacts = []
        old_missing_mols = ['start']
        old_missing_reacts = ['start']
        react_inconcistencies = self.get_all_reaction_inconcistencies()
        while old_missing_mols!=missing_mols and old_missing_reacts!=missing_reacts:
            old_missing_reacts = missing_reacts
            old_missing_mols = missing_mols
            missing_mols = self.solve_reactions_inconcistencies(react_inconcistencies,handler)
            mol_inconcistencies = self.get_all_molecule_inconcistencies()
            self.solve_molecules_inconcistencies(mol_inconcistencies,handler)                        
            react_inconcistencies = self.get_all_reaction_inconcistencies() 
            missing_reacts = sorted(list(set(react_inconcistencies.keys())))
        return react_inconcistencies,mol_inconcistencies
    
    def solve_molecules_inconcistencies(self,inconcistencies,handler):
        for mol,error in inconcistencies.items():
            if error == 'not in reaction':
                self.delete_molecule(mol)
                continue
            elif error == 'biosynthesis missing':                
                if self.molecules[mol].mol_type == 'ProteinComplex':
                    if handler.have_entry(mol):
                        reaction = handler.make_biosynthesis_for_molecule(mol)
                    elif self.have_created_molecule(mol):
                        reaction = self.make_biosynthesis_for_molecule(mol)
                    if reaction == None:
                        print("Warning: no synthesis reaction retrieved for:",mol)
                        continue
                    if not issubclass(type(reaction),list):
                        reaction = [reaction]
                    self.add_reactions(reaction)
                    self.molecules[mol].biosynthesis_reaction = [r.ref_name for r in reaction]
                    continue
                elif self.molecules[mol].mol_type == 'RNAComplex':
                    reaction = None
                    for r in self.reactions.keys():
                        for reactant in self.reactions[r].molecules:
                            if mol==reactant['molecule'] and float(reactant['stoichiometry'])>0:
                                reaction = r
                                break
                        if reaction != None: break
                    self.molecules[mol].biosynthesis_reaction = reaction
                    continue
                else:
                    reaction = self.make_biosynthesis_for_molecule(mol)
                    if reaction == None:
                        #print("Warning: no synthesis reaction retrieved for:",mol)
                        continue
                    if not issubclass(type(reaction),list):
                        reaction = [reaction]
                    self.add_reactions([r for r in reaction if r.ref_name not in self.reactions.keys()])
                    self.molecules[mol].biosynthesis_reaction = [r.ref_name for r in reaction]                
            elif error == 'degradation missing':
                reaction = self.make_degradation_for_molecule(mol)
                if reaction == None:
                    #print("Warning: no synthesis reaction retrieved for:",mol)
                    continue
                if issubclass(type(reaction),str):
                    self.molecules[mol].degradation_reaction = reaction
                    continue
                if not issubclass(type(reaction),list):
                    reaction = [reaction]
                self.add_reactions([r for r in reaction if r.ref_name not in self.reactions.keys()])
                self.molecules[mol].degradation_reaction = [r.ref_name for r in reaction]
            elif error == 'composition missing':
                comp = self.get_composition_for_molecule(mol)
                self.molecules[mol].composition = comp
            elif error == 'molecular-weight missing':
                mw = self.get_mw_for_molecule(mol)
                self.molecules[mol].mw = mw
    
    
    def solve_reactions_inconcistencies(self,inconcistencies,handler):
        missing_mols = []
        mol_created = []
        for inc in inconcistencies.values():
            for mol,error in inc.items():
                if error == 'missing molecule':
                    missing_mols.append(mol)
                    if mol not in mol_created:
                        mol_obj = None
                        if handler.have_entry(mol) or handler.have_created_molecule(mol):
                            mol_obj = handler.make_molecule_from_entry_name(mol)
                        elif self.have_created_molecule(mol):
                            mol_obj = self.make_molecule_from_name(mol)
                        if mol_obj!=None and mol_obj!=False:
                            self.add_molecule(mol_obj)
                            mol_created.append(mol)
        return sorted(list(set(missing_mols)))
    
    
    def set_genome_reference(self,genome_name,genome_ref):
        self.genome_reference[genome_name] = gen_ref
       
    def set_chromosome_region_length(self,regions):
        for chrm in self.chromosomes.keys():
            if chrm not in regions.keys():
                print("ERROR: Missing chromosome's number of regions. Chromosome:",chrm)
                return False
        self.chromosome_regions_length = regions
        return True
        
    def get_chromosomes(self):
        return self.chromosomes.keys()
    
    def add_chromosome_features(self,feats):
        if not issubclass(type(feats),dict):
            print("ERROR: feats must be a 'dict' object. Got:",type(feats))
            return False
        chrms = self.get_chromosomes()
        for feat_name,f in feats.items():
            if f['chromosome'] in chrms:
                self.chromosomes[f['chromosome']].add_feature(feat_name,f)
            else:
                print("WARNING: chromosome not found. Chromosome:",f['chromosome'])
        return True 
    
    def get_replication_terminus_dnas(self):
        return [self.get_chromosome_region_from_coordinate(chrm,self.chromosomes[chrm].features[feat]['start'])+"_"+self.transcription_process.dna_sufix for chrm,feat in self.replication_process.replication_terminus_feature.items()]
        
    def have_created_molecule(self,mol):
        return mol in self.created_molecules.keys()
    
    def get_dna_subsequence(self,chrm,start,length,direction):
        if chrm not in self.chromosomes.keys():
            print("ERROR: Chromosome does not exist. Got:",chrm)
            return ''
        start -=1
        if direction == 'Forward':
            if start+length > len(self.chromosomes[chrm].sequence):
                print("WARNING: invalid range of genome sequence. Got Start:",start,"Length:",length,"Direction:",direction)
                return ''
            return self.chromosomes[chrm].sequence[start:start+length]
        elif direction == 'Reverse':
            if start <= length:
                print("WARNING: invalid range of genome sequence. Got Start:",start,"Length:",length,"Direction:",direction)
                return ''
            return complimentary_sequence(self.chromosomes[chrm].sequence[start:start+length][::-1])
        else:
            print("ERROR: Invalid direction, must be 'Forward' or 'Reverse'. Got:",direction)
            return ''
    
    def create_chromosome_regions(self):
        for chrm_name,chrm in self.chromosomes.items():
            chrm.offset = len(self.dna_regions)
            max_length = self.chromosome_regions_length[chrm_name]
            features = sorted(list(chrm.features.values()), key=lambda k: k['start'])
            for feat in features:
                feat['ref_name'] = [feat['ref_name']]
            i = 1
            while i<len(features):
                if features[i]['start'] < features[i-1]['start']+features[i-1]['length']:
                    features[i-1]['length'] = features[i-1]['length']+features[i]['length']
                    features[i-1]['ref_name'].append(features[i]['ref_name'])
                    features.pop(i)
                else:
                    i += 1                
            chrm_reg_number = 0
            replication_position = "Left"
            pos = 1
            while pos<features[0]['start']-1:                
                reg = ChromosomeRegion()
                reg.ref_name = chrm.ref_name+"_region_"+str(chrm_reg_number)
                reg.chromosome = chrm.ref_name
                reg.start = pos
                reg.length = min(max_length,features[0]['start']-pos)
                reg.end = reg.start-1+reg.length
                pos = reg.end
                reg.sequence = self.get_dna_subsequence(chrm_name,reg.start,reg.length,"Forward")
                reg.replication_position = replication_position
                reg.previous_region = chrm.ref_name+"_region_"+str(chrm_reg_number-1)
                reg.next_region = chrm.ref_name+"_region_"+str(chrm_reg_number+1)
                self.dna_regions.append(reg)
                self.chromosome_regions[self.dna_regions[-1].ref_name] = self.dna_regions[-1]
                chrm_reg_number += 1
            for i in range(len(features)):
                reg = ChromosomeRegion()
                reg.ref_name = chrm.ref_name+"_region_"+str(chrm_reg_number)
                reg.chromosome = chrm.ref_name
                reg.start = features[i]['start']
                reg.length = features[i]['length']
                reg.end = reg.start-1+reg.length
                reg.sequence = self.get_dna_subsequence(chrm_name,reg.start,reg.length,"Forward")
                reg.features = features[i]['ref_name']
                reg.is_replication_origin = self.replication_process.replication_start_feature[chrm.ref_name] in features[i]['ref_name']
                if reg.is_replication_origin:
                    replication_position = "Right"
                reg.is_replication_terminus = self.replication_process.replication_terminus_feature[chrm.ref_name] in features[i]['ref_name']
                reg.replication_position = replication_position
                if reg.is_replication_terminus:
                    reg.replication_position = "End"
                    replication_position = "Left"
                reg.previous_region = chrm.ref_name+"_region_"+str(chrm_reg_number-1)
                reg.next_region = chrm.ref_name+"_region_"+str(chrm_reg_number+1)
                self.dna_regions.append(reg)
                self.chromosome_regions[self.dna_regions[-1].ref_name] = self.dna_regions[-1]
                chrm_reg_number += 1
                pos = reg.end+1
                if i+1<len(features):                    
                    next_stop = features[i+1]['start']
                else:
                    next_stop = chrm.length               
                while pos<next_stop:
                    reg = ChromosomeRegion()
                    reg.ref_name = chrm.ref_name+"_region_"+str(chrm_reg_number)
                    reg.chromosome = chrm.ref_name
                    reg.start = pos
                    reg.length =  next_reg_length = min(next_stop-pos,max_length)
                    reg.end = reg.start-1+reg.length
                    pos = reg.end+1
                    reg.sequence = self.get_dna_subsequence(chrm_name,reg.start,reg.length,"Forward")
                    reg.replication_position = replication_position
                    reg.previous_region = chrm.ref_name+"_region_"+str(chrm_reg_number-1)
                    reg.next_region = chrm.ref_name+"_region_"+str(chrm_reg_number+1)
                    self.dna_regions.append(reg)
                    self.chromosome_regions[self.dna_regions[-1].ref_name] = self.dna_regions[-1]
                    chrm_reg_number += 1
                
            self.dna_regions[-1].next_region = self.dna_regions[chrm.offset].ref_name
            self.dna_regions[chrm.offset].previous_region = self.dna_regions[-1].ref_name
            chrm.n_regions = chrm_reg_number
            
            for i  in range(1,len(self.dna_regions)):
                if self.dna_regions[i].start<=self.dna_regions[i-1].end:
                    print(self.dna_regions[i-1])
                    print(self.dna_regions[i])
            
        for tr in self.transcription_units.values():
            chrm_region_end_transcription = self.get_chromosome_region_from_coordinate(tr.chromosome,tr.transcription_end)
            chrm_region_start_transcription = self.get_chromosome_region_from_coordinate(tr.chromosome,tr.transcription_start)
            if tr.direction == "Reverse":
                tmp = chrm_region_end_transcription
                chrm_region_end_transcription = chrm_region_start_transcription
                chrm_region_start_transcription = tmp
            self.chromosome_regions[chrm_region_end_transcription].transcription_unit_end.append(tr.ref_name)
            self.chromosome_regions[chrm_region_start_transcription].transcription_unit_start.append(tr.ref_name)
        
            
            
    def get_chromosome_region_from_coordinate(self,chrm,pos):
        if chrm not in self.chromosomes.keys():
            print("ERROR: chromosome not found. Chromosome:",chrm)
            return False
        if pos>self.chromosomes[chrm].length:
            print("ERROR: invalid position "+str(pos)+". Chromosome:",chrm,"Length:",self.chromosomes[chrm].length)
            return False
        min_reg = self.chromosomes[chrm].offset
        max_reg = self.chromosomes[chrm].offset+self.chromosomes[chrm].n_regions
        regs = range(min_reg,max_reg,1)
        length = len(regs)
        while length>1:
            idx = int((length-1)/2)+length%2
            if pos < self.dna_regions[regs[idx]].start:
                regs = regs[:idx]
            elif pos > self.dna_regions[regs[idx]].end:
                regs = regs[idx+1:]
            else:
                return self.dna_regions[regs[idx]].ref_name
            length = len(regs)
        
        return self.dna_regions[regs[0]].ref_name
    
    def add_chromosome(self,chrm_obj):
        if not issubclass(type(chrm_obj),Chromosome):
            print("ERROR: object must be 'Chromosome' class. Got:",type(chrm_obj))
            return False
        if chrm_obj.ref_name in self.chromosomes.keys():
            print("ERROR: 'ref_name' already exists:",chrm_obj.ref_name)
            return False
        self.chromosomes[chrm_obj.ref_name] = chrm_obj
        return True
        
    def add_chromosomes(self,chrm_obj_list):
        if type(chrm_obj_list) != list:
            print("ERROR: 'chrm_obj_list' must be a list of 'Chromosome' class objects. Got:",type(chrm_obj_list))
            return False         
        return sum([self.add_chromosome(chrm_obj) for chrm_obj in chrm_obj_list])
    
    def add_transcription_unit(self,tu_obj):
        if not issubclass(type(tu_obj),TranscriptionUnit):
            print("ERROR: object must be 'TranscriptionUnit' class. Got:",type(tu_obj))
            return False
        if tu_obj.ref_name in self.transcription_units.keys():
            print("ERROR: 'ref_name' already exists:",tu_obj.ref_name)
            return False
        tu_obj.chromosome = self.genes[tu_obj.genes[0]].chromosome
        genes = [self.genes[g] for g in tu_obj.genes]
        if tu_obj.is_polycistronic:
            tu_obj.is_cleaved = True
            for g in genes:
                if g.type == 'mRNA':
                    tu_obj.is_cleaved = False
                    break
        tu_obj.direction = genes[0].direction
        genes.sort(key=lambda x: x.start, reverse=False)
        tu_obj.transcription_start = genes[0].start
        tu_obj.transcription_end = genes[-1].start + genes[-1].length
        tu_obj.transcription_length = tu_obj.transcription_end-tu_obj.transcription_start
        tu_obj.genes = [g.ref_name for g in genes]
        self.transcription_units[tu_obj.ref_name] = tu_obj
        return True
    
    def add_transcription_units(self,tu_obj_list):
        if type(tu_obj_list) != list:
            print("ERROR: 'tu_obj_list' must be a list of 'TranscriptionUnit' class objects. Got:",type(tu_obj_list))
            return False         
        return sum([self.add_transcription_unit(tu_obj) for tu_obj in tu_obj_list])
    
    def add_gene(self,gen_obj):
        if not issubclass(type(gen_obj),Gene):
            print("ERROR: object must be 'Gene' class. Got:",type(gen_obj))
            return False
        if gen_obj.ref_name in self.genes.keys():
            print("ERROR: 'ref_name' already exists:",gen_obj.ref_name)
            return False
        self.genes[gen_obj.ref_name] = gen_obj
        return True
        
    def add_genes(self,gen_obj_list):
        if type(gen_obj_list) != list:
            print("ERROR: 'gen_obj_list' must be a list of 'Gene' class objects. Got:",type(gen_obj_list))
            return False         
        return sum([self.add_gene(gen_obj) for gen_obj in gen_obj_list])
    
    def add_all_DnaA_complexes(self,feature_prefix):        
        active_dnaa = self.replication_process.dnaA_complex[-1]['molecule']
        active_sufix = "_"+self.transcription_process.dna_sufix+"_"+active_dnaa
        inactive_dnaa = self.replication_process.inactive_dnaA_complex[-1]['molecule']
        inactive_sufix = "_"+self.transcription_process.dna_sufix+"_"+active_dnaa
        sufixes = [active_sufix,inactive_sufix]
        regions = []
        for reg in self.chromosome_regions.values():
            have_dnaa_box = False
            for feat in reg.features:
                if feature_prefix in feat:
                    have_dnaa_box = True
                    break
            if have_dnaa_box:
                regions.append(reg.ref_name)
        for sufix in sufixes:        
            dnaa_complexes = [self.make_dna_protein_complex(reg+sufix) for reg in regions if reg+sufix not in self.molecules.keys()]
            for m in dnaa_complexes:
                self.created_molecules[m.ref_name] = "DNA DnaA Complex"
            self.add_molecules(dnaa_complexes)
        return True
    
    def add_all_dna_condensed_states(self,minimal_region_length):
        prot = self.replication_process.get_dna_condensation_protein()['molecule']
        sufix = "_"+self.transcription_process.dna_sufix+"_"+prot
        regions = []
        for reg in self.chromosome_regions.values():
            if reg.length >= minimal_region_length:
                regions.append(reg.ref_name)
        condensed_complexes = [self.make_dna_protein_complex(reg+sufix) for reg in regions]
        for m in condensed_complexes:
            self.created_molecules[m.ref_name] = "Condensed DNA Complex"
        self.add_molecules(condensed_complexes)
        return True
            
        
    def add_molecule(self,mol_obj):
        if not issubclass(type(mol_obj),Molecule):
            print("ERROR: object must inherit 'Molecule' class. Got:",type(mol_obj))
            return False
        if mol_obj.ref_name in self.molecules.keys():
            print("ERROR: 'ref_name' already exists:",mol_obj.ref_name)
            return False
        self.molecules[mol_obj.ref_name] = mol_obj
        return True
        
    def add_molecules(self,mol_obj_list):
        if type(mol_obj_list) != list:
            print("ERROR: 'mol_obj_list' must be a list of 'Molecule' class inherited objects. Got:",type(mol_obj_list))
            return False         
        return sum([self.add_molecule(mol_obj) for mol_obj in mol_obj_list])
    
    def delete_molecule(self,m):
        if m in self.molecules.keys():
            del self.molecules[m]
        else:
            print('WARNING: Tried to delete molecule, molecule does not exist. Mol:',m)
        
    def add_reaction(self,react_obj):
        if not issubclass(type(react_obj),Reaction):
            print("ERROR: object must inherit 'Reaction' class. Got:",type(react_obj))
            return False
        if react_obj.ref_name in self.reactions.keys():
            print("ERROR: 'ref_name' already exists:",react_obj.ref_name)
            return False
        self.solve_reaction_redundancy(react_obj)
        self.reactions[react_obj.ref_name] = react_obj
        return True
        
    def add_reactions(self,react_obj_list):
        if type(react_obj_list) != list:
            print("ERROR: 'react_obj_list' must be a list of 'Reaction' class inherited objects. Got:",type(react_obj_list))
            return False         
        return sum([self.add_reaction(react_obj) for react_obj in react_obj_list])
        
    def save(self,path,compressed=False):
        try:
            with open(path,'w') as f:
                if compressed:
                    f.write(json.dumps(self.reprJSON(), cls=ComplexEncoder))
                else:
                    f.write(json.dumps(self.reprJSON(), cls=ComplexEncoder, indent=4))
                return True
        except:
            print("ERROR: Failed to write data in:",path)
            return False 
    
    def generate_network(self,file=""):
        G =  nx.DiGraph()
        for r in self.reactions.keys():
            n,attr = self.reactions[r].get_node()
            G.add_node(n,**attr)
            if self.reactions[r].is_reversible:
                G.add_node(n+"_rev",**attr)
            for m,loc,mloc in self.reactions[r].get_molecules_and_locations():
                attr = self.molecules[m].get_node_attr()
                attr['compartment'] = loc
                if mloc not in G.nodes:
                    G.add_node(mloc,**attr)
            G.add_edges_from([e for e in self.reactions[r].get_edges() if e[0] in G.nodes and e[1] in G.nodes])
            if self.reactions[r].is_reversible:
                G.add_edges_from([e for e in self.reactions[r].get_edges(rev=True) if e[0] in G.nodes and e[1] in G.nodes])
        for n in G.nodes():
            if G.nodes[n]['type']=='m':
                if G.nodes[n]['moltype'] == "ProteinComplex" and "DNA" in n:
                    G.nodes[n]['moltype'] = "DNA-Protein Complex"
                    continue
                word = re.sub(r"([a-z])([A-Z])", r"\1 \2", G.nodes[n]['moltype'])
                word = re.sub(r"([A-Z])([A-Z])([a-z])", r"\1 \2\3", word)
                G.nodes[n]['moltype'] = word
        if file == "":
            return G
        else:
            strs = file.split(".")
            if len(strs)>=2:
                if strs[-1] == "gml":
                    return nx.write_gml(G,file)                    
                elif strs[-1] == "graphml":
                    return nx.write_graphml(G,file)
                elif strs[-1] == "sbml":
                    return generate_sbml_model(G,file)
                else:
                    print("Extension not suported.")
                    return False
            else:
                print("Extension not informed.")
                return False
    
    def solve_reaction_redundancy(self,react):
        enz = list(set([e['location']+"__"+e['molecule'] for e in react.enzymes]))
        if len(enz) < len(react.enzymes):
            new_enz = []
            for e in enz:
                info = e.split("__")
                sto = 0
                for old_enz in react.enzymes:
                    if old_enz['location'] == info[0] and old_enz['molecule'] == info[1]:
                        sto += old_enz['stoichiometry']
                new_enz.append(react_element(info[1],info[0],sto))
            react.enzymes = new_enz
        additional_aux = []
        mols = list(set([e['location']+"__"+e['molecule'] for e in react.molecules]))
        if len(mols) < len(react.molecules):
            new_mols = []
            for m in mols:
                info = m.split("__")
                sto_reactant = 0
                sto_product = 0
                sto_total = 0
                for old_mol in react.molecules:
                    if old_mol['location'] == info[0] and old_mol['molecule'] == info[1]:
                        if old_mol['stoichiometry']<0:
                            sto_reactant += old_mol['stoichiometry']
                        elif old_mol['stoichiometry']>0:
                            sto_product += old_mol['stoichiometry']
                        sto_total += abs(old_mol['stoichiometry'])
                new_mols.append(react_element(info[1],info[0],sto_reactant+sto_product))
                left_mol = sto_total-abs(sto_reactant+sto_product)
                if left_mol:
                    additional_aux.append(react_element(info[1],info[0],left_mol))
            react.molecules = new_mols
        total_aux = additional_aux+react.auxiliar_molecules
        auxs = list(set([e['location']+"__"+e['molecule'] for e in total_aux]))
        if len(auxs) < len(total_aux):
            new_aux = []
            for a in auxs:
                info = a.split("__")
                sto = 0
                for old_aux in total_aux:
                    if old_aux['location'] == info[0] and old_aux['molecule'] == info[1]:
                        sto += old_aux['stoichiometry']
                new_aux.append(react_element(info[1],info[0],sto))
            react.auxiliar_molecules = new_aux
        return True
        
        
    def get_reaction_inconcistencies(self,r):
        inc = {}
        react = self.reactions[r]
        for mol in react.molecules:
            if not issubclass(dict,type(mol)):
                print(react)
            if mol['molecule'] not in self.molecules.keys():
                inc[mol['molecule']] = 'missing molecule'
        for mol in react.enzymes:
            if mol['molecule'] not in self.molecules.keys():
                inc[mol['molecule']] = 'missing molecule'
        for mol in react.auxiliar_molecules:
            if mol['molecule'] not in self.molecules.keys():
                inc[mol['molecule']] = 'missing molecule'
        return inc
                
    def get_all_reaction_inconcistencies(self):
        inc = {r:self.get_reaction_inconcistencies(r) for r in self.reactions.keys()}
        return {key:inc[key] for key in inc if len(inc[key])}    
    
    def get_molecule_inconcistencies(self,m):
        mol = self.molecules[m]
        if mol.get_type() == "SmallMolecule":
            in_reaction = False
            for r in self.reactions.keys():
                for rm in self.reactions[r].molecules:
                    if m == rm['molecule']:
                        in_reaction = True
                        break
                if in_reaction: break
            if not in_reaction: return 'not in reaction'
        if hasattr(mol, 'biosynthesis_reaction'):
            if mol.biosynthesis_reaction == None:
                return 'biosynthesis missing'
        if hasattr(mol, 'composition'):
            if not len(mol.composition):
                return 'composition missing'
        if hasattr(mol, 'degradation_reaction'):
            if mol.degradation_reaction == None:
                return 'degradation missing'
        if hasattr(mol, 'mw'):
            if not mol.mw:
                return 'molecular-weight missing'
        return None
    
    
    def get_all_molecule_inconcistencies(self):
        inc = {m:self.get_molecule_inconcistencies(m) for m in self.molecules.keys()}
        return {key:inc[key] for key in inc if inc[key]!=None}
    
    
    def compute_reaction_mass_balance(self,r):
        if len(self.get_reaction_inconcistencies(r)):
            return False
        react = self.reactions[r]
        w = sum([self.molecules[m['molecule']].mw*m['stoichiometry'] for m in react.molecules])
        react.mass_balance = w
        return True
        
    def compute_all_reaction_mass_balance(self):
        return all([self.compute_reaction_mass_balance(r) for r in self.reactions.keys()])
    
    def set_translation_process(self,trans):
        if not issubclass(type(trans),ProteinTranslation):
            print("ERROR: translation must be a 'ProteinTranslation' object. Got:",type(trans))
            return False
        self.translation_process = trans
        return True
    
    def set_transcription_process(self,trans):
        if not issubclass(type(trans),RNATranscription):
            print("ERROR: trans must be a 'RNATranscription' object. Got:",type(trans))
            return False
        self.transcription_process = trans
        return True
    
    def set_replication_process(self,repl):
        if not issubclass(type(repl),DNAReplication):
            print("ERROR: repl must be a 'DNAReplication' object. Got:",type(repl))
            return False
        self.replication_process = repl
        return True
    
    def set_protein_degradation_process(self,degr):
        if not issubclass(type(degr),ProteinDegradation):
            print("ERROR: degr must be a 'ProteinDegradation' object. Got:",type(degr))
            return False
        self.protein_degradation_process = degr
        return True
    
    def set_rna_degradation_process(self,degr):
        if not issubclass(type(degr),RNADegradation):
            print("ERROR: degr must be a 'RNADegradation' object. Got:",type(degr))
            return False
        self.rna_degradation_process = degr
        return True
    
    def remove_molecules_from_composition(self,mol1,mol2):
        if mol1 not in self.molecules.keys():
            print("ERROR: Molecule not found:",mol1)
            return False
        mol = self.molecules[mol1]
        if not hasattr(mol,"composition"):
            print("ERROR: Molecule is not a macromolecule:",mol1)
            return False
        if type(mol2) == str:
            mol2 = [mol2]
        for m in mol2:
            if m not in self.molecules.keys():
                print("ERROR: Molecule not found:",mol_name)
                return False
            m = '-'+m
            mol.composition.append(m)
        mol.mw = self.get_mw_for_molecule(mol1)
        mol.degradation_reaction = None
        return True
        
    
    def get_mw_for_molecule(self,mol):
        mol = self.molecules[mol]
        if hasattr(mol, 'composition'):
            mws = []
            for m in mol.composition:
                factor=1
                mol_name = m
                if m[0] == '-':
                    factor = -1
                    mol_name = mol_name[1:]
                if mol_name not in self.molecules.keys():
                    break
                if self.molecules[mol_name].mw == 0.0:
                    break
                mws.append(self.molecules[mol_name].mw*factor)
            if len(mws)==len(mol.composition):
                mw = sum(mws)
                if mol.mol_type in ["ProteinMonomer","ImatureProteinMonomer","Peptide"]:
                    n_prost = 0
                    if hasattr(mol,"prosthetic_groups"):
                        n_prost = sum([m['stoichiometry'] for m in mol.prosthetic_groups])
                    mw -= (len(mol.composition)-1-n_prost)*self.molecules['H2O'].mw
                if mol.mol_type in ["RNA","ImatureRNA"]:
                    mw -= (len(mol.composition))*self.molecules['OH'].mw
                if mol.mol_type in ["DNA"]:
                    mw -= (len(mol.composition))*self.molecules['OH'].mw
                return mw
        return 0.0
    
    def get_composition_for_molecule(self,mol):
        if self.molecules[mol].mol_type in ['RNA']:
            return self.get_rna_composition(mol)
        return None
        
    def get_rna_composition(self,mol):
        if mol in self.created_molecules.keys() and self.created_molecules[mol] == "Unknown Function RNA":
            seq = mol.replace("RNA_","")
            return self.transcription_process.get_rna_composition(seq)
        gene_seq = None
        direction = None
        if self.transcription_process.tu_transcript_sufix in mol:
            tu_name = mol.replace("_"+self.transcription_process.tu_transcript_sufix,"")
            tu = self.transcription_units[tu_name]
            gene_seq = self.get_dna_subsequence(*tu.get_coordinates())
            direction = tu.direction
        else:
            gene = self.genes[mol]
            gene_seq = self.get_dna_subsequence(*gene.get_gene_coordinates())
            direction = gene.direction
        return self.transcription_process.get_rna_composition(self.transcription_process.transcribe_gene(gene_seq))
    
    def make_degradation_for_molecule(self,mol):
        if hasattr(self.molecules[mol],'localization'):
            if self.molecules[mol].localization=='e':
                return "No degradation"
        if self.molecules[mol].degradation_reaction != None:
            return self.molecules[mol].degradation_reaction
        if self.molecules[mol].mol_type in ['ProteinMonomer','ImatureProteinMonomer','Peptide']:
            return self.make_degradation_reaction_for_protein(mol)
        if self.molecules[mol].mol_type in ['RNA','ImatureRNA','SmallRNA','RNAComplex']:
            return self.make_degradation_reaction_for_rna(mol)
        if self.molecules[mol].mol_type in ['TranscriptionComplex']:
            return self.make_transcription_stall_reaction(mol)
        if self.molecules[mol].mol_type in ['TranslationComplex']:
            if self.created_molecules[mol] == "Translation Complex":
                return self.make_translation_stall_reaction(mol)
            if self.created_molecules[mol] == "Stalled Translation Complex":
                return self.make_stalled_translation_elongation_reaction(mol)
        if self.molecules[mol].mol_type in ['ProteinComplex']:
            if mol in self.created_molecules.keys():
                if self.created_molecules[mol] == "Condensed DNA Complex":
                    return self.make_dna_relaxation_reaction(mol)   
            if self.molecules[mol].biosynthesis_reaction!=None:
                if self.reactions[self.molecules[mol].biosynthesis_reaction[0]].is_reversible:
                    return "No degradation"
                else:                
                    return self.make_complex_desintegration_reaction(mol)
                         
        return None
    
    def make_complex_desintegration_reaction(self,mol):
        prot = self.molecules[mol]
        react = ComplexBiosynthesisReaction()
        react.ref_name = mol+"_Desintegration_Reaction"
        react.name = react.ref_name.replace("_"," ")
        react.process = "Macromolecular Complexation"
        react.annotations = ['Protein Complex Desintegration Reaction']
        react.molecules.append(react_element(mol,'c',-1))
        for mol_name in prot.composition:
            loc = 'c'
            local_mol_name = mol_name
            factor = 1
            if mol_name[0] == '-':
                factor = -1
                local_mol_name = local_mol_name[1:]
            if local_mol_name in self.molecules.keys() and self.molecules[local_mol_name].mol_type=="ProteinMonomer":
                loc = self.molecules[local_mol_name].localization
            react.molecules.append(react_element(local_mol_name,loc,factor))
        react.is_reversible = False
        return react
    
    def make_dna_relaxation_reaction(self,mol_name):
        txt = mol_name.split("_"+self.transcription_process.dna_sufix+"_")
        chrm_region_name = txt[0]
        prot_name = txt[1]
        dna = chrm_region_name+"_"+self.transcription_process.dna_sufix
        react = ProteinDNABinding()
        react.ref_name = dna+"_Relaxation_Reaction"
        react.name = react.ref_name.replace("_"," ")
        react.process = "Protein-DNA Interaction"
        react.annotations = ['DNA Relaxation Reaction']
        react.molecules.append(react_element(dna,'c'))
        react.molecules.append(react_element(mol_name,'c',-1))
        for m in self.replication_process.get_dna_condensation_unbounded_protein():
            react.molecules.append(m)
        react.chromosome_region = chrm_region_name
        react.binder = prot_name
        react.is_reversible = False
        return react
    
    def make_stalled_translation_elongation_reaction(self,mol_name):
        tmRNA_aa = self.translation_process.get_translation_stall_tmRNA()[0]['molecule']
        aa = tmRNA_aa.split("_")[-1]
        tmRNA = tmRNA_aa.replace("_"+aa,"")
        chrm = self.chromosomes[self.genes[tmRNA].chromosome]
        feat = chrm.features[self.translation_process.tmRNA_proteolysis_tag_chromosome_feature]
        gene_seq = self.get_dna_subsequence(chrm.ref_name,feat['start'],feat['length'],feat['direction'])
        trans = self.translation_process.translate_gene(gene_seq,trans_out='aatrna',mark_initial_codon=False)
        aa_trna = Counter(trans)
        trans = self.translation_process.translate_gene(gene_seq,trans_out='trna',mark_initial_codon=False)
        trna = Counter(trans)
        aminoacids = [aa]+self.translation_process.translate_gene(gene_seq,trans_out='aa',mark_initial_codon=False)
        
        react = TranslationElongationReaction()
        react.ref_name = mol_name+"_Translation_Elongation"
        react.name = react.ref_name.replace("_"," ")
        react.process = "Translation"
        react.annotations = ['Translation Elongation Reaction']
        react.is_reversible = False
        stalled_complex = self.translation_process.translation_complex_prefix+"_"+ self.translation_process.stalled_translation_complex_sufix
        react.molecules.append(react_element(stalled_complex,'c',-1))
        react.molecules.append(react_element(tmRNA,'c'))
        for aat,sto in aa_trna.items():
            react.molecules.append(react_element(aat,'c',-sto))
        for t,sto in trna.items():
            react.molecules.append(react_element(t,'c',sto))
        for m in self.translation_process.get_elongation_auxiliaries():
            react.auxiliar_molecules.append(m)
        for m in self.translation_process.get_translation_complex_constitution():
            react.molecules.append(m)
        for m in self.translation_process.get_elongation_energy_molecules():
            m['stoichiometry'] *= len(trans)
            react.molecules.append(m)
        react.molecules.append(react_element('H2O','c',len(trans)))
        
        
        pep = Peptide()
        pep.annotations.append('Proteolysis Tagged Peptide')
        pep.composition = aminoacids
        pep.ref_name = self.translation_process.peptide_prefix+"_"+"".join([self.translation_process.aminoacid_letter[aa] for aa in pep.composition])
        pep.name = pep.ref_name.replace("_"," ")
        pep.is_secreted = False
        pep.localization = 'c'
        pep.biosynthesis_reaction = [react.ref_name]
        
        self.add_molecule(pep)
        self.created_molecules[pep.ref_name] = "Proteolysis Tagged Peptide"
        react.molecules.append(react_element(pep.ref_name,'c'))
            
        return react
    
    def make_translation_stall_reaction(self,mol_name):
        prot_name = mol_name.replace(self.translation_process.translation_complex_prefix+"_","")
        prot = self.molecules[prot_name]
        react = ProteinDegradationReaction()
        react.ref_name = mol_name+"_Stall"
        react.name = react.ref_name.replace("_"," ")
        react.process = "Translation"
        react.annotations = ['Translation Stall Reaction']
        tu = self.transcription_units[self.genes[prot.gene].transcription_unit]
        if tu.is_polycistronic:
            rna = tu.ref_name+"_"+self.transcription_process.tu_transcript_sufix
        else:
            rna = prot.gene
        react.molecules.append(react_element(rna,'c'))
        for m in self.translation_process.get_translation_stall_auxiliaries():
            react.auxiliar_molecules.append(m)
        for m in self.translation_process.get_translation_stall_tmRNA():
            m['stoichiometry'] *= -1
            react.molecules.append(m)
        react.molecules.append(react_element(mol_name,'c',-1))
        stalled_complex = self.translation_process.translation_complex_prefix+"_"+self.translation_process.stalled_translation_complex_sufix
        react.molecules.append(react_element(stalled_complex,'c'))
        if stalled_complex not in self.molecules.keys():
            mol = TranslationComplex()
            mol.ref_name = stalled_complex
            mol.name = mol.ref_name.replace("_"," ")
            for m in self.translation_process.get_translation_complex_constitution():
                mol.composition.append(m['molecule'])
            for m in self.translation_process.get_translation_stall_tmRNA():
                mol.composition.append(m['molecule'])
            mol.biosynthesis_reaction = [react.ref_name]
            self.add_molecule(mol)
            self.created_molecules[stalled_complex] = "Stalled Translation Complex"
        else:
            self.molecules[stalled_complex].biosynthesis_reaction.append(react.ref_name)
        return react
        
        
    
    def make_transcription_stall_reaction(self,mol):
        if not mol in self.molecules.keys():
            return None
        txt = mol.replace("_"+self.transcription_process.transcription_complex_sufix,"").replace("_"+self.transcription_process.active_transcription_complex_sufix,"").split("_"+self.transcription_process.dna_sufix+"_")
        chrm_region_name = txt[0]
        tu_name = txt[1]
        chrm_region = self.chromosome_regions[chrm_region_name]
        tu = self.transcription_units[tu_name]
        react = RNADegradationReaction()
        react.ref_name = mol+"_Stall"
        react.name = react.ref_name.replace("_"," ")
        react.process = "Transcription"
        react.annotations = ['Transcription Stall Reaction']
        react.molecules.append(react_element(mol,'c',-1))
        if tu.direction == "Forward":
            next_region = chrm_region.next_region
        else:
            next_region = chrm_region.previous_region
        react.auxiliar_molecules.append(react_element(next_region+ "_"+ self.replication_process.active_replication_complex_sufix,'c'))
        
        for m in self.molecules[mol].composition:
            react.molecules.append(react_element(m,'c'))
            if "RNA_" in m and m in self.created_molecules.keys() and self.created_molecules[m] == "Unknown Function RNA" and m not in self.molecules.keys():                
                rna = self.make_rna(m)
                rna.biosynthesis_reaction = react.ref_name
                self.add_molecule(rna)                
        return react
        
        
    
    def make_degradation_reaction_for_rna(self,mol):
        rna = self.molecules[mol]
        aa = None
        if rna.mol_type == 'RNAComplex':
            for m in rna.composition:
                if self.molecules[m].mol_type == 'RNA':
                    trna = m
                else:
                    aa = m
            rna = self.molecules[trna]
        if not len(rna.composition): return None
        react = RNADegradationReaction()
        react.ref_name = mol+"_Degradation"
        react.name = react.ref_name.replace("_"," ")
        react.process = "RNA Degradation"
        react.annotations = ['RNA Degradation Reaction']
        react.molecules.append(react_element(mol,'c',-1))
        if aa != None:
            react.molecules.append(react_element(aa,'c'))
            for enz in self.rna_degradation_process.get_aminoacylated_trna_degradation_enzymes():
                react.enzymes.append(enz)
        comp = Counter(rna.composition)
        for m,s in comp.items():
            react.molecules.append(react_element(m,'c',s))
        react.molecules.append(react_element("H2O",'c',-len(rna.composition)))
        react.molecules.append(react_element("H",'c',len(rna.composition)))
        for enz in self.rna_degradation_process.get_degradation_enzymes():
                react.enzymes.append(enz)
        return react        
        
    
    def make_degradation_reaction_for_protein(self,mol):
        prot = self.molecules[mol]
        react = ProteinDegradationReaction()
        react.ref_name = prot.ref_name+"_Degradation"
        react.name = react.ref_name.replace("_"," ")
        react.process = "Protein Degradation"
        react.annotations = ['Protein Degradation Reaction']
        react.molecules.append(react_element(mol,prot.localization,-1))
        comp = Counter(prot.composition)
        proteases = []
        for m,s in comp.items():
            factor = 1
            if m[0] == '-':
                m = m[1:]
                factor=-1
            react.molecules.append(react_element(m,'c',s*factor))
        n_prost = 0
        if hasattr(prot,"prosthetic_groups"):
            n_prost = sum([m['stoichiometry'] for m in prot.prosthetic_groups])
        react.molecules.append(react_element("H2O",'c',-len(prot.composition)+1+n_prost))
        if mol in self.created_molecules.keys():    
            if self.created_molecules[mol] == "Proteolysis Tagged Peptide":
                for enz in self.protein_degradation_process.get_signaled_degradation_enzymes():
                    react.enzymes.append(enz)
                    proteases.append(enz['molecule'])
        elif prot.localization == 'c':
            for enz in self.protein_degradation_process.get_intracelullar_degradation_enzymes():
                react.enzymes.append(enz)
                proteases.append(enz['molecule'])
        elif prot.localization == 'm':
            for enz in self.protein_degradation_process.get_membrane_degradation_enzymes():
                react.enzymes.append(enz)
                proteases.append(enz['molecule'])
        for enz in self.protein_degradation_process.get_peptidases():
                react.enzymes.append(enz)
        energy = 0
        for enz in proteases:
            energy += int(len(prot.composition)/self.protein_degradation_process.mean_cleavage_length[enz]* self.protein_degradation_process.mean_cleavage_cost[enz])        
        for m in self.protein_degradation_process.get_energy_molecules():
            m['stoichiometry'] *= energy
            react.molecules.append(m)
        return react
        
        
    def make_biosynthesis_for_molecule(self,mol):
        if self.molecules[mol].mol_type == 'ProteinMonomer':
            return self.make_maturation_reaction_for_protein(mol)
        if self.molecules[mol].mol_type == 'ImatureProteinMonomer':
            return self.make_translation_elongation_reaction_for_imature_protein(mol)
        if self.molecules[mol].mol_type == 'ProteinComplex':
            if mol in self.created_molecules.keys():
                if self.created_molecules[mol] in ['Protein DNA Complex']:
                    return self.make_protein_dna_bind_reaction(mol)
                if self.created_molecules[mol] in ['Condensed DNA Complex']:
                    react = self.make_protein_dna_bind_reaction(mol)
                    react.is_reversible = False
                    return react
                if self.created_molecules[mol] in ['DNA DnaA Complex']:
                    return self.make_dnaA_polymerization_reaction(mol)
        if self.molecules[mol].mol_type == 'TranslationComplex':
            return self.make_translation_complex_biosynthesis(mol)
        if self.molecules[mol].mol_type in ['RNA','ImatureRNA']:
            if len(self.molecules[mol].modification_reactions):
                return self.make_maturation_reaction_for_rna(mol)
            if self.transcription_process.tu_transcript_sufix in mol or not self.transcription_units[self.molecules[mol].transcription_unit].is_cleaved:
                return self.make_rna_synthesis_reaction(mol)
            if self.molecules[mol].transcription_unit+"_"+self.transcription_process.tu_transcript_sufix+"_Processing" not in self.reactions.keys():
                return self.make_rna_cleavage_reaction(mol)
            else:
                return self.reactions[self.molecules[mol].transcription_unit+"_"+
                                      self.transcription_process.tu_transcript_sufix+"_Processing"]            
        if self.molecules[mol].mol_type in ['TranscriptionComplex']:
            if self.created_molecules[mol] == 'Transcribing Complex':
                return self.make_transcription_elongation_reaction(mol)
            if self.created_molecules[mol] == 'Transcription Complex':
                return self.make_transcription_complex_formation_reaction(mol)
        if self.molecules[mol].mol_type in ['DNA']:
            if self.created_molecules[mol] == 'DNA':
                if mol+"_Replication_Reaction" in self.reactions.keys():
                    return self.reactions[mol+"_Replication_Reaction"] 
                return self.make_dna_replication_reaction(mol)
        if self.molecules[mol].mol_type in ['ReplicationComplex']:
            if self.created_molecules[mol] == 'Replicating Complex':
                reg = self.chromosome_regions[mol.replace("_"+self.replication_process.active_replication_complex_sufix,"")]
                if reg.replication_position == 'End':
                    previous_region = [reg.next_region,reg.previous_region]
                if reg.replication_position == 'Left':
                    previous_region = [reg.next_region]
                elif reg.replication_position == 'Right':
                    previous_region = [reg.previous_region]                
                dna_names = [pr+"_"+self.transcription_process.dna_sufix for pr in previous_region]
                reacts = []
                for dna in dna_names:
                    if dna+"_Replication_Reaction" in self.reactions.keys():
                        reacts.append(self.reactions[dna+"_Replication_Reaction"])
                    else:
                        reacts.append(self.make_dna_replication_reaction(dna))
                return reacts
            else:
                reg = self.chromosome_regions[mol.replace("_"+self.replication_process.replication_complex_sufix,"")]
                if reg.ref_name+"_Replication_Initiation" in self.reactions.keys():
                    return self.reactions[chrm+"_Replication_Initiation"]
                return self.make_replication_complex_formation_reaction(mol)
        return None
    

    
    def make_dnaA_polymerization_reaction(self,mol):
        if self.replication_process.active_dnaA['molecule'] in mol and len(mol)>len(self.replication_process.active_dnaA['molecule']):
            return self.make_protein_dna_bind_reaction(mol)
        else:
            reacts = []
            
            dnaA_complexes = [m['molecule'] for m in self.replication_process.dnaA_complex]
            inactive_dnaA_complexes = [m['molecule'] for m in self.replication_process.inactive_dnaA_complex]
            txt = mol.split("_"+self.transcription_process.dna_sufix+"_")
            chrm_region = txt[0]
            dnaA = txt[1]
            is_inactive = dnaA in inactive_dnaA_complexes            
            
            if is_inactive:
                previous_dnaA = dnaA_complexes[inactive_dnaA_complexes.index(dnaA)]
            else:
                previous_dnaA = dnaA_complexes[dnaA_complexes.index(dnaA)-1]
            
            react = ComplexBiosynthesisReaction()
            react.ref_name = mol+"_Complex_Formation"
            react.name = mol+" Complex Formation"
            react.process = "Replication"
            react.annotations = ['DnaA Complex Formation']
            previous_dnaA_complex = mol.replace(dnaA,previous_dnaA)
            react.molecules.append(react_element(mol,'c'))
            react.molecules.append(react_element(previous_dnaA_complex,'c',-1))
            if is_inactive:
                react.molecules.append(self.replication_process.get_inactive_dnaA())
            else:
                react.molecules.append(self.replication_process.get_active_dnaA())
            react.molecules[-1]['stoichiometry'] *= -1
            self.created_molecules[previous_dnaA_complex] = "DNA DnaA Complex"            
            reacts.append(react)
            
            react = ComplexBiosynthesisReaction()
            react.ref_name = mol+"_Complex_Dissociation"
            react.name = mol+" Complex Dissociation"
            react.process = "Replication"
            react.annotations = ['DnaA Complex Dissociation']
            react.molecules.append(react_element(mol,'c',-1))
            react.molecules.append(react_element(dnaA,'c'))
            react.molecules.append(react_element(chrm_region+"_"+self.transcription_process.dna_sufix,'c'))
            reacts.append(react)
            return reacts
            
    
    def make_replication_complex_formation_reaction(self,mol):
        chrm_region_name = mol.replace("_"+self.replication_process.replication_complex_sufix,"")
        chrm_region = self.chromosome_regions[chrm_region_name]
        chrm = chrm_region.chromosome
        react = ComplexBiosynthesisReaction()
        react.ref_name = chrm+"_"+"Replication_Initiation"
        react.name = react.ref_name.replace("_"," ")
        react.process = "Replication"
        react.annotations = ['Replication Initiation Reaction']
        react.annotations.append(chrm)
        react.is_reversible = False
        region_0 = chrm+"_region_0"
        region_last = chrm+"_region_"+str(self.chromosomes[chrm].n_regions-1)
        complex_0 = region_0+"_"+self.replication_process.replication_complex_sufix
        complex_last = region_last+"_"+self.replication_process.replication_complex_sufix
        react.molecules.append(react_element(complex_0,'c'))
        react.molecules.append(react_element(complex_last,'c'))
        dna_0 = region_0+"_"+self.transcription_process.dna_sufix
        dna_last = region_last+"_"+self.transcription_process.dna_sufix
        react.molecules.append(react_element(dna_0,'c',-1))
        react.molecules.append(react_element(dna_last,'c',-1))
        self.created_molecules[dna_0] = "DNA"
        self.created_molecules[dna_last] = "DNA"
        for prot in self.replication_process.get_replication_complex_composition():
            react.molecules.append(prot)
            react.molecules[-1]['stoichiometry'] *= -2
        react.molecules.append(self.replication_process.get_last_dnaA_complex())
        react.molecules[-1]['stoichiometry'] *= len(self.replication_process.dnaA_replication_boxes)
        dnaA_boxes = Counter([self.chromosomes[chrm].features[box]['start'] for box in self.replication_process.dnaA_replication_boxes])
        for box,sto in dnaA_boxes.items():
            reg = self.get_chromosome_region_from_coordinate(chrm,box)
            react.molecules.append(react_element(reg+"_"+self.transcription_process.dna_sufix +"_"+self.replication_process.dnaA_complex[-1]['molecule'],'c',-sto))
            self.created_molecules[react.molecules[-1]['molecule']] = "DNA DnaA Complex"
            react.molecules.append(react_element(reg+"_"+self.transcription_process.dna_sufix,'c'))
    
        return react
    
    def make_dna_replication_reaction(self,mol):
        chrm_region_name = mol.replace("_"+self.transcription_process.dna_sufix,"")
        chrm_region = self.chromosome_regions[chrm_region_name]
        react = DNAReplicationReaction()
        react.ref_name = mol+"_"+"Replication_Reaction"
        react.name = react.ref_name.replace("_"," ")
        react.process = "Replication"
        react.annotations = ['Replication Reaction']
        react.chromosome_region = chrm_region_name
        dna = mol
        if chrm_region.replication_position != 'End':
            react.molecules.append(react_element(dna,'c',2))
            
            if chrm_region.is_replication_origin or self.chromosome_regions[chrm_region.next_region].is_replication_origin:
                react.molecules.append(react_element(chrm_region_name+ "_"+ self.replication_process.replication_complex_sufix, 'c',-1))
                self.created_molecules[chrm_region_name+"_"+self.replication_process.replication_complex_sufix] = "Replication Complex"
            else:
                react.molecules.append(react_element(chrm_region_name+ "_"+ self.replication_process.active_replication_complex_sufix, 'c',-1))
                self.created_molecules[chrm_region_name+"_"+self.replication_process.active_replication_complex_sufix] = "Replicating Complex"
            
            next_region = None
            if chrm_region.replication_position == 'Right':
                next_region = chrm_region.next_region
            elif chrm_region.replication_position == 'Left':
                next_region = chrm_region.previous_region
            react.next_chromosome_region = next_region
            react.molecules.append(react_element(next_region+ "_"+self.replication_process.active_replication_complex_sufix,'c'))
            self.created_molecules[next_region+"_"+self.replication_process.active_replication_complex_sufix] = "Replicating Complex"
            if self.chromosome_regions[next_region].is_replication_terminus:
                react.auxiliar_molecules.append(react_element(next_region+"_"+self.transcription_process.dna_sufix,'c'))
            else:
                react.molecules.append(react_element(next_region+"_"+self.transcription_process.dna_sufix,'c',-1))
            self.created_molecules[next_region+"_"+self.transcription_process.dna_sufix] = "DNA"
        else:
            react.molecules.append(react_element(dna,'c'))
            react.molecules.append(react_element(chrm_region_name+ "_"+ self.replication_process.active_replication_complex_sufix, 'c',-2))
            for prot in self.replication_process.get_replication_complex_composition():
                prot['stoichiometry'] *= 2
                react.molecules.append(prot)
        for prot in self.replication_process.get_dna_repair_enzymes():
            react.auxiliar_molecules.append(prot)        
        ntp = Counter(self.replication_process.get_sequence_composition(chrm_region.sequence,n_phosphates=3)+self.replication_process.get_sequence_composition(complimentary_sequence(chrm_region.sequence),n_phosphates=3))
        for nn,val in ntp.items():
            react.molecules.append(react_element(nn,'c',-val))
        react.molecules.append(react_element('PI','c',4*chrm_region.length))

        
        react.molecules.append(react_element('H2O','c',-2*(chrm_region.length)))
        react.molecules.append(react_element('H','c',2*(chrm_region.length)))
        for m in self.replication_process.get_energy_molecules():
            react.molecules.append(m)
        return react
    
    
    def make_protein_dna_bind_reaction(self,mol):
        txt = mol.split("_"+self.transcription_process.dna_sufix+"_")
        chrm_region_name = txt[0]
        dna_name = txt[0]+"_"+self.transcription_process.dna_sufix
        prot_name = txt[1]
        react = ProteinDNABinding()
        react.ref_name = prot_name+"_"+dna_name+"_"+"Binding_Reaction"
        react.name = react.ref_name.replace("_"," ")
        react.process = "Protein-DNA Interaction "
        react.annotations = ['Protein DNA Binding Reaction']
        react.molecules.append(react_element(mol,'c'))
        react.molecules.append(react_element(dna_name,'c',-1))
        react.molecules.append(react_element(prot_name,'c',-1))
        react.chromosome_region = chrm_region_name
        react.binder = prot_name
        return react
    
    def make_transcription_complex_formation_reaction(self,mol):
        txt = mol.replace("_"+self.transcription_process.transcription_complex_sufix,"").split("_"+self.transcription_process.dna_sufix+"_")
        chrm_region_name = txt[0]
        tu_name = txt[1]
        chrm_region = self.chromosome_regions[chrm_region_name]
        tu = self.transcription_units[tu_name]
        dna = chrm_region_name+"_"+self.transcription_process.dna_sufix
        reacts = []
        for b in tu.transcription_regulators:
            react = TranscriptionComplexFormation()
            react.ref_name = tu.ref_name+"_"+"Transcription_Complex_Formation"
            react.name = react.ref_name.replace("_"," ")
            react.process = "Transcription"
            react.annotations = ['RNA Transcription Complex Formation Reaction']     
            react.chromosome_region = chrm_region_name
            react.transcription_unit = tu_name
            react.molecules.append(react_element(mol,'c'))
            react.molecules.append(self.transcription_process.get_rna_polymerase_holoenzyme())
            react.molecules[-1]['stoichiometry'] *= -1
            if b == None:
                react.molecules.append(react_element(dna,'c',-1))
            else:
                react.molecules.append(react_element(dna+"_"+b,'c',-1))
                self.created_molecules[dna+"_"+b] = "Protein DNA Complex"
                react.molecules.append(react_element(b,'c'))
                react.transcription_factor = tu.transcription_regulators[b]
                react.ref_name = b+"_"+react.ref_name
            reacts.append(react)
        return reacts
                
    
    def make_transcription_elongation_reaction(self,mol):
        txt = mol.replace("_"+self.transcription_process.active_transcription_complex_sufix,"").split("_"+self.transcription_process.dna_sufix+"_")
        chrm_region_name = txt[0]
        tu_name = txt[1]
        chrm_region = self.chromosome_regions[chrm_region_name]
        tu = self.transcription_units[tu_name]
        if tu.direction == "Reverse":
            previous_region = chrm_region.next_region
            next_region = chrm_region.previous_region  
        else:
            previous_region = chrm_region.previous_region
            next_region = chrm_region.next_region
        transcription_start = False
        if tu.ref_name in self.chromosome_regions[previous_region].transcription_unit_start:
            transcription_complex = previous_region+"_"+self.transcription_process.dna_sufix+"_"+tu.ref_name+"_"+self.transcription_process.transcription_complex_sufix
            self.created_molecules[transcription_complex] = "Transcription Complex"
            transcription_start = True
        else:
            transcription_complex = previous_region+"_"+self.transcription_process.dna_sufix+"_"+tu.ref_name+"_"+self.transcription_process.active_transcription_complex_sufix
            self.created_molecules[transcription_complex] = "Transcribing Complex"
        react = RNASynthesisReaction()
        react.ref_name = previous_region+"_"+tu.ref_name+"_"+"Transcription_Elongation"
        react.name = tu.ref_name+" "+"Transcription Elongation"
        react.process = "Transcription"
        react.annotations = ['RNA Transcription Elongation Reaction']      
        react.chromosome_region = chrm_region_name
        react.next_chromosome_region = next_region
        react.molecules.append(react_element(mol,'c'))
        react.molecules.append(react_element(transcription_complex,'c',-1))
        react.molecules.append(react_element(previous_region+"_"+self.transcription_process.dna_sufix, 'c'))
        self.created_molecules[previous_region+"_"+self.transcription_process.dna_sufix] = "DNA"
        if transcription_start:
            for ef in self.transcription_process.get_transcription_elongation_factors():
                react.molecules.append(ef)
                react.molecules[-1]['stoichiometry'] *= -1
            react.molecules.append(self.transcription_process.get_sigma_factor())
        react.molecules.append(react_element(chrm_region.ref_name+"_"+self.transcription_process.dna_sufix, 'c',-1))
        seq = self.get_dna_subsequence(*tu.get_coordinates_inside_chromosome_region(self.chromosome_regions[previous_region]))
        ntp = Counter(self.transcription_process.get_rna_composition( self.transcription_process.transcribe_gene(seq),n_phosphates=3))
        for nn,val in ntp.items():
            react.molecules.append(react_element(nn,'c',-val))
        react.molecules.append(react_element('PI','c',2*len(seq)))
        react.molecules.append(react_element('H2O','c',-len(seq)))
        react.molecules.append(react_element('H','c',len(seq)))
            
        return react
        
        
    
    def make_rna_synthesis_reaction(self,mol):
        rna = self.molecules[mol]
        tu = self.transcription_units[rna.transcription_unit]
        react = RNASynthesisReaction()
        react.ref_name = tu.ref_name+"_"+"Transcription_End"
        react.name = tu.ref_name+" "+"Transcription End"
        react.process = "Transcription"
        react.annotations = ['RNA Transcription Reaction']
        react.transcript = rna.ref_name
        react.molecules.append(react_element(rna.ref_name,'c'))
        react.molecules.append(self.transcription_process.get_rna_polymerase())
        for rf in self.transcription_process.transcription_release_factors:
            react.auxiliar_molecules.append(rf)
        if tu.direction == "Forward":
            chrm_region = self.get_chromosome_region_from_coordinate(tu.chromosome,tu.transcription_end)
        elif tu.direction == "Reverse":
            chrm_region = self.get_chromosome_region_from_coordinate(tu.chromosome,tu.transcription_start)
        react.chromosome_region = chrm_region
        react.molecules.append(react_element(chrm_region+"_"+self.transcription_process.dna_sufix,'c'))
        self.created_molecules[chrm_region+"_"+self.transcription_process.dna_sufix] = "DNA"
        if tu.ref_name not in self.chromosome_regions[chrm_region].transcription_unit_start:            
            transcription_complex = chrm_region+"_"+self.transcription_process.dna_sufix+"_"+tu.ref_name+"_"+ self.transcription_process.active_transcription_complex_sufix        
            react.molecules.append(react_element(transcription_complex,'c',-1))
            self.created_molecules[transcription_complex] = "Transcribing Complex"
            for ef in self.transcription_process.get_transcription_elongation_factors():
                react.molecules.append(ef)
        else:            
            transcription_complex = chrm_region+"_"+self.transcription_process.dna_sufix + "_" + tu.ref_name + "_" + self.transcription_process.transcription_complex_sufix        
            react.molecules.append(react_element(transcription_complex,'c',-1))
            self.created_molecules[transcription_complex] = "Transcription Complex"
            for ef in self.transcription_process.get_transcription_elongation_factors():
                react.auxiliar_molecules.append(ef)
            react.molecules.append(self.transcription_process.get_sigma_factor())
        seq = self.get_dna_subsequence(*tu.get_coordinates_inside_chromosome_region(self.chromosome_regions[chrm_region]))
        ntp = Counter(self.transcription_process.get_rna_composition( self.transcription_process.transcribe_gene(seq),n_phosphates=3))
        for nn,val in ntp.items():
            react.molecules.append(react_element(nn,'c',-val))
        react.molecules.append(react_element('PI','c',2*len(seq)))
        react.molecules.append(react_element('H2O','c',-len(seq)))
        react.molecules.append(react_element('H','c',len(seq)))
        return react        
        
    
    def make_rna_cleavage_reaction(self,mol):
        rna = self.molecules[mol]
        tu = self.transcription_units[rna.transcription_unit]
        react = RNAProcessingReaction()
        tu_transcript = tu.ref_name+"_"+self.transcription_process.tu_transcript_sufix
        react.ref_name = tu_transcript+"_"+"Processing"
        react.name = tu_transcript+" "+"Processing Reaction"
        react.process = "RNA Processing"
        react.annotations = ['RNA Processing Reaction']
        react.transcript = tu_transcript
        react.rnas = tu.genes.copy()
        react.molecules.append(react_element(tu_transcript,'c',-1))
        self.created_molecules[tu_transcript] = 'Polycistronic Transcript'
        n_cleavages = 0
        for i in range(len(react.rnas)):
            rna_name = react.rnas[i]
            gene = self.genes[rna_name]
            gene_coordinates = gene.get_gene_coordinates()
            if len(gene.modification_reactions):
                rna_name = self.transcription_process.imature_rna_prefix+"_"+rna_name
            react.molecules.append(react_element(rna_name,'c'))
            if i<len(react.rnas)-1:
                n_cleavages += 1
                next_gene_coordinates = self.genes[react.rnas[i+1]].get_gene_coordinates()
                start = gene_coordinates[1]+gene_coordinates[2]
                length = next_gene_coordinates[1]-start
                inter_gene_coordinates = (gene_coordinates[0],start,length,gene_coordinates[3])
                if length>0:
                    n_cleavages += 1
                    oligo_seq = self.transcription_process.transcribe_gene(self.get_dna_subsequence(*inter_gene_coordinates))
                    if length == 1:
                        react.molecules.append(react_element(self.transcription_process.get_rna_composition(oligo_seq)[0],'c'))
                        react.molecules.append(react_element('H2O','c',-1))
                        react.molecules.append(react_element('H','c'))
                    else:
                        oligo_name = "RNA_"+oligo_seq
                        if oligo_name not in self.molecules.keys():
                            self.created_molecules[oligo_name] = "Unknown Function RNA"
                            oligo = self.make_rna(oligo_name)
                            oligo.biosynthesis_reaction =  react.ref_name
                            self.add_molecule(oligo)
                        react.molecules.append(react_element(oligo_name,'c'))
        for rna_type, enzymes in self.transcription_process.get_rna_cleavage_enzymes().items():
            if rna_type in rna.annotations:
                for enz in enzymes:
                    react.enzymes.append(enz)
                if self.transcription_process.cost_per_cleavage[rna_type]:
                    for m in self.transcription_process.get_cleavage_energy_molecules():
                        m['stoichiometry'] *= int(n_cleavages*self.transcription_process.cost_per_cleavage[rna_type])
                        react.molecules.append(m)
        
        return react        
        
        
    def make_maturation_reaction_for_rna(self,mol):      
        rna = self.molecules[mol]
        gene = self.genes[rna.gene]
        gene_seq = self.get_dna_subsequence(*gene.get_gene_coordinates())
        rna.composition = self.transcription_process.get_rna_composition(self.transcription_process.transcribe_gene(gene_seq))
        react = RNAMaturationReaction()
        imature_rna = react_element(self.transcription_process.imature_rna_prefix+"_"+rna.ref_name,'c',-1)
        self.created_molecules[imature_rna['molecule']] = 'Imature RNA'
        mature_rna = react_element(rna.ref_name,'c')
        react.ref_name = rna.ref_name+"_Maturation"
        react.name = rna.name+" Maturation"
        react.process = "RNA Maturation"
        react.annotations = ['RNA Maturation Reaction']
        react.imature_rna = imature_rna['molecule']
        react.mature_rna = mature_rna['molecule']
        react.molecules.append(imature_rna)
        react.molecules.append(mature_rna)
        for mod_react in sorted(rna.modification_reactions):
            for ann in self.transcription_process.get_modification_annotations(mod_react):
                if ann not in react.annotations: react.annotations.append(ann)
            if self.transcription_process.modification_reactions[mod_react].del_base != None:
                self.transcription_process.apply_modification_reaction(rna.composition,mod_react)
            for enzyme in self.transcription_process.get_modification_enzymes(mod_react):
                if enzyme['molecule'] not in [enz['molecule'] for enz in react.enzymes]:
                    react.enzymes.append(enzyme)
            for auxiliar_mol in self.transcription_process.get_auxiliar_molecules(mod_react):
                react.auxiliar_molecules.append(auxiliar_mol)
            molecules = self.transcription_process.get_modification_molecules(mod_react)
            for m in molecules:
                react.molecules.append(m)
        return react
        
    def make_maturation_reaction_for_protein(self,mol):      
        prot = self.molecules[mol]
        gene = self.genes[prot.gene]
        gene_seq = self.get_dna_subsequence(*gene.get_gene_coordinates())
        prot.composition = self.translation_process.translate_gene(gene_seq)
        react = ProteinMaturationReaction()
        imature_protein = react_element(self.translation_process.imature_protein_prefix+"_"+prot.ref_name,'c',-1)
        self.created_molecules[imature_protein['molecule']] = 'Imature Protein'
        mature_protein = react_element(prot.ref_name,prot.localization)
        react.ref_name = prot.ref_name+"_Maturation"
        react.name = prot.name+" Maturation"
        react.process = "Protein Maturation"
        react.annotations = ['Protein Maturation Reaction']
        react.imature_protein = imature_protein['molecule']
        react.mature_protein = mature_protein['molecule']
        react.molecules.append(imature_protein)
        react.molecules.append(mature_protein)
        for mod_react in prot.modification_reactions:
            for ann in self.translation_process.get_modification_annotations(mod_react):
                if ann not in react.annotations: react.annotations.append(ann)
            if self.translation_process.modification_reactions[mod_react].del_aa != None:
                self.translation_process.apply_modification_reaction(prot.composition,mod_react)
            for enzyme in self.translation_process.get_modification_enzymes(mod_react):
                if enzyme['molecule'] not in [enz['molecule'] for enz in react.enzymes]:
                    react.enzymes.append(enzyme)
            for auxiliar_mol in self.translation_process.get_auxiliar_molecules(mod_react):
                react.auxiliar_molecules.append(auxiliar_mol)
            molecules = self.translation_process.get_modification_molecules(mod_react)
            for m in molecules:
                react.molecules.append(m)
        for m in prot.prosthetic_groups:
            prtg = m.copy()
            for i in range(int(prtg['stoichiometry'])):
                prot.composition.append(prtg['molecule'])
            prtg['stoichiometry'] *= -1
            react.molecules.append(prtg)            
        if prot.localization in ['e','m']:
            for p in self.translation_process.get_translocation_proteins():
                react.auxiliar_molecules.append(p)
            for enz in self.translation_process.get_translocation_enzymes():
                react.enzymes.append(enz)
            for m in self.translation_process.get_elongation_energy_molecules():
                m['stoichiometry'] *= self.translation_process.get_gtp_cost_for_translocation()
                react.molecules.append(m)
        if prot.signal_sequence!=None:
            react.molecules.append(react_element('H2O','c',-1))
            prot.is_secreted = prot.localization=='e'
            pep = self.make_signal_peptide_from_protein(prot.ref_name)
            pep.biosynthesis_reaction = react.ref_name
            prot.signal_sequence['composition'] = pep.composition
            if prot.signal_sequence['type'] == 'N-Terminal':
                prot.composition = prot.composition[int(prot.signal_sequence['length']):]
            else:
                prot.composition = prot.composition[:-int(prot.signal_sequence['length'])]
            signal_peptide = react_element(pep.ref_name,'e')
            react.molecules.append(signal_peptide)
            self.created_molecules[signal_peptide['molecule']] = 'Signal Peptide'
            self.add_molecule(pep)
            
        return react
    
    def make_translation_elongation_reaction_for_imature_protein(self,mol_name):
        improt = self.molecules[mol_name]
        prot_name = mol_name.replace(self.translation_process.imature_protein_prefix+"_","")
        prot = self.molecules[prot_name]
        gene = self.genes[prot.gene]
        tu = self.transcription_units[gene.transcription_unit]
        gene_seq = self.get_dna_subsequence(*gene.get_gene_coordinates())
        trans = self.translation_process.translate_gene(gene_seq,trans_out='aatrna')
        if not len(trans): print("Length error:",gene.ref_name)
        aa_trna = Counter(trans)
        trans = self.translation_process.translate_gene(gene_seq,trans_out='trna')
        trna = Counter(trans)
        react = TranslationElongationReaction()
        react.ref_name = prot.ref_name+"_Translation_Elongation"
        react.name = prot.name+" Translation Elongation"
        react.process = "Translation"
        react.annotations = ['Translation Elongation Reaction']
        react.is_reversible = False
        react.imature_protein = improt.ref_name
        react.molecules.append(react_element(improt.ref_name,'c'))
        if tu.is_polycistronic:
            rna = tu.ref_name+"_"+self.transcription_process.tu_transcript_sufix
            self.created_molecules[rna] = 'Poly mRNA'
        else:
            rna = prot.gene
            self.created_molecules[rna] = 'mRNA'
        react.molecules.append(react_element(rna,'c'))        
        react.translation_complex = self.translation_process.translation_complex_prefix+'_'+prot.ref_name
        self.created_molecules[react.translation_complex] = 'Translation Complex'
        react.molecules.append(react_element(react.translation_complex,'c',-1))
        for aat,sto in aa_trna.items():
            react.molecules.append(react_element(aat,'c',-sto))
        for t,sto in trna.items():
            react.molecules.append(react_element(t,'c',sto))
        for m in self.translation_process.get_elongation_auxiliaries():
            react.auxiliar_molecules.append(m)
        for m in self.translation_process.get_translation_complex_constitution():
            react.molecules.append(m)
        for m in prot.chaperones:
            react.auxiliar_molecules.append(react_element(m,'c'))
        for m in self.translation_process.get_elongation_energy_molecules():
            m['stoichiometry'] *= len(trans)
            react.molecules.append(m)
        react.molecules.append(react_element('H2O','c',len(trans)-1))
        return react
    
    def make_translation_complex_biosynthesis(self,mol_name):
        prot_name = mol_name.replace(self.translation_process.translation_complex_prefix+"_","")
        prot = self.molecules[prot_name]
        react = ComplexBiosynthesisReaction()
        react.ref_name = prot.ref_name+"_Translation_Complex_Formation"
        react.name = prot.name+" Translation Complex Formation"
        react.process = "Translation"
        react.annotations = ['Translation Complex Formation']
        react.is_reversible = False
        tu = self.transcription_units[self.genes[prot.gene].transcription_unit]
        if tu.is_polycistronic:
            rna = tu.ref_name+"_"+self.transcription_process.tu_transcript_sufix
            self.created_molecules[rna] = 'Poly mRNA'
        else:
            rna = prot.gene
            self.created_molecules[rna] = 'mRNA'
        react.molecules.append(react_element(rna,'c',-1))
        for m in self.translation_process.get_translation_complex_constitution():
            m['stoichiometry'] *= -1
            react.molecules.append(m)
        for m in self.translation_process.get_initiation_auxiliaries():
            react.auxiliar_molecules.append(m)
        react.molecules.append(react_element(mol_name,'c'))
        return react
        
    def make_molecule_from_name(self,mol):
        if mol not in self.created_molecules.keys():
            return None
        if self.created_molecules[mol] == 'Imature Protein':
            return self.make_imature_protein(mol)
        if self.created_molecules[mol] == 'Imature RNA':
            return self.make_imature_rna(mol)
        if self.created_molecules[mol] == 'Translation Complex':
            return self.make_translation_complex(mol)
        if self.created_molecules[mol] in ['Poly mRNA','mRNA','Polycistronic Transcript']:
            if mol not in self.molecules.keys():
                return self.make_rna(mol)
        if self.created_molecules[mol] in ["Transcribing Complex","Transcription Complex"]:
            return self.make_transcription_complex(mol)
        if self.created_molecules[mol] == "DNA":
            return self.make_dna(mol)
        if self.created_molecules[mol] in ["Protein DNA Complex","DNA DnaA Complex"]:
            return self.make_dna_protein_complex(mol)
        if self.created_molecules[mol] in ["Replicating Complex","Replication Complex"]:
            return self.make_replication_complex(mol)
        return False
            
    
    def make_signal_peptide_from_protein(self,prot_name):
        prot = self.molecules[prot_name]
        pep = Peptide()
        pep.annotations.append('Signal Peptide')
        pep.annotations.append(prot_name+' Signal Peptide')
        pep.gene = prot.gene
        if prot.signal_sequence['type'] == 'N-Terminal':
            pep.composition = prot.composition[:int(prot.signal_sequence['length'])]
        else:
            pep.composition = prot.composition[-int(prot.signal_sequence['length']):]
        pep.ref_name = self.translation_process.peptide_prefix+"_"+"".join([self.translation_process.aminoacid_letter[aa] for aa in pep.composition])
        pep.name = pep.ref_name.replace("_"," ")
        pep.is_secreted = True
        pep.localization = 'e'
        return pep
        
    def make_imature_protein(self,mol_name):
        prot_name = mol_name.replace(self.translation_process.imature_protein_prefix+"_","")
        prot = self.molecules[prot_name]
        improt = ImatureProteinMonomer(prot,self.translation_process.imature_protein_prefix)
        gene = self.genes[prot.gene]
        gene_seq = self.get_dna_subsequence(*gene.get_gene_coordinates())
        improt.composition = self.translation_process.translate_gene(gene_seq)
        return improt
    
    def make_imature_rna(self,mol_name):
        rna_name = mol_name.replace(self.transcription_process.imature_rna_prefix+"_","")
        rna = self.molecules[rna_name]
        imrna = ImatureRNA(rna,self.transcription_process.imature_rna_prefix)
        gene = self.genes[rna.gene]
        gene_seq = self.get_dna_subsequence(*gene.get_gene_coordinates())
        imrna.composition = self.transcription_process.get_rna_composition( self.transcription_process.transcribe_gene(gene_seq))
        return imrna
    
    def make_rna(self,mol_name):
        rna = RNA()
        rna.ref_name = mol_name
        rna.name = mol_name.replace("_"," ")
        rna.annotations.append(self.created_molecules[mol_name])
        if self.created_molecules[mol_name] == "Unknown Function RNA":
            seq = mol_name.replace("RNA_","")
            rna.composition = self.transcription_process.get_rna_composition(seq)
        else:
            if self.transcription_process.tu_transcript_sufix in mol_name:
                tu_name = mol_name.replace("_"+self.transcription_process.tu_transcript_sufix,"")
            else:
                tu_name = self.genes[mol_name].transcription_unit
            rna.transcription_unit = tu_name
            tu = self.transcription_units[tu_name]
            gene_seq = self.get_dna_subsequence(*tu.get_coordinates())
            if tu.ref_name not in rna.ref_name:
                gene_seq = self.get_dna_subsequence(self.genes[mol_name].get_gene_coordinates())        
            rna.composition = self.transcription_process.get_rna_composition( self.transcription_process.transcribe_gene(gene_seq),n_phosphates=1)
        return rna
    
    def make_translation_complex(self,mol_name):
        prot_name = mol_name.replace(self.translation_process.translation_complex_prefix+"_","")
        prot = self.molecules[prot_name]
        mol = TranslationComplex()
        mol.ref_name = mol_name
        mol.name = prot.name+" "+self.translation_process.translation_complex_prefix
        tu = self.transcription_units[self.genes[prot.gene].transcription_unit]
        if tu.is_polycistronic:
            rna = tu.ref_name+"_"+self.transcription_process.tu_transcript_sufix
        else:
            rna = prot.gene
        mol.composition.append(rna)
        for m in self.translation_process.get_translation_complex_constitution():
            mol.composition.append(m['molecule'])
        return mol
    
    
    def make_transcription_complex(self,mol_name):
        txt = mol_name.replace("_"+self.transcription_process.transcription_complex_sufix,"").replace("_"+self.transcription_process.active_transcription_complex_sufix,"").split("_"+self.transcription_process.dna_sufix+"_")
        chrm_region_name = txt[0]
        tu_name = txt[1]
        chrm_region = self.chromosome_regions[chrm_region_name]
        tu = self.transcription_units[tu_name]
        mol = TranscriptionComplex()
        mol.ref_name = mol_name
        mol.name = mol.ref_name.replace("_"," ")        
        if self.created_molecules[mol_name] == "Transcribing Complex":
            mol.composition.append(self.transcription_process.get_rna_polymerase()['molecule'])
            for ef in self.transcription_process.get_transcription_elongation_factors():
                mol.composition.append(ef['molecule'])
            mol.transcribed_sequence = ''
            if tu.ref_name not in chrm_region.transcription_unit_start:                
                current_region = chrm_region
                while tu.ref_name not in current_region.transcription_unit_start:
                    current_region = self.chromosome_regions[current_region.get_previous_region(tu.direction)]
                    mol.transcribed_sequence = self.get_dna_subsequence(*tu.get_coordinates_inside_chromosome_region(current_region)) + mol.transcribed_sequence
                    
            if len(mol.transcribed_sequence):
                rna = "RNA_"+self.transcription_process.transcribe_gene(mol.transcribed_sequence)
                self.created_molecules[rna] = "Unknown Function RNA"
                mol.composition.append(rna)     
        else:
            mol.composition.append(self.transcription_process.get_rna_polymerase_holoenzyme()['molecule'])
        dna = chrm_region_name+"_"+self.transcription_process.dna_sufix
        mol.composition.append(dna)
        
            
        return mol
    
    def make_dna(self,mol_name):
        chrm_region = self.chromosome_regions[mol_name.replace("_"+self.transcription_process.dna_sufix,"")]
        mol = DNA()
        mol.ref_name = mol_name
        mol.name = mol.ref_name.replace("_"," ")
        mol.chromosome_region = chrm_region.ref_name
        mol.composition = self.replication_process.get_sequence_composition(chrm_region.sequence,n_phosphates=1)+self.replication_process.get_sequence_composition(complimentary_sequence(chrm_region.sequence),n_phosphates=1)
        return mol
    
    def make_dna_protein_complex(self,mol_name):
        txt = mol_name.split("_"+self.transcription_process.dna_sufix+"_")
        dna_name = txt[0]+"_"+self.transcription_process.dna_sufix
        prot_name = txt[1]
        mol = ProteinComplex()
        mol.ref_name = mol_name
        mol.name = prot_name+" bound to "+dna_name
        mol.composition.append(dna_name)
        mol.composition.append(prot_name)
        return mol
    
    def make_replication_complex(self,mol_name):
        chrm_region_name = mol_name.replace("_"+self.replication_process.replication_complex_sufix,"").replace("_"+self.replication_process.active_replication_complex_sufix,"")
        chrm_region = self.chromosome_regions[chrm_region_name]
        mol = ReplicationComplex()
        mol.ref_name = mol_name
        mol.name = mol.ref_name.replace("_"," ")
        mol.chromosome_region = chrm_region_name
        if not chrm_region.is_replication_terminus:
            dna = chrm_region_name+"_"+self.transcription_process.dna_sufix
            mol.composition.append(dna)
        for prot in self.replication_process.get_replication_complex_composition():
            for i in range(prot['stoichiometry']):
                mol.composition.append(prot['molecule'])
        return mol
    
    
class Chromosome(Jsonable):
    def __init__(self):
        self.ref_name = None
        self.name = None
        self.sequence = None
        self.length = None
        self.features = {}
        self.genes = []
        self.original_database = None
        self.offset = None
        self.n_regions = None
        
    def add_feature(self,feat_name,feat):
        if not issubclass(type(feat),dict):
            print("ERROR: feat must be a 'dict' object. Got:",type(feat))
            return False
        self.features[feat_name] = feat
        return True 
    

class ChromosomeRegion(Jsonable):
    def __init__(self):
        self.ref_name = None
        self.sequence = None
        self.start = None
        self.length = None
        self.end = None
        self.transcription_unit_start = []
        self.transcription_unit_end = []
        self.is_replication_origin = False
        self.is_replication_terminus = False
        self.previous_region = None
        self.next_region = None
        self.chromosome = None
        self.features = []
        self.replication_position = None
        
    def get_next_region(self,direction):
        if direction == "Forward":
            return self.next_region
        elif direction == "Reverse":
            return self.previous_region
        else:
            return False
        
    def get_previous_region(self,direction):
        if direction == "Forward":
            return self.previous_region
        elif direction == "Reverse":
            return self.next_region
        else:
            return False
        
        
class TranscriptionUnit(Jsonable):
    def __init__(self):
        self.ref_name = None
        self.name = None
        self.genes = None
        self.promoter_10_start = None
        self.promoter_10_length = None
        self.promoter_35_start = None
        self.promoter_35_length = None
        self.transcription_start = None
        self.transcription_end = None
        self.transcription_length = None
        self.direction = None
        self.type = None
        self.original_database = None
        self.is_cleaved = False
        self.is_polycistronic = None
        self.chromosome = None
        self.transcription_regulators = {None:1}
        
    def get_coordinates(self):
        return self.chromosome,self.transcription_start,self.transcription_length,self.direction
    
    def get_coordinates_inside_chromosome_region(self,chrm_reg):
        if self.chromosome != chrm_reg.chromosome:
            print("ERROR: Transcription Unit not in "+chrm_reg.chromosome+". TU: "+self.ref_name+" Chromosome: "+self.chromosome)
            return False        
        if self.transcription_start<chrm_reg.start and self.transcription_end>chrm_reg.end:
            return self.chromosome,chrm_reg.start,chrm_reg.length,self.direction
        if self.transcription_start>=chrm_reg.start and self.transcription_end<=chrm_reg.end:
            return self.get_coordinates()
        if self.transcription_start<chrm_reg.start and self.transcription_end<=chrm_reg.end:
            return self.chromosome,chrm_reg.start,self.transcription_end-chrm_reg.start,self.direction
        if self.transcription_start>=chrm_reg.start and self.transcription_end>chrm_reg.end:
            return self.chromosome,self.transcription_start,chrm_reg.end-self.transcription_start+1,self.direction
        
        
        print("ERROR: Transcription Unit not in "+chrm_reg.ref_name+". TU: "+self.ref_name)
        return False
            
    
class Gene(Jsonable):
    def __init__(self):
        self.ref_name = None
        self.name = None
        self.symbol = None
        self.chromosome = None
        self.transcription_unit = None
        self.start = None
        self.length = None
        self.direction = None
        self.original_database = None
        self.modification_reactions = []
        
    def get_gene_coordinates(self):
        return self.chromosome,self.start,self.length,self.direction
    
class Molecule:
    def __init__(self):
        self.ref_name = None
        self.name = None
        self.mw = 0.
        self.volume = 0.
        self.mol_type = None
        self.cross_references = []
        self.annotations = []
        self.original_database = ''
    
    def get_type(self):
        return self.mol_type
    
    def get_node_attr(self):
        attr = {'type':'m'}
        for at,value in self.__dict__.items():
            at = at.replace("_","")
            if type(value) in [str,float,int]:
                if at == "name":
                    at = "usualname"
                if at == "refname":
                    at = "name"
                attr[at] = value
            else:
                attr[at] = json.dumps(value, ensure_ascii=False)        
        return attr

    def get_node(self):
        return self.ref_name,self.get_node_attr()
    
class SmallMolecule(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "SmallMolecule"
        self.formula = None
        self.smiles = None
        self.iupac_name = None
        self.is_hydrophobic = None
        self.charge = None
        self.synonyms = []
        self.log_d = None
        self.log_p = None
        self.deltag_formation = None
        
    def compute_mw_from_formula(self):
        if self.formula == None:
            return False
        self.mw = float(str(molmass.Formula(self.formula).mass))
        return True
        
class ProteinMonomer(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "ProteinMonomer"
        self.synonyms = []
        self.gene = None
        self.biosynthesis_reaction = None
        self.degradation_reaction = None
        self.composition = []
        self.chaperones = []
        self.localization = None
        self.modification_reactions = []
        self.prosthetic_groups = []
        self.signal_sequence = None
        self.is_secreted = False

class Peptide(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "Peptide"
        self.synonyms = []
        self.gene = None
        self.biosynthesis_reaction = None
        self.degradation_reaction = None
        self.composition = []
        self.localization = None
        self.modification_reactions = []
        self.is_secreted = False
        
class ImatureProteinMonomer(ProteinMonomer,Jsonable):
    def __init__(self):
        ProteinMonomer.__init__(self)
        self.mol_type = "ImatureProteinMonomer"
        
    def __init__(self,prot,prefix):
        ProteinMonomer.__init__(self)
        self.ref_name = prefix+'_'+prot.ref_name
        self.name = prefix+' '+prot.name
        self.volume = prot.volume
        self.mol_type = "ImatureProteinMonomer"
        self.cross_references = prot.cross_references
        self.annotations = prot.annotations+['Imature Protein Monomer']
        self.original_database = 'Created by model'
        self.synonyms = prot.synonyms
        self.gene = prot.gene
        self.biosynthesis_reaction = None
        self.composition = []
        self.chaperones = prot.chaperones
        self.localization = 'c'
        self.modification_reactions = []
        self.prosthetic_groups = []
        

class ProteinComplex(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "ProteinComplex"
        self.synonyms = []
        self.biosynthesis_reaction = None
        self.degradation_reaction = None
        self.composition = []
    

class RNA(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "RNA"
        self.synonyms = []
        self.biosynthesis_reaction = None
        self.degradation_reaction = None
        self.composition = []
        self.transcription_unit = None
        self.modification_reactions = []
        
class ImatureRNA(RNA,Jsonable):
    def __init__(self):
        RNA.__init__(self)
        self.mol_type = "ImatureRNA"
        
    def __init__(self,rna,prefix):
        RNA.__init__(self)
        self.ref_name = prefix+'_'+rna.ref_name
        self.name = prefix+' '+rna.name
        self.mw = rna.mw
        self.volume = rna.volume
        self.mol_type = "ImatureRNA"
        self.cross_references = rna.cross_references
        self.annotations = rna.annotations+['Imature RNA']
        self.original_database = 'Created by model'
        self.synonyms = rna.synonyms
        self.gene = rna.gene
        self.transcription_unit = rna.transcription_unit
    

    
class RNAComplex(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "RNAComplex"
        self.synonyms = []
        self.biosynthesis_reaction = None
        self.degradation_reaction = None
        self.composition = []
        
class TranscriptionComplex(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "TranscriptionComplex"
        self.synonyms = []
        self.biosynthesis_reaction = None
        self.degradation_reaction = None
        self.composition = []
        self.chromosome_region = None
        self.transcription_unit = None
        self.transcribed_sequence = None

class TranslationComplex(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "TranslationComplex"
        self.synonyms = []
        self.biosynthesis_reaction = None
        self.degradation_reaction = None
        self.composition = []
    
class Stimulus(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "Stimulus"
        self.synonyms = []
        self.value = None
        self.mw = 0.001
        
class DNA(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "DNA"
        self.synonyms = []
        self.biosynthesis_reaction = None
        self.composition = []
        self.chromosome_region = None
        self.modification_reactions = []

class ReplicationComplex(Molecule,Jsonable):
    def __init__(self):
        Molecule.__init__(self)
        self.mol_type = "ReplicationComplex"
        self.synonyms = []
        self.biosynthesis_reaction = None
        self.composition = []
        self.chromosome_region = None

class Reaction:
    def __init__(self):
        self.ref_name = None
        self.name = None
        self.react_type = None
        self.cross_references = []
        self.annotations = []
        self.original_database = None
        self.process = None
        self.enzymes = []
        self.auxiliar_molecules = []
        self.molecules = []
        self.is_reversible = None
        self.mass_balance = None
        
    def get_type(self):
        return self.react_type
        
    def get_node(self):
        attr = {'type':'r'}
        for at,value in self.__dict__.items():
            if at in ['molecules','auxiliar_molecules','enzymes']: continue
            at = at.replace("_","")
            if type(value) in [str,float,int]:
                if at == "name":
                    at = "usualname"
                if at == "refname":
                    at = "name"
                attr[at] = value
            else:
                attr[at] = json.dumps(value, ensure_ascii=False)        
        return self.ref_name,attr
        
    def get_molecules_and_locations(self):
        mols = []
        for m in self.molecules:
            mols.append((m['molecule'],m['location'], m['location']+"_"+m['molecule']))
        for m in self.auxiliar_molecules:
            mols.append((m['molecule'],m['location'], m['location']+"_"+m['molecule']))
        for m in self.enzymes:
            mols.append((m['molecule'],m['location'], m['location']+"_"+m['molecule']))
        return mols
        
    def get_edges(self,rev=False):
        edges = []
        react_name = self.ref_name
        if rev:
            react_name += "_rev"
        for m in self.molecules:
            mol_name = m['location']+'_'+m['molecule']            
            attr = {'sto':abs(m['stoichiometry'])}
            if (m['stoichiometry']<0 and not rev) or (m['stoichiometry']>0 and rev):
                attr['type']='r'
                edges.append((mol_name,react_name,attr))
            else:
                attr['type']='p'
                edges.append((react_name,mol_name,attr))
        for m in self.enzymes:
            mol_name = m['location']+'_'+m['molecule']
            attr = {'sto':abs(m['stoichiometry']),'type':'m'}
            edges.append((mol_name,react_name,attr))
        for m in self.auxiliar_molecules:
            mol_name = m['location']+'_'+m['molecule']
            attr = {'sto':abs(m['stoichiometry']),'type':'m'}
            edges.append((mol_name,react_name,attr))
        return edges

    
    
class ChemicalReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "ChemicalReaction"
        self.rate_law_forward = None
        self.rate_law_backward = None
        self.parameters_forward = {}
        self.parameters_backward = {}
        self.optimal_ph = None
        self.optimal_temperature = None
        self.activators = {}
        self.inhibitors = {}
        self.pathways = []


class RNADegradationReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "RNADegradationReaction"
        self.is_reversible = False
        
class ProteinDegradationReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "ProteinDegradationReaction"
        self.is_reversible = False
        
class ComplexBiosynthesisReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "ComplexBiosynthesisReaction"
        self.is_reversible = True



class TranslationElongationReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "TranslationElongationReaction"
        self.translation_complex = None
        self.imature_protein = None
        
        
class ProteinMaturationReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "ProteinMaturationReaction"
        self.imature_protein = None
        self.mature_protein = None
        self.is_reversible = False
        
class ProteinModificationReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "ProteinModificationReaction"
        self.del_aa = None
        self.add_aa = None
        self.position = None
        
class RNAModificationReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "RNAModificationReaction"
        self.del_base = None
        self.add_base = None
        self.position = None
        
class RNAMaturationReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "RNAMaturationReaction"
        self.imature_rna = None
        self.mature_rna = None
        self.is_reversible = False   
        
class RNAProcessingReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "RNAProcessingReaction"
        self.transcript = None
        self.rnas = []
        self.is_reversible = False   
        
class RNASynthesisReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "RNASynthesisReaction"
        self.transcript = None
        self.chromosome_region = None
        self.next_chromosome_region = None
        self.is_reversible = False   
        
class TranscriptionComplexFormation(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "TranscriptionComplexFormation"
        self.transcription_unit = None
        self.chromosome_region = None
        self.is_reversible = False
        self.transcription_factor = 1.

class ProteinDNABinding(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "ProteinDNABinding"
        self.chromosome_region = None
        self.is_reversible = True
        self.binder = None

class DNAReplicationReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "DNAReplicationReaction"
        self.chromosome_region = None
        self.next_chromosome_region = None
        self.is_reversible = False
        
class CellDivisionReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "CellDivisionReaction"
        self.is_reversible = False
        self.ref_name = "cellular_division_reaction"
        self.name = "Cellular Division Reaction"
        
class TransportReaction(Reaction,Jsonable):
    def __init__(self):
        Reaction.__init__(self)
        self.react_type = "TransportReaction"
        self.from_compartment = None
        self.to_compartment = None
        

class RNADegradation(Jsonable):
    def __init__(self):
        self.degradation_enzymes = []
        self.aminoacylated_trna_degradation_enzymes = []
        
    def set_degradation_enzymes(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.degradation_enzymes = prot_list
        return True
    
    def set_aminoacylated_trna_degradation_enzymes(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.aminoacylated_trna_degradation_enzymes = prot_list
        return True
    
    def get_degradation_enzymes(self):
        return [e.copy() for e in self.degradation_enzymes]
    
    def get_aminoacylated_trna_degradation_enzymes(self):
        return [e.copy() for e in self.aminoacylated_trna_degradation_enzymes]
        
class ProteinDegradation(Jsonable):
    def __init__(self):
        self.intracellullar_degradation_enzymes = []
        self.membrane_degradation_enzymes = []
        self.signaled_degradation_enzymes = []
        self.peptidases = []
        self.mean_cleavage_cost = {}
        self.mean_cleavage_length = {}
        self.energy_molecules = []
        
    def set_intracellullar_degradation_enzymes(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.intracellullar_degradation_enzymes = prot_list
        return True
    
    def set_membrane_degradation_enzymes(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.membrane_degradation_enzymes = prot_list
        return True
    
    def set_signaled_degradation_enzymes(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.signaled_degradation_enzymes = prot_list
        return True
    
    def set_peptidases(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.peptidases = prot_list
        return True
    
    def get_intracelullar_degradation_enzymes(self):
        return [e.copy() for e in self.intracellullar_degradation_enzymes]
    
    def get_membrane_degradation_enzymes(self):
        return [e.copy() for e in self.membrane_degradation_enzymes]
    
    def get_signaled_degradation_enzymes(self):
        return [e.copy() for e in self.signaled_degradation_enzymes]
    
    def get_peptidases(self):
        return [e.copy() for e in self.peptidases]
    
    def set_mean_cleavage_cost(self,base_dict):
        if not issubclass(type(base_dict),dict):
            print("ERROR: base_list must be a 'dict' object. Got:",type(base_dict))
            return False
        self.mean_cleavage_cost = base_dict
        return True
    
    def set_mean_cleavage_length(self,base_dict):
        if not issubclass(type(base_dict),dict):
            print("ERROR: base_list must be a 'dict' object. Got:",type(base_dict))
            return False
        self.mean_cleavage_length = base_dict
        return True
    
    def set_energy_molecules(self,mol_list):
        if not issubclass(type(mol_list),list):
            print("ERROR: mol_list must be a 'list' object. Got:",type(mol_list))
            return False
        self.energy_molecules = mol_list
        return True
    
    def get_energy_molecules(self):
        return [c.copy() for c in self.energy_molecules]
        
class DNAReplication(Jsonable):
    def __init__(self):
        self.deoxiribonucleic_acids = []
        self.replication_complex_sufix = ''
        self.active_replication_complex_sufix = ''
        self.replication_start_feature = {}
        self.replication_terminus_feature = {}
        self.replication_complex_composition = []
        self.ppi = None
        self.active_dnaA = None
        self.inactive_dnaA = None
        self.dnaA_complex = []
        self.inactive_dnaA_complex = []
        self.dnaA_replication_boxes = []
        self.energy_molecules = []
        self.dna_repair_enzymes = []
        self.dna_condensation_protein = None
        self.dna_condensation_unbounded_protein = []
        self.dna_condensation_minimal_region_length = None
        
    def set_replication_complex_sufix(self,sufix):
        if not issubclass(type(sufix),str):
            print("ERROR: sufix must be a 'str' object. Got:",type(sufix))
            return False
        self.replication_complex_sufix = sufix
        return True
    
    def set_active_replication_complex_sufix(self,sufix):
        if not issubclass(type(sufix),str):
            print("ERROR: sufix must be a 'str' object. Got:",type(sufix))
            return False
        self.active_replication_complex_sufix = sufix
        return True
                               
    def set_deoxiribonucleic_acids(self,base_dict):
        if not issubclass(type(base_dict),dict):
            print("ERROR: base_list must be a 'dict' object. Got:",type(base_dict))
            return False
        self.deoxiribonucleic_acids = base_dict
        return True
    
    def get_sequence_composition(self,seq,n_phosphates=1):
        dna = []
        for i in range(len(seq)):
            dbase = seq[i]
            if dbase not in ['A','T','C','G']:
                print("WARNING: Base must be 'A','T','C',or 'G'. Got: "+dbase)
                break
            base = self.deoxiribonucleic_acids[dbase][n_phosphates-1]
            dna.append(base)
        return dna
    
    def set_replication_complex_composition(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.replication_complex_composition = prot_list
        return True
    
    def set_dna_repair_enzymes(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.dna_repair_enzymes = prot_list
        return True
    
    def set_ppi(self,ppi):
        if not issubclass(type(ppi),dict):
            print("ERROR: ppi must be a 'dict' object. Got:",type(ppi))
            return False
        self.ppi = ppi
        return True

    def set_sorted_dnaA_complexes(self,dnaa):
        if not issubclass(type(dnaa),list):
            print("ERROR: dnaa must be a 'list' object. Got:",type(dnaa))
            return False
        self.dnaA_complex = dnaa
        return True
    
    def set_sorted_inactive_dnaA_complexes(self,dnaa):
        if not issubclass(type(dnaa),list):
            print("ERROR: dnaa must be a 'list' object. Got:",type(dnaa))
            return False
        self.inactive_dnaA_complex = dnaa
        return True
    
    def set_active_dnaA(self,dnaa):
        if not issubclass(type(dnaa),dict):
            print("ERROR: dnaa must be a 'dict' object. Got:",type(dnaa))
            return False
        self.active_dnaA = dnaa
        return True
    
    def set_inactive_dnaA(self,dnaa):
        if not issubclass(type(dnaa),dict):
            print("ERROR: dnaa must be a 'dict' object. Got:",type(dnaa))
            return False
        self.inactive_dnaA = dnaa
        return True
    
    def set_dna_condensation_protein(self,prot):
        if not issubclass(type(prot),dict):
            print("ERROR: prot must be a 'dict' object. Got:",type(prot))
            return False
        self.dna_condensation_protein = prot
        return True
    
    def set_dna_condensation_unbounded_protein(self,mol_list):
        if not issubclass(type(mol_list),list):
            print("ERROR: mol_list must be a 'list' object. Got:",type(mol_list))
            return False
        self.dna_condensation_unbounded_protein = mol_list
        return True
    
    def set_dna_condensation_minimal_region_length(self,val):
        if not issubclass(type(val),int):
            print("ERROR: val must be a 'int' object. Got:",type(val))
            return False
        self.dna_condensation_minimal_region_length = val
        return True
    
    def set_replication_start_feature(self,start_features):
        if not issubclass(type(start_features),dict):
            print("ERROR: start_features must be a 'dict' object. Got:",type(start_features))
            return False
        self.replication_start_feature = start_features
        return True
    
    def set_replication_terminus_feature(self,terminus_features):
        if not issubclass(type(terminus_features),dict):
            print("ERROR: terminus_features must be a 'dict' object. Got:",type(terminus_features))
            return False
        self.replication_terminus_feature = terminus_features
        return True
    
    def set_dnaA_replication_boxes(self,feat_list):
        if not issubclass(type(feat_list),list):
            print("ERROR: feat_list must be a 'list' object. Got:",type(feat_list))
            return False
        self.dnaA_replication_boxes = feat_list
        return True
    
    def set_energy_molecules(self,mol_list):
        if not issubclass(type(mol_list),list):
            print("ERROR: mol_list must be a 'list' object. Got:",type(mol_list))
            return False
        self.energy_molecules = mol_list
        return True
    
    def get_active_dnaA(self):
        return self.active_dnaA.copy()
    
    def get_inactive_dnaA(self):
        return self.inactive_dnaA.copy()
    
    def get_last_dnaA_complex(self):
        return self.dnaA_complex[-1].copy()
    
    def get_last_inactive_dnaA_complex(self):
        return self.inactive_dnaA_complex[-1].copy()
    
    def get_dna_condensation_protein(self):
        return self.dna_condensation_protein.copy()
    
    def get_dna_condensation_unbounded_protein(self):
        return [c.copy() for c in self.dna_condensation_unbounded_protein]    
    
    def get_ppi(self):
        return self.ppi.copy()
    
    def get_replication_complex_composition(self):
        return [c.copy() for c in self.replication_complex_composition]
    
    def get_energy_molecules(self):
        return [c.copy() for c in self.energy_molecules]
    
    def get_dna_repair_enzymes(self):
        return [c.copy() for c in self.dna_repair_enzymes]
    
        
class RNATranscription(Jsonable):
    def __init__(self):
        self.ribonucleic_acids = []
        self.imature_rna_prefix = ''
        self.tu_transcript_sufix = ''
        self.oligonucleotide_prefix = ''
        self.modification_reactions = {}
        self.rna_cleavage_enzymes = {}
        self.cleavage_energy_molecules = []
        self.cost_per_cleavage = {}
        self.dna_sufix = ''
        self.transcription_complex_sufix = ''
        self.active_transcription_complex_sufix = ''
        self.transcription_release_factors = []
        self.transcription_elongation_factors = []
        self.rna_polymerase = None
        self.rna_polymerase_holoenzyme = None
        self.sigma_factor = None
        self.ppi = None
        
    def set_imature_rna_prefix(self,prefix):
        if not issubclass(type(prefix),str):
            print("ERROR: prefix must be a 'str' object. Got:",type(prefix))
            return False
        self.imature_rna_prefix = prefix
        return True
    
    def set_oligonucleotide_prefix(self,prefix):
        if not issubclass(type(prefix),str):
            print("ERROR: prefix must be a 'str' object. Got:",type(prefix))
            return False
        self.oligonucleotide_prefix = prefix
        return True
    
    def set_tu_transcript_sufix(self,sufix):
        if not issubclass(type(sufix),str):
            print("ERROR: sufix must be a 'str' object. Got:",type(sufix))
            return False
        self.tu_transcript_sufix = sufix
        return True
    
    def set_dna_sufix(self,sufix):
        if not issubclass(type(sufix),str):
            print("ERROR: sufix must be a 'str' object. Got:",type(sufix))
            return False
        self.dna_sufix = sufix
        return True
    
    def set_transcription_complex_sufix(self,sufix):
        if not issubclass(type(sufix),str):
            print("ERROR: sufix must be a 'str' object. Got:",type(sufix))
            return False
        self.transcription_complex_sufix = sufix
        return True
    
    def set_active_transcription_complex_sufix(self,sufix):
        if not issubclass(type(sufix),str):
            print("ERROR: sufix must be a 'str' object. Got:",type(sufix))
            return False
        self.active_transcription_complex_sufix = sufix
        return True
                               
    def set_ribonucleic_acids(self,base_dict):
        if not issubclass(type(base_dict),dict):
            print("ERROR: base_list must be a 'dict' object. Got:",type(base_dict))
            return False
        self.ribonucleic_acids = base_dict
        return True
        
    def add_modification_reaction(self,react_obj):
        if not issubclass(type(react_obj),RNAModificationReaction):
            print("ERROR: object must inherit 'RNAModificationReaction' class. Got:",type(react_obj))
            return False
        if react_obj.ref_name in self.modification_reactions.keys():
            print("ERROR: 'ref_name' already exists:",react_obj.ref_name)
            return False
        self.modification_reactions[react_obj.ref_name] = react_obj
        return True
    
    def add_modification_reactions(self,react_obj_list): 
        if type(react_obj_list) != list:
            print("ERROR: 'react_obj_list' must be a list of 'RNAModificationReaction' class inherited objects. Got:",type(react_obj_list))
            return False         
        return sum([self.add_modification_reaction(react_obj) for react_obj in react_obj_list])
        
    def is_base(self,mol):
        if mol in [v[0] for v in self.ribonucleic_acids.values()]: return True
    
    def is_modified_base(self,mol):
        if not 'MP' in mol: return False
        flag = False
        for b in self.ribonucleic_acids.keys():
            if b in mol:
                flag=True
                break
        if not flag: return False
        if not len(mol)>3:
            return False
        return True
        
    def apply_modification_reaction(self,base,mod):
        if mod not in self.modification_reactions.keys():
            print("WARNING: Modification reaction not found. Modification Reaction:",mod)
            return False
        del_base = self.modification_reactions[mod].del_base
        add_base = self.modification_reactions[mod].add_base
        position = self.modification_reactions[mod].position
        if position != None:
            position = int(position)-1
            if base[position]==del_base:
                if add_base!=None:
                    base[position] = add_base
                    return True
            else:
                print("WARNING: base "+del_base+" not found in position. Modification:",mod,"Base:",base[position],"Position:",position+1)
                return False
        return False
    
    def get_modification_enzymes(self,mod):
        return [e.copy() for e in self.modification_reactions[mod].enzymes]
    
    def get_auxiliar_molecules(self,mod):
        return [e.copy() for e in self.modification_reactions[mod].auxiliar_molecules]
        
    def get_modification_molecules(self,mod):
        return [e.copy() for e in self.modification_reactions[mod].molecules]
    
    def get_modification_annotations(self,mod):
        return self.modification_reactions[mod].annotations
    
    def set_rna_cleavage_enzymes(self,prot_dict):
        if not issubclass(type(prot_dict),dict):
            print("ERROR: prot_dict must be a 'dict' object. Got:",type(prot_dict))
            return False
        self.rna_cleavage_enzymes = prot_dict
        return True
    
    def get_rna_cleavage_enzymes(self):
        return {key:self.rna_cleavage_enzymes[key].copy() for key in self.rna_cleavage_enzymes.keys()}
    
    def set_cleavage_energy_molecules(self,mol_list):
        if not issubclass(type(mol_list),list):
            print("ERROR: mol_list must be a 'list' object. Got:",type(mol_list))
            return False
        self.cleavage_energy_molecules = mol_list
        return True
    
    
    def set_cost_per_cleavage(self,cost):
        if not issubclass(type(cost),dict):
            print("ERROR: cost must be a 'int' object. Got:",type(cost))
            return False
        self.cost_per_cleavage = cost
        return True
    
    def get_cleavage_energy_molecules(self):
        return [m.copy() for m in self.cleavage_energy_molecules]
                               
    def set_transcription_release_factors(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.transcription_release_factors = prot_list
        return True
    
    def get_transcription_release_factors(self):
        return [e.copy() for e in self.transcription_release_factors]
    
    def set_transcription_elongation_factors(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.transcription_elongation_factors = prot_list
        return True
    
    def get_transcription_elongation_factors(self):
        return [e.copy() for e in self.transcription_elongation_factors]
    
    def set_rna_polymerase(self,rnapoly):
        if not issubclass(type(rnapoly),dict):
            print("ERROR: rnapoly must be a 'dict' object. Got:",type(rnapoly))
            return False
        self.rna_polymerase = rnapoly
        return True
    
    def get_rna_polymerase(self):
        return self.rna_polymerase.copy()
    
    def set_rna_polymerase_holoenzyme(self,rnapoly):
        if not issubclass(type(rnapoly),dict):
            print("ERROR: rnapoly must be a 'dict' object. Got:",type(rnapoly))
            return False
        self.rna_polymerase_holoenzyme = rnapoly
        return True
    
    def get_rna_polymerase_holoenzyme(self):
        return self.rna_polymerase_holoenzyme.copy()
    
    def set_sigma_factor(self,sigma):
        if not issubclass(type(sigma),dict):
            print("ERROR: sigma must be a 'dict' object. Got:",type(sigma))
            return False
        self.sigma_factor = sigma
        return True
    
    def get_sigma_factor(self):
        return self.sigma_factor.copy()
    
    def set_ppi(self,ppi):
        if not issubclass(type(ppi),dict):
            print("ERROR: ppi must be a 'dict' object. Got:",type(ppi))
            return False
        self.ppi = ppi
        return True
    
    def get_ppi(self):
        return self.ppi.copy()
                               
    def transcribe_gene(self,seq):
        return complimentary_sequence(seq,target='RNA')
        
    
    def get_rna_composition(self,seq,n_phosphates=1):
        rna = []
        for i in range(len(seq)):
            base = seq[i]
            if base not in ['A','U','C','G']:
                print("WARNING: Base must be 'A','U','C',or 'G'. Got: "+base)
                break
            rna.append(self.ribonucleic_acids[base][n_phosphates-1])
        return rna

            
    
class ProteinTranslation(Jsonable):
    def __init__(self,_translation_table):
        self.ribosome = None
        self.translation_table = _translation_table
        self.aminoacids_names = list(set([codon['aminoacid'] for codon in self.translation_table.values()]))
        self.translation_complex_constitution = []
        self.translation_complex_prefix = ''
        self.stalled_translation_complex_sufix = ''
        self.initiation_auxiliaries = []
        self.elongation_auxiliaries = []
        self.elongation_energy_molecules = []
        self.imature_protein_prefix = ''
        self.mature_protein_prefix = ''
        self.modification_reactions = {}
        self.translocation_proteins = []
        self.translocation_enzymes = []
        self.get_cost_for_translocation = None
        self.translation_stall_auxiliaries = []
        self.translation_stall_tmRNA = []
        self.peptide_prefix = ''
        self.tmRNA_proteolysis_tag_chromosome_feature = ''
        self.aminoacid_letter = {
            'GLY' : 'G',
            'ALA' : 'A',
            'LEU' : 'L',
            'MET' : 'M',
            'PHE' : 'F',
            'TRP' : 'W',
            'LYS' : 'K',
            'GLN' : 'Q',
            'GLU' : 'E',
            'SER' : 'S',
            'PRO' : 'P',
            'VAL' : 'V',
            'ILE' : 'I',
            'CYS' : 'C',
            'TYR' : 'Y',
            'HIS' : 'H',
            'ARG' : 'R',
            'ASN' : 'N',
            'ASP' : 'D',
            'THR' : 'T'
        }
        
    def set_ribosome(self,rib):
        if not issubclass(type(rib),dict):
            print("ERROR: rib must be a 'dict' object. Got:",type(rib))
            return False
        self.ribosome = rib
        return True
        
    def set_translation_complex_constitution(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.translation_complex_constitution = prot_list
        return True
    
    def set_tmRNA_proteolysis_tag_chromosome_feature(self,feat):
        if not issubclass(type(feat),str):
            print("ERROR: feat must be a 'str' object. Got:",type(feat))
            return False
        self.tmRNA_proteolysis_tag_chromosome_feature = feat
        return True
    
    def set_translation_complex_prefix(self,prefix):
        if not issubclass(type(prefix),str):
            print("ERROR: prefix must be a 'str' object. Got:",type(prefix))
            return False
        self.translation_complex_prefix = prefix
        return True
                                   
    def set_stalled_translation_complex_sufix(self,sufix):
        if not issubclass(type(sufix),str):
            print("ERROR: sufix must be a 'str' object. Got:",type(sufix))
            return False
        self.stalled_translation_complex_sufix = sufix
        return True                               
    
    def set_initiation_auxiliaries(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.initiation_auxiliaries = prot_list
        return True
    
    def set_elongation_auxiliaries(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.elongation_auxiliaries = prot_list
        return True
    
    def set_elongation_energy_molecules(self,mol_list):
        if not issubclass(type(mol_list),list):
            print("ERROR: mol_list must be a 'list' object. Got:",type(mol_list))
            return False
        self.elongation_energy_molecules = mol_list
        return True
     
    def set_imature_protein_prefix(self,prefix):
        if not issubclass(type(prefix),str):
            print("ERROR: prefix must be a 'str' object. Got:",type(prefix))
            return False
        self.imature_protein_prefix = prefix
        return True
    
    def set_peptide_prefix(self,prefix):
        if not issubclass(type(prefix),str):
            print("ERROR: prefix must be a 'str' object. Got:",type(prefix))
            return False
        self.peptide_prefix = prefix
        return True
    
    def set_translocation_proteins(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.translocation_proteins = prot_list
        return True
    
    def set_translocation_enzymes(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.translocation_enzymes = prot_list
        return True
    
    def set_gtp_cost_for_translocation(self,cost):
        if not issubclass(type(cost),int):
            print("ERROR: cost must be a 'int' object. Got:",type(cost))
            return False
        self.gtp_cost_for_translocation = cost
        return True
    
    def set_translation_stall_auxiliaries(self,prot_list):
        if not issubclass(type(prot_list),list):
            print("ERROR: prot_list must be a 'list' object. Got:",type(prot_list))
            return False
        self.translation_stall_auxiliaries = prot_list
        return True
    
    def set_translation_stall_tmRNA(self,tmRNA_list):
        if not issubclass(type(tmRNA_list),list):
            print("ERROR: tmRNA_list must be a 'list' object. Got:",type(tmRNA_list))
            return False
        self.translation_stall_tmRNA = tmRNA_list
        return True
    
    def get_translation_stall_tmRNA(self):
        return [m.copy() for m in self.translation_stall_tmRNA]
    
    def get_translation_stall_auxiliaries(self):
        return [m.copy() for m in self.translation_stall_auxiliaries]
    
    def get_translocation_proteins(self):
        return [m.copy() for m in self.translocation_proteins]
    
    def get_translocation_enzymes(self):
        return [m.copy() for m in self.translocation_enzymes]
    
    def get_gtp_cost_for_translocation(self):
        return self.gtp_cost_for_translocation
    
    def get_elongation_auxiliaries(self):
        return [m.copy() for m in self.elongation_auxiliaries]
    
    def get_initiation_auxiliaries(self):
        return [m.copy() for m in self.initiation_auxiliaries]
    
    def get_translation_complex_constitution(self):
        return [m.copy() for m in self.translation_complex_constitution]
    
    def get_elongation_energy_molecules(self):
        return [m.copy() for m in self.elongation_energy_molecules]
    
    def get_ribosome(self):
        return self.ribosome.copy()
        
    def add_modification_reaction(self,react_obj):
        if not issubclass(type(react_obj),ProteinModificationReaction):
            print("ERROR: object must inherit 'ProteinModificationReaction' class. Got:",type(react_obj))
            return False
        if react_obj.ref_name in self.modification_reactions.keys():
            print("ERROR: 'ref_name' already exists:",react_obj.ref_name)
            return False
        self.modification_reactions[react_obj.ref_name] = react_obj
        return True
    
    def add_modification_reactions(self,react_obj_list): 
        if type(react_obj_list) != list:
            print("ERROR: 'react_obj_list' must be a list of 'ProteinModificationReaction' class inherited objects. Got:",type(react_obj_list))
            return False         
        return sum([self.add_modification_reaction(react_obj) for react_obj in react_obj_list])
        
    def translate_gene(self,seq,trans_out='aa',mark_initial_codon=True): #trans_out = 'aa' or 'aatrna' or 'trna'
        aa = []
        if len(seq)%3 != 0:
            print("ERROR: Sequence length is not multiple of 3. Length:",len(seq))
            return aa
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if i==0 and mark_initial_codon:
                codon = 'i'+codon
            if codon not in self.translation_table.keys():
                if i<len(seq)-3:
                    print("WARNING: Stop codon reached before end of sequence. Codon:",codon,"Total Length:",len(seq),"Stopped at:",i)
                    break
                else:
                    continue
            if trans_out == 'aa':
                aminoacid = self.translation_table[codon]['aminoacid']
            elif trans_out == 'aatrna':
                aminoacid = self.translation_table[codon]['tRNA']+'_'+self.translation_table[codon]['aminoacid']
            elif trans_out == 'trna':
                aminoacid = self.translation_table[codon]['tRNA']
            aa.append(aminoacid)
        return aa

    def is_aminoacid(self,a):
        return a in self.aminoacids_names
    
    def is_modified_aminoacid(self,a):
        for aa in self.aminoacids_names:
            if a[-len(aa)::].upper() == aa and len(a)>len(aa): return True           
        return False
    
    def apply_modification_reaction(self,aa,mod):
        if mod not in self.modification_reactions.keys():
            print("WARNING: Modification reaction not found. Modification Reaction:",mod)
            return False
        del_aa = self.modification_reactions[mod].del_aa
        add_aa = self.modification_reactions[mod].add_aa
        position = self.modification_reactions[mod].position
        if position != None:
            position = int(position)-1
            if aa[position]==del_aa:
                if add_aa!=None:
                    aa[position] = add_aa
                else:
                    aa.pop(position)
            else:
                print("WARNING:aminoacid not found in position. Modification:",mod,"Aminoacid:",del_aa,"Position:",position)
        else:
            aa_found = False
            for i in range(len(aa)):
                if aa[i]==del_aa:
                    if add_aa!=None:
                        aa[i] = add_aa
                    else:
                        aa.pop(i)
                    aa_found = True
                    break
            if not aa_found: print("WARNING:aminoacid not found in sequence. Aminoacid:",del_aa)
        return True
    
    def get_modification_enzymes(self,mod):
        return [e.copy() for e in self.modification_reactions[mod].enzymes]
    
    def get_auxiliar_molecules(self,mod):
        return [m.copy() for m in self.modification_reactions[mod].auxiliar_molecules]
        
    def get_modification_molecules(self,mod):
        return [m.copy() for m in self.modification_reactions[mod].molecules]
    
    def get_modification_annotations(self,mod):
        return self.modification_reactions[mod].annotations

        
        
        
class WholeCellKBHandlerJSON:
    def __init__(self,file_path):
        try:
            self.f = open(file_path)
        except:
            print("ERROR: File could not be opened:",path)
        try:
            self.wcdata = json.load(self.f)
        except:
            print("ERROR: Failed to read json file.")
            
        self.ignore_list = []
        self.name_mapping = {}
        self.entry_index = {self.wcdata['data'][i]['wid']:i for i in range(len(self.wcdata['data']))}
        self.database_name = "WholeCellKB"
        self.created_molecules = {}
        
    def set_ignore_list(self,path):
        with open(path) as f:
            for line in f.readlines():
                if line[0]!='#':
                    self.ignore_list.append(line.strip())
                
    def set_name_mapping(self,path):
        with open(path) as f:
            for line in f.readlines():
                if line[0]!='#' and line!='':
                    vals = line.split('\t')
                    self.name_mapping[vals[0]]=vals[1].strip()
    
    def get_database_name(self):
        return self.database_name
    
    def get_models(self):
        return list(set([entry['model'] for entry in self.wcdata['data']]))
    
    def get_entry(self,idx):
        return self.wcdata['data'][idx]
    
    def have_entry(self,entry_name):
        return entry_name in self.entry_index.keys()
    
    def have_created_molecule(self,mol):
        return mol in self.created_molecules.keys()
    
    def find_entry(self,e):
        if e not in self.entry_index.keys():
            return None
        return self.get_entry(self.entry_index[e])
    
    def get_chemical_reactions(self,reacts_names):
        reacts = []
        for r in reacts_names:
            react_entry = self.find_entry(r)
            if react_entry['model']!='Reaction':
                print("ERROR: entry 'model' must be 'Reaction'. Got:",react_entry['model'])
                return False
            react = ChemicalReaction()
            react.ref_name = react_entry['wid']
            react.name = react_entry['name']
            react.cross_references = {cr["source"]:cr["xid"] for cr in react_entry['cross_references']}
            react.annotations = ['Chemical Reaction']
            for ann in react_entry['type']:
                react.annotations.append(ann)
            react.original_database = self.get_database_name()
            react.process = react_entry['processes']
            if react_entry['processes']=="Process_ReplicationInitiation":
                react.process = "Replication Initiation"
            if react_entry['enzyme']!=None:
                mol = react_element(react_entry['enzyme']['protein'],react_entry['enzyme']['compartment'])
                react.enzymes.append(mol)
            for m in react_entry['stoichiometry']:
                mol = react_element(m['molecule'],m['compartment'],int(float(m['coefficient'])))
                if m['molecule']in self.name_mapping.keys():
                    mol['molecule'] = self.name_mapping[m['molecule']]
                react.molecules.append(mol)
                if mol['molecule'] == "MG_287_MONOMER_ACP" and mol['stoichiometry']>0:                
                    react.molecules.append(react_element("H2O",'c'))
                    react.molecules.append(react_element("H",'c',-1))
            for m in react_entry['coenzymes']:
                mol = react_element(m['metabolite'],m['compartment'])
                if m['metabolite']in self.name_mapping.keys():
                    mol['molecule'] = self.name_mapping[m['metabolite']]
                react.auxiliar_molecules.append(mol)
            react.is_reversible = react_entry['direction'] == 'Reversible'
            if react_entry['kinetics_forward'] != None:
                rate,param = self.get_kinetic_law(react_entry['kinetics_forward'])
                react.rate_law_forward = rate
                react.parameters_forward = param
            if react_entry['kinetics_backward'] != None:
                rate,param = self.get_kinetic_law(react_entry['kinetics_backward'])
                react.rate_law_backward = rate
                react.parameters_backward = param
            if react_entry['optimal_ph']!=None:
                react.optimal_ph = react_entry['optimal_ph']['value']
            if react_entry['optimal_temperature']!=None:
                react.optimal_temperature = " ".join([react_entry['optimal_temperature']['value'],react_entry['optimal_temperature']['units']])
            react.pathways = react_entry['pathways']
            reacts.append(react)
        return reacts
    
    def get_trna_list(self):
        glist = [entry for entry in self.wcdata['data'] if entry['model']=='Gene' and entry['wid'] not in self.ignore_list]
        return list(set([entry['wid'] for entry in glist if entry['amino_acid']!=None]))
    
    def get_translation_table(self):
        table = {}
        trna_list = self.get_trna_list()
        for trna in trna_list:
            trna_entry = self.get_entry(self.entry_index[trna])
            for codon in trna_entry['codons']:
                if trna_entry['amino_acid']=='FMET':
                    codon['sequence'] = 'i'+codon['sequence']
                if codon['sequence'] not in table:
                    table[codon['sequence']] = {'tRNA':trna_entry['wid'],'aminoacid':trna_entry['amino_acid']}
                else:
                    print("ERROR: Multiple entry of codon:",codon['sequence'],table[codon['sequence']],trna_entry['wid']+'_'+trna_entry['amino_acid'])
        return table
    
    def get_chromosome_features(self):
        feats = {}
        for feat in [entry for entry in self.wcdata['data'] if entry['model']=='ChromosomeFeature' and entry['wid'] not in self.ignore_list]:
            f = {}
            f['ref_name'] = feat['wid']
            f['name'] = feat['name']
            f['start'] = int(feat['coordinate'])
            f['length'] = int(feat['length'])
            f['direction'] = feat['direction']
            f['chromosome'] = feat['chromosome']
            f['type'] = None
            if len(feat['type']):
                f['type'] = feat['type'][0]
            feats[f['ref_name']] = f
        return feats
    
    def get_chromosome_transcriptional_regulation(self):
        regs = []
        for reg in [entry for entry in self.wcdata['data'] if entry['model']=='ChromosomeFeature' and entry['wid'] not in self.ignore_list]:
            f = {}
            f['ref_name'] = entry['wid']
            f['name'] = entry['name']
            f['start'] = entry['coordinate']
            f['length'] = entry['length']
            f['direction'] = entry['direction']
            f['type'] = None
            if len(entry['type']):
                f['type'] = entry['type'][0]
            feats.apend(f)
        return feats
    
    
    def get_chromosome_list(self):
        return list(set([entry['wid'] for entry in self.wcdata['data'] if entry['model']=='Chromosome' and entry['wid'] not in self.ignore_list]))
    
    def make_chromosome(self,chrm_entry):
        if chrm_entry['model']!='Chromosome':
            print("ERROR: entry 'model' must be 'Chromosome'. Got:",chrm_entry['model'])
            return False
        chrm = Chromosome()
        chrm.ref_name = chrm_entry['wid']
        chrm.name = chrm_entry['name']
        chrm.sequence = chrm_entry['sequence']
        chrm.length = int(chrm_entry['length'])
        chrm.genes = chrm_entry['genes']
        chrm.original_database = self.get_database_name()
        return chrm
    
    def create_all_chromosomes(self):
        chrms = self.get_chromosome_list()
        chrm_objs = []
        for i in [self.entry_index[c] for c in chrms]:
            chrm_objs.append(self.make_chromosome(self.get_entry(i)))
        return chrm_objs
    
    def get_transcription_unit_list(self):
        return list(set([entry['wid'] for entry in self.wcdata['data'] if entry['model']=='TranscriptionUnit' and entry['wid'] not in self.ignore_list]))
    
    def get_transcriptional_regulations_for_tu(self,tu):
        return [entry for entry in self.wcdata['data'] if entry['model']=='TranscriptionalRegulation' and entry['transcription_unit']==tu and entry['wid'] not in self.ignore_list]
    
    def make_transcription_unit(self,tu_entry):
        if tu_entry['model']!='TranscriptionUnit':
            print("ERROR: entry 'model' must be 'TranscriptionUnit'. Got:",gen_entry['model'])
            return False
        tu = TranscriptionUnit()
        tu.ref_name = tu_entry['wid']
        tu.name = tu_entry['name']
        tu.genes = tu_entry['genes']
        if tu_entry['promoter_10_coordinate']!= 'None':
            tu.promoter_10_start = int(tu_entry['promoter_10_coordinate'])
        if tu_entry['promoter_10_length']!= 'None':
            tu.promoter_10_length = int(tu_entry['promoter_10_length'])
        if tu_entry['promoter_35_coordinate']!= 'None':
            tu.promoter_35_start = int(tu_entry['promoter_35_coordinate'])
        if tu_entry['promoter_35_length']!= 'None':
            tu.promoter_35_length = int(tu_entry['promoter_35_length'])
        if tu_entry['tss_coordinate']!= 'None':
            tu.transcription_start = int(tu_entry['tss_coordinate'])
        tu.original_database = self.get_database_name()
        tu.is_polycistronic = len(tu.genes)>1
        for reg in self.get_transcriptional_regulations_for_tu(tu.ref_name):
            tu.transcription_regulators[reg['transcription_factor']] = float(reg['activity']['value'])
        return tu
    
    def create_all_transcription_units(self):
        tus = self.get_transcription_unit_list()
        tu_objs = []
        for i in [self.entry_index[tu] for tu in tus]:
            tu_objs.append(self.make_transcription_unit(self.get_entry(i)))
        return tu_objs   
    
    
    def get_gene_list(self):
        return list(set([entry['wid'] for entry in self.wcdata['data'] if entry['model']=='Gene' and entry['wid'] not in self.ignore_list]))
       
        
    def make_gene(self,gen_entry):
        if gen_entry['model']!='Gene':
            print("ERROR: entry 'model' must be 'Gene'. Got:",gen_entry['model'])
            return False
        gen = Gene()
        gen.ref_name = gen_entry['wid']
        gen.name = gen_entry['name']
        gen.symbol = gen_entry['symbol']
        gen.chromosome = gen_entry['chromosome']
        gen.transcription_unit = gen_entry['transcription_units'][0]
        gen.start = int(gen_entry['coordinate'])
        gen.length = int(gen_entry['length'])
        gen.direction = gen_entry['direction']
        gen.type = gen_entry['type'][0]
        gen.original_database = self.get_database_name()
        #for mod in gen_entry['modification_reactions']:
        #    if len(mod['reactions'])>0:
        #        for r in mod['reactions']:
        #            if 'MODIFICATION' in r:
        #                gen.modification_reactions.append(r)
        return gen
        
    def create_all_genes(self):
        genes = self.get_gene_list()
        gen_objs = []
        for i in [self.entry_index[g] for g in genes]:
            gen_objs.append(self.make_gene(self.get_entry(i)))
        return gen_objs
    
    def get_terminal_organelle_assembly_list(self):
        return list(set([entry['wid'] for entry in self.wcdata['data'] if entry['model']=='Reaction' and entry['processes']=='Process_TerminalOrganelleAssembly']))
    
    
    def make_terminal_organelle_assembly_reaction(self,react_entry):
        if react_entry['model']!='Reaction':
            print("ERROR: entry 'model' must be 'Reaction'. Got:",react_entry['model'])
            return False
        if react_entry['processes']!='Process_TerminalOrganelleAssembly':
            print("ERROR: entry 'processes' must be 'Process_TerminalOrganelleAssembly'. Got:",react_entry['processes'])
            return False
        react = TransportReaction()
        react.ref_name = react_entry['wid']
        react.name = react_entry['name']
        react.cross_references = {cr["source"]:cr["xid"] for cr in react_entry['cross_references']}
        react.annotations = ['Transport Reaction']
        for ann in react_entry['type']:
            react.annotations.append(ann)
        react.original_database = self.get_database_name()
        react.process = 'Terminal Organelle Assembly'
        if react_entry['enzyme']!=None:
            mol = react_element(react_entry['enzyme']['protein'],react_entry['enzyme']['compartment'])
            react.enzymes.append(mol)
        react.is_reversible = len(react.enzymes)==0
        prot_from = react_element(react_entry['modification']['molecule'],react_entry['modification']['compartment'],-1)
        prot_to = react_element(react_entry['modification']['molecule'],'t'+react_entry['modification']['compartment'])
        react.molecules = [prot_from,prot_to]
        react.from_compartment = prot_from['location']
        react.to_compartment = prot_to['location']
        return react
    
    def create_all_terminal_organelle_assembly_reactions(self):
        reactions = self.get_terminal_organelle_assembly_list()
        react_objs = []
        for i in [self.entry_index[r] for r in reactions]:
            react_objs.append(self.make_terminal_organelle_assembly_reaction(self.get_entry(i)))
        return react_objs    
        
    def make_metabolite(self,met_entry):
        if met_entry['model']!='Metabolite':
            print("ERROR: entry 'model' must be 'Metabolite'. Got:",met_entry['model'])
            return False
        met = SmallMolecule()
        met.ref_name = met_entry['wid']
        met.name = met_entry['name']
        met.volume = met_entry['volume']
        met.cross_references = {cr["source"]:cr["xid"] for cr in met_entry['cross_references']}
        met.annotations = ['Metabolite']
        for ann in met_entry['type']:
            met.annotations.append(ann)
        met.original_database = self.get_database_name()
        met.formula = met_entry['empirical_formula']
        met.compute_mw_from_formula()
        met.smiles = met_entry['smiles']
        met.iupac_name = met_entry['iupac_name']
        met.is_hydrophobic = met_entry['is_hydrophobic']
        met.charge = met_entry['charge']
        met.synonyms = met_entry['synonyms']
        met.log_d = met_entry['log_d']
        met.log_p = met_entry['log_p']
        met.deltag_formation = met_entry['deltag_formation']
        return met
        
    
    def get_metabolic_reactions_list(self):
        return list(set([entry['wid'] for entry in self.wcdata['data'] if entry['model']=='Reaction' and entry['processes']=='Process_Metabolism' and entry['wid'] not in self.ignore_list]))
                
    def make_metabolic_reaction(self,react_entry):
        if react_entry['model']!='Reaction':
            print("ERROR: entry 'model' must be 'Reaction'. Got:",react_entry['model'])
            return False
        if react_entry['processes']!='Process_Metabolism':
            print("ERROR: entry 'processes' must be 'Process_Metabolism'. Got:",react_entry['processes'])
            return False
        react = ChemicalReaction()
        react.ref_name = react_entry['wid']
        react.name = react_entry['name']
        react.cross_references = {cr["source"]:cr["xid"] for cr in react_entry['cross_references']}
        react.annotations = ['Metabolic Reaction']
        for ann in react_entry['type']:
            react.annotations.append(ann)
        react.original_database = self.get_database_name()
        react.process = 'Metabolism'
        if react_entry['enzyme']!=None:
            mol = react_element(react_entry['enzyme']['protein'],react_entry['enzyme']['compartment'])
            react.enzymes.append(mol)
        for m in react_entry['stoichiometry']:
            mol = react_element(m['molecule'],m['compartment'],int(float(m['coefficient'])))
            if m['molecule']in self.name_mapping.keys():
                mol['molecule'] = self.name_mapping[m['molecule']]
            react.molecules.append(mol)
            if mol['molecule'] == "MG_287_MONOMER_ACP":                
                react.molecules.append(react_element("H2O",'c',mol['stoichiometry']))
                react.molecules.append(react_element("H",'c',-mol['stoichiometry']))
        for m in react_entry['coenzymes']:
            mol = react_element(m['metabolite'],m['compartment'])
            if m['metabolite']in self.name_mapping.keys():
                mol['metabolite'] = self.name_mapping[m['metabolite']]
            react.auxiliar_molecules.append(mol)
        react.is_reversible = react_entry['direction'] == 'Reversible'
        if react_entry['kinetics_forward'] != None:
            rate,param = self.get_kinetic_law(react_entry['kinetics_forward'])
            react.rate_law_forward = rate
            react.parameters_forward = param
        if react_entry['kinetics_backward'] != None:
            rate,param = self.get_kinetic_law(react_entry['kinetics_backward'])
            react.rate_law_backward = rate
            react.parameters_backward = param
        if react_entry['optimal_ph']!=None:
            react.optimal_ph = react_entry['optimal_ph']['value']
        if react_entry['optimal_temperature']!=None:
            react.optimal_temperature = " ".join([react_entry['optimal_temperature']['value'],react_entry['optimal_temperature']['units']])
        react.pathways = react_entry['pathways']
        return react
        
    def get_kinetic_law(self, k):
        rate_law = ''
        parameters = {}
        if k['rate_law'] != '':
            rate_law = k['rate_law']
            parameters['vmax'] = ' '.join([k['vmax'],k['vmax_unit']])
            km = k['km'].split(',')
            if len(km)>1:
                for i in range(len(km)):
                    parameters['Km'+str(i)] = ' '.join([km[i].strip(),'uM'])
            else:
                parameters['Km'] = ' '.join([km[0].strip(),'uM'])
        elif len(k['evidence'])>0:
            text = k['evidence'][0]['value'].split(';')
            if len(text)==2:
                rate_law = text[0][11:]
                key = ''
                counter = 1
                for val in text[1].replace('=',',').split(','):
                    val = val.strip()
                    if val == 'Km' or val == 'Vmax':
                        key = val
                        continue
                    if key == 'Vmax':
                        parameters[key] = val
                        continue
                    if key == 'Km':
                        parameters[key+str(counter)] = val
                        counter += 1
                        continue                
                    
        return rate_law,parameters
        
    def create_all_metabolic_reactions(self):
        reactions = self.get_metabolic_reactions_list()
        react_objs = []
        for i in [self.entry_index[r] for r in reactions]:
            react_objs.append(self.make_metabolic_reaction(self.get_entry(i)))
        return react_objs
    
    def get_rna_modification_reactions_list(self):
        return list(set([entry['wid'] for entry in self.wcdata['data'] if entry['model']=='Reaction' and entry['processes'] in ['Process_RNAModification'] and entry['wid'] not in self.ignore_list]))
    
    
    def make_rna_modification_reaction(self,react_entry,transcription):
        if react_entry['model']!='Reaction':
            print("ERROR: entry 'model' must be 'Reaction'. Got:",react_entry['model'])
            return False
        if react_entry['processes']not in ['Process_RNAModification']:
            print("ERROR: entry 'processes' must be 'Process_RNAModification'. Got:",react_entry['processes'])
            return False
        react = RNAModificationReaction()
        react.ref_name = react_entry['wid']
        react.name = react_entry['name']
        react.cross_references = {cr["source"]:cr["xid"] for cr in react_entry['cross_references']}
        react.annotations = ['RNA Modification Reaction']
        for ann in react_entry['type']:
            react.annotations.append(ann)
        react.original_database = self.get_database_name()
        react.process = 'RNA Modification'
        if react_entry['modification']!=None:
            react.position = react_entry['modification']['position']
        if react_entry['enzyme']!=None:
            mol = react_element(react_entry['enzyme']['protein'],react_entry['enzyme']['compartment'])
            react.enzymes.append(mol)
        for m in react_entry['stoichiometry']:
            mol = react_element(m['molecule'],m['compartment'],int(float(m['coefficient'])))
            if m['molecule']in self.name_mapping.keys():
                mol['molecule'] = self.name_mapping[m['molecule']]
            react.molecules.append(mol)

            if float(mol['stoichiometry'])<0 and (transcription.is_base(mol['molecule']) or transcription.is_modified_base(mol['molecule'])):
                react.del_base = mol['molecule']
                continue
            if mol['stoichiometry']>0 and transcription.is_modified_base(mol['molecule']):
                react.add_base = mol['molecule']
                continue
            react.molecules.append(mol)
        for m in react_entry['coenzymes']:
            mol = react_element(m['metabolite'],m['compartment'])
            if m['metabolite']in self.name_mapping.keys():
                mol['metabolite'] = self.name_mapping[m['metabolite']]
            react.auxiliar_molecules.append(mol)
        react.is_reversible = react_entry['direction'] == 'Reversible'
        if react_entry['kinetics_forward'] != None:
            rate,param = self.get_kinetic_law(react_entry['kinetics_forward'])
            react.rate_law_forward = rate
            react.parameters_forward = param
        if react_entry['kinetics_backward'] != None:
            rate,param = self.get_kinetic_law(react_entry['kinetics_backward'])
            react.rate_law_backward = rate
            react.parameters_backward = param
        if react_entry['optimal_ph']!=None:
            react.optimal_ph = react_entry['optimal_ph']['value']
        if react_entry['optimal_temperature']!=None:
            react.optimal_temperature = " ".join([react_entry['optimal_temperature']['value'],react_entry['optimal_temperature']['units']])
        react.pathways = react_entry['pathways']
        return react
    
    
    def create_all_rna_modification_reactions(self,transcription):
        reactions = self.get_rna_modification_reactions_list()
        react_objs = []
        for i in [self.entry_index[r] for r in reactions]:
            react_objs.append(self.make_rna_modification_reaction(self.get_entry(i),transcription))
        return react_objs
    
    
    def get_protein_modification_reactions_list(self):
        return list(set([entry['wid'] for entry in self.wcdata['data'] if entry['model']=='Reaction' and entry['processes'] in ['Process_ProteinModification','Process_ProteinProcessingI','Process_ProteinProcessingII','Process_ProteinTranslocation'] and entry['wid'] not in self.ignore_list]))
    
    def make_protein_modification_reaction(self,react_entry,translation):
        if react_entry['model']!='Reaction':
            print("ERROR: entry 'model' must be 'Reaction'. Got:",react_entry['model'])
            return False
        if react_entry['processes']not in ['Process_ProteinModification','Process_ProteinProcessingI','Process_ProteinProcessingII','Process_ProteinTranslocation']:
            print("ERROR: entry 'processes' must be 'Process_ProteinModification' or 'Process_ProteinProcessingI' or 'Process_ProteinProcessingII' or 'Process_ProteinTranslocation'. Got:",react_entry['processes'])
            return False
        react = ProteinModificationReaction()
        react.ref_name = react_entry['wid']
        react.name = react_entry['name']
        react.cross_references = {cr["source"]:cr["xid"] for cr in react_entry['cross_references']}
        react.annotations = ['Protein Modification Reaction']
        for ann in react_entry['type']:
            react.annotations.append(ann)
        react.original_database = self.get_database_name()
        react.process = 'Protein Modification'
        if react_entry['modification']!=None:
            react.position = react_entry['modification']['position']
        if react_entry['enzyme']!=None:
            mol = react_element(react_entry['enzyme']['protein'],react_entry['enzyme']['compartment'])
            react.enzymes.append(mol)
        for m in react_entry['stoichiometry']:
            if m['molecule'] in self.name_mapping.keys():
                m['molecule'] = self.name_mapping[m['molecule']]
            mol = react_element(m['molecule'],m['compartment'],int(float(m['coefficient']))) 
            if float(mol['stoichiometry'])<0 and translation.is_aminoacid(mol['molecule']):
                if 'ligation' in react.annotations:
                    continue
                react.del_aa = mol['molecule']
                continue
            if mol['stoichiometry']>0 and translation.is_aminoacid(mol['molecule']) and 'Deformylation' in react.ref_name:
                react.add_aa = mol['molecule']
                continue
            if mol['stoichiometry']>0 and translation.is_modified_aminoacid(mol['molecule']):
                react.add_aa = mol['molecule']
                continue
            react.molecules.append(mol)
        for m in react_entry['coenzymes']:
            if m['metabolite'] in self.name_mapping.keys():
                m['molecule'] = self.name_mapping[m['metabolite']]
            mol = react_element(m['metabolite'],m['compartment'])
            react.auxiliar_molecules.append(mol)
        react.is_reversible = react_entry['direction'] == 'Reversible'
        if react_entry['kinetics_forward'] != None:
            rate,param = self.get_kinetic_law(react_entry['kinetics_forward'])
            react.rate_law_forward = rate
            react.parameters_forward = param
        if react_entry['kinetics_backward'] != None:
            rate,param = self.get_kinetic_law(react_entry['kinetics_backward'])
            react.rate_law_backward = rate
            react.parameters_backward = param
        if react_entry['optimal_ph']!=None:
            react.optimal_ph = react_entry['optimal_ph']['value']
        if react_entry['optimal_temperature']!=None:
            react.optimal_temperature = " ".join([react_entry['optimal_temperature']['value'],react_entry['optimal_temperature']['units']])
        react.pathways = react_entry['pathways']
        return react
    
    def create_all_protein_modification_reactions(self,translation):
        reactions = self.get_protein_modification_reactions_list()
        react_objs = []
        for i in [self.entry_index[r] for r in reactions]:
            react_objs.append(self.make_protein_modification_reaction(self.get_entry(i),translation))
        return react_objs
    
    def get_aminoacylation_reactions_list(self):
        return list(set([entry['wid'] for entry in self.wcdata['data'] if entry['model']=='Reaction' and entry['processes']=='Process_tRNAAminoacylation' and entry['wid'] not in self.ignore_list]))
    
    def make_aminoacylation_reaction(self,react_entry,translation):
        if react_entry['model']!='Reaction':
            print("ERROR: entry 'model' must be 'Reaction'. Got:",react_entry['model'])
            return False
        if react_entry['processes']!='Process_tRNAAminoacylation':
            print("ERROR: entry 'processes' must be 'Process_tRNAAminoacylation'. Got:",react_entry['processes'])
            return False
        react = ChemicalReaction()
        react.ref_name = react_entry['wid']
        react.name = react_entry['name']
        react.cross_references = {cr["source"]:cr["xid"] for cr in react_entry['cross_references']}
        react.annotations = ['tRNA Aminoacylation']
        for ann in react_entry['type']:
            react.annotations.append(ann)
        react.original_database = self.get_database_name()
        react.process = 'tRNA Aminoacylation'
        if react_entry['enzyme']!=None:
            mol = react_element(react_entry['enzyme']['protein'],react_entry['enzyme']['compartment'])
            react.enzymes.append(mol)
        aminoacid = ''
        if react_entry['type'][0]=='transfer':
            gene_entry = self.find_entry(react_entry['modification']['molecule'])
            aminoacid = gene_entry['amino_acid']
            original_aminoacid = ''        
        for m in react_entry['stoichiometry']:
            mol = react_element(m['molecule'],m['compartment'],int(float(m['coefficient'])))
            if m['molecule']in self.name_mapping.keys():
                mol['molecule'] = self.name_mapping[m['molecule']]
            if translation.is_aminoacid(mol['molecule']):
                if react_entry['type'][0]!='transfer':
                    aminoacid = mol['molecule']
                else:
                    if float(mol['stoichiometry'])<0 and mol['molecule']!=aminoacid:
                        original_aminoacid = mol['molecule']
                        continue
                    elif mol['stoichiometry']>0 and original_aminoacid != mol['molecule']:
                        continue
            react.molecules.append(mol)
        for m in react_entry['coenzymes']:
            mol_name = m['metabolite']
            if m['metabolite'] in self.name_mapping.keys():
                mol_name = self.name_mapping[m['metabolite']]
            mol = react_element(mol_name,m['compartment'])
            react.auxiliar_molecules.append(mol)
        #add aminoacylated tRNA
        aa_tRNA = react_element(react_entry['modification']['molecule']+'_'+aminoacid, react_entry['modification']['compartment'])
        react.molecules.append(aa_tRNA)
        self.created_molecules[aa_tRNA['molecule']] = 'aminoacylated tRNA'
        #add tRNA
        tRNA = {}
        tRNA['molecule'] = react_entry['modification']['molecule']
        if react_entry['type'][0]=='transfer':
            tRNA['molecule'] = tRNA['molecule']+'_'+original_aminoacid
        tRNA['stoichiometry'] = -1
        tRNA['location'] = react_entry['modification']['compartment']
        react.molecules.append(tRNA)
        react.molecules.append(react_element("H2O",'c',-1))
        react.molecules.append(react_element("H",'c'))
        
        react.is_reversible = react_entry['direction'] == 'Reversible'    
        if react_entry['kinetics_forward'] != None:
            rate,param = self.get_kinetic_law(react_entry['kinetics_forward'])
            react.rate_law_forward = rate
            react.parameters_forward = param
        if react_entry['kinetics_backward'] != None:
            rate,param = self.get_kinetic_law(react_entry['kinetics_backward'])
            react.rate_law_backward = rate
            react.parameters_backward = param
        if react_entry['optimal_ph']!=None:
            react.optimal_ph = react_entry['optimal_ph']['value']
        if react_entry['optimal_temperature']!=None:
            react.optimal_temperature = " ".join([react_entry['optimal_temperature']['value'],react_entry['optimal_temperature']['units']])
        react.pathways = react_entry['pathways']
        return react

    
    def create_all_aminoacylation_reactions(self,translation):
        reactions = self.get_aminoacylation_reactions_list()
        react_objs = []
        for i in [self.entry_index[r] for r in reactions]:
            react_objs.append(self.make_aminoacylation_reaction(self.get_entry(i),translation))
        return react_objs
    
    
    def make_protein_monomer(self,prot_entry):
        if prot_entry['model']!='ProteinMonomer':
            print("ERROR: entry 'model' must be 'ProteinMonomer'. Got:",prot_entry['model'])
            return False
        prot = ProteinMonomer()
        prot.ref_name = prot_entry['wid']
        prot.name = prot_entry['name']
        prot.cross_references = {cr["source"]:cr["xid"] for cr in prot_entry['cross_references']}
        prot.annotations = ['Protein Monomer']
        for ann in prot_entry['type']:
            prot.annotations.append(ann)
        prot.original_database = self.get_database_name()
        prot.gene = prot_entry['gene']
        gene = self.find_entry(prot.gene)
        gene_len = 3
        if gene!=None:
            gene_len = int(gene['length'])
        prot.chaperones = prot_entry['chaperones']
        prot.localization = prot_entry['localization'].replace('t','')
        if prot_entry['signal_sequence'] != None:
            prot.signal_sequence = {}
            prot.signal_sequence['length'] = prot_entry['signal_sequence']['length']
            prot.signal_sequence['location'] = prot_entry['signal_sequence']['location']
            prot.signal_sequence['type'] = prot_entry['signal_sequence']['type']
            if prot.signal_sequence['type']=='Lipoprotein':
                prot.modification_reactions.append('ProLipoproteinDiacylglycerylTransferase')
            for i in range(int((gene_len/3)/35)):
                prot.modification_reactions.append('ATPHs_ProteinTranslocation_SecA')                 
        for p in prot_entry['prosthetic_groups']:
            pros = {}
            if p['coefficient'] == 'None':
                p['coefficient'] = 1
            pros['stoichiometry'] = int(float(p['coefficient']))
            pros['location'] = p['compartment']
            pros['molecule'] = p['metabolite']
            prot.prosthetic_groups.append(pros)
        for mod in prot_entry['modification_reactions']:
            if len(mod['reactions'])>0:
                if mod['reactions'][0] not in self.ignore_list:
                    prot.modification_reactions.append(mod['reactions'][0])
        prot.modification_reactions.append('Deformylation')
        if prot_entry['is_n_terminal_methionine_cleaved']['value']=='1':
            prot.modification_reactions.append('NTerminalMethionineClevage')
        return prot
    
    def make_protein_complex(self,compl_entry):
        if compl_entry['model']!='ProteinComplex':
            print("ERROR: entry 'model' must be 'ProteinComplex'. Got:",compl_entry['model'])
            return False
        compl = ProteinComplex()
        compl.ref_name = compl_entry['wid']
        compl.name = compl_entry['name']
        compl.cross_references = {cr["source"]:cr["xid"] for cr in compl_entry['cross_references']}
        compl.annotations = ['Protein Complex']
        for ann in compl_entry['type']:
            compl.annotations.append(ann)
        compl.original_database = self.get_database_name()
        all_molecules = [m['molecule'] for m in compl_entry['biosynthesis']]
        for m in compl_entry['biosynthesis']:
            local_m = m.copy()
            sto = int(float(m['coefficient']))
            if sto<0:
                if m['molecule'] in self.name_mapping.keys():
                    m['molecule'] = self.name_mapping[m['molecule']]
                local_m = m.copy()
                if local_m['molecule'] == "ATP" and ("AMP" in all_molecules or "ADP" in all_molecules):
                    continue
                if local_m['molecule'] == "COA":
                    local_m['molecule'] = "phosphopantetheine"
                compl.composition = compl.composition + [local_m['molecule'] for i in range(-sto)]
        if compl.ref_name == 'MG_224_9MER_GDP':
            compl.composition = ['MG_224_MONOMER_GDP' for i in range(9)]
        return compl
        
    def make_stimulus(self,sti_entry):
        if sti_entry['model']!='Stimulus':
            print("ERROR: entry 'model' must be 'Stimulus'. Got:",sti_entry['model'])
            return False
        sti = Stimulus()
        sti.ref_name = sti_entry['wid']
        sti.name = sti_entry['name']
        sti.cross_references = {cr["source"]:cr["xid"] for cr in sti_entry['cross_references']}
        sti.annotations = ['Stimulus']
        for ann in sti_entry['type']:
            sti.annotations.append(ann)
        sti.original_database = self.get_database_name()
        sti.value = sti_entry['value']['value']+' '+sti_entry['value']['units']
        return sti
    
    def make_aminoacylated_tRNA(self,mol):
        compl = RNAComplex()
        compl.ref_name = mol
        compl.name = "Aminoacylated tRNA "+mol
        compl.annotations = ['Aminoacylated tRNA']
        compl.original_database = "Created by WCKB handler"
        aa = mol.split('_')[-1]
        rna = mol.replace("_"+aa,"")
        compl.composition.append(aa)
        compl.composition.append(rna)
        return compl
    
    def make_rna(self,rna_entry):
        if rna_entry['model']!='Gene':
            print("ERROR: entry 'model' must be 'Gene'. Got:",rna_entry['model'])
            return False
        rna = RNA()
        rna.ref_name = rna_entry['wid']
        rna.name = rna_entry['name']+' RNA'
        rna.symbol = rna_entry['symbol']+' RNA'
        rna.gene = rna_entry['wid']
        rna.annotations = ['RNA']
        for ann in rna_entry['type']:
            rna.annotations.append(ann)
        rna.transcription_unit = rna_entry['transcription_units'][0]
        #rna.type = rna_entry['type'][0]
        rna.original_database = self.get_database_name()
        #for mod in rna_entry['modification_reactions']:
        #    if len(mod['reactions'])>0:
        #        for r in mod['reactions']:
        #            if 'MODIFICATION' in r:
        #                rna.modification_reactions.append(r)
        return rna
        
    def make_molecule_from_entry_name(self,mol):
        mol_obj = None
        entry = self.find_entry(mol)
        if entry!=None:
            if entry['wid'] in self.name_mapping.keys():
                entry = self.find_entry(self.name_mapping[entry['wid']])
            model = entry['model']
            if model == 'ProteinComplex':
                mol_obj = self.make_protein_complex(entry)
            elif model == 'ProteinMonomer':
                mol_obj = self.make_protein_monomer(entry)
            elif model == 'Stimulus':
                mol_obj = self.make_stimulus(entry)
            elif model == 'Metabolite':
                mol_obj = self.make_metabolite(entry)
            elif model == 'Gene':
                mol_obj = self.make_rna(entry)
        elif mol in self.created_molecules.keys():
            if self.created_molecules[mol] == 'aminoacylated tRNA':
                mol_obj = self.make_aminoacylated_tRNA(mol)
            
        return mol_obj
            
    
    def get_ribosome_assembly_reactions_list(self):
        return list(set([entry['wid'] for entry in self.wcdata['data'] if entry['model']=='Reaction' and entry['processes']=='Process_RibosomeAssembly' and entry['wid'] not in self.ignore_list]))
        
    def make_complex_biosynthesis_reaction(self,compl_entry):
        if compl_entry['model']!='ProteinComplex':
            print("ERROR: entry 'model' must be 'ProteinComplex'. Got:",compl_entry['model'])
            return False
        compl = ComplexBiosynthesisReaction()
        compl.ref_name = compl_entry['wid']+'_BIOSYNTHESIS'
        compl.name = compl_entry['name']+' Biosynthesis'
        compl.cross_references = {cr["source"]:cr["xid"] for cr in compl_entry['cross_references']}
        compl.annotations = ['Protein Complex Biosynthesis']
        compl.process = 'Macromolecular Complexation'
        for ann in compl_entry['type']:
            compl.annotations.append(ann)
        compl.original_database = self.get_database_name()
        product_counter = 0
        for m in compl_entry['biosynthesis']:
            mol = react_element(m['molecule'],m['compartment'],int(float(m['coefficient'])))
            if m['molecule'] in self.name_mapping.keys():
                mol['molecule'] = self.name_mapping[m['molecule']]
            compl.molecules.append(mol)            
            if mol['molecule'] == "MG_287_MONOMER" and mol['stoichiometry']<0:
                compl.molecules.append(react_element("H2O",'c',-2))
            if mol['molecule'] == "MG_287_MONOMER_ACP" and mol['stoichiometry']>0:
                compl.molecules.append(react_element("H2O",'c',1))
                compl.molecules.append(react_element("H",'c',-1))
            product_counter += int(mol['stoichiometry']>0)
        if "RIBOSOME" in compl_entry['wid'] and not "70" in compl_entry['wid']:
            rib_assembly_react = self.get_ribosome_assembly_reactions_list()
            for r in rib_assembly_react:
                react_entry = self.find_entry(r)
                if react_entry['enzyme']!=None:
                    mol = react_element(react_entry['enzyme']['protein'],react_entry['enzyme']['compartment'])
                    compl.enzymes.append(mol)
                for m in react_entry['stoichiometry']:
                    mol = react_element(m['molecule'],m['compartment'],int(float(m['coefficient'])))
                    if m['molecule'] in self.name_mapping.keys():
                        mol['molecule'] = self.name_mapping[m['molecule']]
                    compl.molecules.append(mol)
                    product_counter += int(mol['stoichiometry']>0)
        compl.is_reversible = product_counter == 1
        if 'ox' in compl.ref_name:
            compl.is_reversible = True
        return compl
        
        
    def make_biosynthesis_for_molecule(self,mol):
        mol_obj = None
        entry = self.find_entry(mol)
        if entry!=None:
            model = entry['model']
            if model == 'ProteinComplex':
                mol_obj = self.make_complex_biosynthesis_reaction(entry)
        return mol_obj            
        
        
def generate_sbml_model(g,path):
    import libsbml as sbml
    
    try:
        document = sbml.SBMLDocument(3, 1)
    except ValueError:
        raise SystemExit('Could not create SBMLDocumention object')

 
    model = document.createModel()
    model.setTimeUnits("second")
    
    compartments = list(set([x.split('_')[0] for x,y in g.nodes(data=True) if y['type']=='m']))
    
    for comp in compartments:
        c1 = model.createCompartment()   
        c1.setId(comp)
        c1.setConstant(True)
        c1.setSize(1)
        c1.setSpatialDimensions(3)
        c1.setUnits('litre')
    
    
    molecules = [x for x,y in g.nodes(data=True) if y['type']=='m']
    for m in molecules:
        s1 = model.createSpecies()
        s1.setId(m)
        compartment = m.split('_')[0]
        s1.setCompartment(compartment)
        s1.setConstant(False)
        s1.setBoundaryCondition(False)
        s1.setHasOnlySubstanceUnits(True)
    
    
    reactions = [x for x,y in g.nodes(data=True) if y['type']=='r']
    
    for react in reactions:
        r1 = model.createReaction()
        r1.setId(react)
        r1.setReversible(False)
        r1.setFast(False)
 
        reactants = g.predecessors(react)
        for m in reactants:
            sto = g[m][react]['sto']
            t = g[m][react]['type']
            if t == 'r':
                sp = r1.createReactant()
                sp.setConstant(True)
                sp.setStoichiometry(sto)
            elif t == 'm':
                sp = r1.createModifier()
            sp.setSpecies(m)            
            
            
        products = g.successors(react)
        for m in products:
            sto = g[react][m]['sto']
            sp = r1.createProduct()
            sp.setSpecies(m)
            sp.setConstant(True)
            sp.setStoichiometry(sto)
            
    with open(path,'w') as f:
        f.write(sbml.writeSBMLToString(document))