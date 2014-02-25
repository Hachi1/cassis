import re
from astropy import units as u, constants as const
from cassis.pyparsing import Word, Literal, Group, Keyword, OneOrMore, Optional, Combine, alphas, nums, alphanums
from collections import OrderedDict
import pandas as pd


plusorminus = Literal( '+' ) | Literal( '-' )
#decimal pattern = float number, optional +/- beginning, optional "," thousands separators within
decimal_pattern = Combine(Optional(plusorminus) + Word(nums, nums+",") + Optional("." + OneOrMore(Word(nums))))
#gives string representation, gives float representation
decimal_pattern.setParseAction(lambda token: ''.join(token), lambda token: float(token[0]))

int_pattern = Word(nums).addParseAction(lambda token: int(token[0]))

#atomic number / ion number pattern == species pattern -- gives integer atom/ion numbers from 1-3 digit numbers.
atomic_number_pattern = Word(nums, min=1, max=3).setResultsName('atomic_number').setParseAction(lambda i: int(i[0]) )
ion_number_pattern = Word(nums, min=1, max=3).setResultsName('ion_number').setParseAction(lambda i: int(i[0]) )
#suprresses parsing of parts
species_pattern = atomic_number_pattern + Literal('.').suppress() + ion_number_pattern

header_pattern = species_pattern + Word(nums).setResultsName('lines_saved').setParseAction(lambda i: int(i[0])) + Literal('lines saved').suppress()
header_pattern += Word(nums).setResultsName('positive_lines_saved').setParseAction(lambda i: int(i[0])) + Literal('positive lines saved').suppress()
header_pattern += Word(nums).setResultsName('even_lines').setParseAction(lambda i: int(i[0])) + Literal('even').suppress()
header_pattern += Word(nums).setResultsName('odd_lines').setParseAction(lambda i: int(i[0])) + Literal('odd levels').suppress()
header_pattern += Word(nums).setResultsName('ionization_potential').setParseAction(lambda i: ((float(i[0])/u.cm) * const.c * const.h).to('eV')  )  + Keyword('ion pot cm-1 eve').suppress()

term_pattern = Group(OneOrMore(Word(alphanums, min=2))).setResultsName('term')
term_label_pattern = Word(alphanums, min=1, max=1).setResultsName('label') + term_pattern

level_energy = decimal_pattern.copy()
level_energy.addParseAction(lambda token: (float(token[0])/u.cm * const.c * const.h).to('eV'))

level_data_pattern = species_pattern + (Literal('EVE') ^ Literal('ODD')).setResultsName('parity') 
level_data_pattern += int_pattern.setResultsName('index') 
level_data_pattern += decimal_pattern.setResultsName('energy') 
level_data_pattern += decimal_pattern.setResultsName('J')


def read_gf_gam_file(fname):
    """
    Reading a gfxxxx.gam file
    
    Parameters:
    -----------
    
    fname: str
        filename to read in 
        
    Returns:
    --------
    
    ~dict containing the header
    ~pandas.Dataframe object containing the tabular level data 
    ~pandas.Dataframe object containing the tabular level data second part (nobody knows what it is ;) )
    
    """
    fhandle = open(fname)

    #reading the header
    linestr=fhandle.readline()         
    header = header_pattern.parseString(linestr).asDict()

    print 'header', header

    #reading the term label
    while True:
        linestr=fhandle.readline()    
        if linestr.strip().startswith('level'): break
        
    while True:
        linestr=fhandle.readline()    
        if linestr.strip().startswith('level'): continue
        if linestr.strip().startswith('ELEM'): break
    
    
    levels_data = OrderedDict( [  ('atomic_number',[]), ('ion_number',[]), ('parity',[]), ('index',[]), ('energy',[]), ('J',[]), ('levdescriptor',[]), \
                                  ('glande',[]), ('suma',[]), ('c4',[]), ('c6',[]), ('sumf',[]), ('lec_1_coeff',[]), ('lec_1_label',[]), ('lec_1_counter',[]), \
                                  ('lec_2_coeff',[]), ('lec_2_label',[]), ('lec_2_counter',[]), ('lec_3_coeff',[]), ('lec_3_label',[]), ('lec_3_counter',[])  ] )

    levels_data_secpart = OrderedDict( [  ('atomic_number',[]), ('ion_number',[]), ('parity',[]), ('index',[]), ('energy',[]), ('J',[]), ('levdescriptor',[]), \
                                  ('unkown_coeff',[])  ] )
     
    #first part of table
    while True:
            linestr=fhandle.readline()
            if len(linestr)<90: break
            #use parser for all fields up to J
            parsed_level_data_dict = level_data_pattern.parseString(linestr).asDict()
            for key in parsed_level_data_dict: levels_data[key].append(parsed_level_data_dict[key])
            levels_data['levdescriptor'].append(linestr[30:41].strip())
            levels_data['glande'].append(float(linestr[42:47]))
            levels_data['suma'].append(float(linestr[47:56]))
            levels_data['c4'].append(float(linestr[56:65]))
            levels_data['c6'].append(float(linestr[65:74]))
            levels_data['sumf'].append(float(linestr[74:80]))
            levels_data['lec_1_coeff'].append(float(linestr[81:87]))
            levels_data['lec_1_label'].append(linestr[88:96].strip())
            levels_data['lec_1_counter'].append(int(linestr[96:98]))
            #the following may or may not be there
            if len(linestr[98:104].strip())>0:
                    levels_data['lec_2_coeff'].append(float(linestr[98:104]))
            else:
                    levels_data['lec_2_coeff'].append(0.0)
            levels_data['lec_2_label'].append(linestr[105:113])
            if len(linestr[113:115].strip())>0:
                    levels_data['lec_2_counter'].append(int(linestr[113:115]))
            else:
                    levels_data['lec_2_counter'].append(0)
            if len(linestr[115:121].strip())>0:
                    levels_data['lec_3_coeff'].append(float(linestr[115:121]))
            else:
                    levels_data['lec_3_coeff'].append(0.0)
            levels_data['lec_3_label'].append(linestr[122:130])
            if len(linestr[130:132].strip())>0:
                    levels_data['lec_3_counter'].append(int(linestr[130:132]))
            else:
                    levels_data['lec_3_counter'].append(0)

            
    #second part of table where line size ~ 80
    while True:
            linestr=fhandle.readline()
            if linestr=="": break
            #use parser for all fields up to J
            parsed_level_data_dict = level_data_pattern.parseString(linestr).asDict()
            for key in parsed_level_data_dict: levels_data_secpart[key].append(parsed_level_data_dict[key])
            levels_data_secpart['levdescriptor'].append(linestr[29:43].strip())
            levels_data_secpart['unkown_coeff'].append(float(linestr[43:64]))
 
    
    levels_data = pd.DataFrame(levels_data)
    levels_data['observed'] = levels_data['energy'] >= 0
    levels_data['energy'] = abs((levels_data['energy'].values/u.cm * const.c * const.h).to('eV')).value
  
    levels_data_secpart = pd.DataFrame(levels_data_secpart)
    levels_data_secpart['observed'] = levels_data_secpart['energy'] >= 0
    levels_data_secpart['energy'] = abs((levels_data_secpart['energy'].values/u.cm * const.c * const.h).to('eV')).value
    
    return header, levels_data, levels_data_secpart
            
