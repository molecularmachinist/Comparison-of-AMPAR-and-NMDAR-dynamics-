from modeller import *
from modeller.automodel import *
# sequence from uniprot

seq_ampa_full='MQKIMHISVLLSPVLWGLIFGVSSNSIQIGGLFPRGADQEYSAFRVGMVQFSTSEFRLTPHIDNLEVANSFAVTNAFCSQFSRGVYAIFGFYDKKSVNTITSFCGTLHVSFITPSFPTDGTHPFVIQMRPDLKGALLSLIEYYQWDKFAYLYDSDRGLSTLQAVLDSAAEKKWQVTAINVGNINNDKKDEMYRSLFQDLELKKERRVILDCERDKVNDIVDQVITIGKHVKGYHYIIANLGFTDGDLLKIQFGGANVSGFQIVDYDDSLVSKFIERWSTLEEKEYPGAHTTTIKYTSALTYDAVQVMTEAFRNLRKQRIEISRRGNAGDCLANPAVPWGQGVEIERALKQVQVEGLSGNIKFDQNGKRINYTINIMELKTNGPRKIGYWSEVDKMVVTLTELPSGNDTSGLENKTVVVTTILESPYVMMKKNHEMLEGNERYEGYCVDLAAEIAKHCGFKYKLTIVGDGKYGARDADTKIWNGMVGELVYGKADIAIAPLTITLVREEVIDFSKPFMSLGISIMIKKPQKSKPGVFSFLDPLAYEIWMCIVFAYIGVSVVLFLVSRFSPYEWHTEEFEDGRETQSSESTNEFGIFNSLWFSLGAFMQQGCDISPRSLSGRIVGGVWWFFTLIIISSYTANLAAFLTVERMVSPIESAEDLSKQTEIAYGTLDSGSTKEFFRRSKIAVFDKMWTYMRSAEPSVFVRTTAEGVARVRKSKGKYAYLLESTMNEYIEQRKPCDTMKVGGNLDSKGYGIATPKGSSLRNAVNLAVLKLNEQGLLDKLKNKWWYDKGECGSGGGDSKEKTSALSLSNVAGVFYILVGGLGLAMLVALIEFCYKSRAEAKRMKVAKNAQNINPSSSQNSQNFATYKEGYNVYGIESVKI'


seq_ampa = seq_ampa_full[24:841]         # correct end and begining
seq_ampa =seq_ampa[:374]+seq_ampa[380:]  # remove deleted residues from the middle


three2one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

one2three = {v: k for k, v in three2one.items()}

def seq_in_chains(pdb, chains=['A','B','C','D']):
    result = dict(zip(chains, ['']*len(chains)))
    with open(pdb,'r') as file:
        res_id_old = 0
        for line in file:
            if line[:4] == 'ATOM':
                res, chain,res_id = line.split()[3:6]
                if len(chain) > 1:
                    res_id = chain[1:]
                    chain = chain[0]
                if chain in chains:
                    res_id = int(res_id)
                    if res_id_old < res_id or result[chain]=='':
                        res_id_old = res_id
                        result[chain] += three2one[res]
    return result

def diff_seq(seq1,seq2):
    
    diff = []
    for i, aa_pair in enumerate( zip(seq1,seq2) ):
        if aa_pair[0] != aa_pair[1]:
            diff.append([i,aa_pair[0],aa_pair[1]])
    return diff

def print_diff_seq(seq1, seq2, i0=1):
    print('Resid\tseq1       seq2')
    for i, aa1, aa2 in diff_seq(seq1,seq2):
        print('{},\t{}({}) !=  {}({})'.format(i+i0,one2three[aa1],aa1,one2three[aa2],aa2))


def seq_from_modeller(code_in, seq_file_out):
    
    e = Environ()
    m = Model(e, file=code_in)
    aln = Alignment(e)
    aln.append_model(m, align_codes=code_in)
    aln.write(file=seq_file_out)


    
def aligment_file(full_seqs, missing_seqs, modeller_seq_file_in, file_out):
    
    # read header from modeller .seq file
    with open(modeller_seq_file_in, 'r') as file:
        header = [ next(file) for _ in range(3)]
    
    seq_missing = ''
    for chain in missing_seqs:
        seq_missing += chain + '/'
    seq_missing += '*'
    
    
    seq_full = ''
    for chain in full_seqs:
        seq_full += chain + '/'
    seq_full += '*'
    
    
    
    # Formate for the file
    seq_missing_for_file = ''
    seq_full_for_file = ''
    
    start, end, length = 0, 75, len(seq_missing)
    
    if length != len(seq_full):
        raise ValueError('All chains of full_seqs and missing_seqs should have same number of amino acids')
    
    while end < length + 75:
        seq_missing_for_file += seq_missing[start:end]+'\n'
        seq_full_for_file += seq_full[start:end]+'\n'
        start += 75
        end += 75
        
    with open(file_out, 'w') as file:
        for line in header:
            file.write(line)
            
        file.write(seq_missing_for_file)
        
        file.write(header[1][:-1]+'_fill\nsequence:::::::::\n')
        
        file.write(seq_full_for_file)



def add_residues_with_modeller(code,alnfile, path='.'):
    log.verbose()
    env = Environ()

    # directories for input atom files
    env.io.atom_files_directory = [path]


    a = AutoModel(env, alnfile = alnfile,
                knowns = code, sequence = code + '_fill')
    a.starting_model= 1
    a.ending_model  = 1

    a.make()
    