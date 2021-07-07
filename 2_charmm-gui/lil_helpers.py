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