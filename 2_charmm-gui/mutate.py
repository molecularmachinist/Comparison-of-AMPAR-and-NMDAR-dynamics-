import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] 

import sys, time, os
import itertools
import pymol
import lil_helpers as lh


def mutate_protein_pdb(pdbIn,pdbOut,arr_iRes,arr_mutRes,arr_chain=['A']):
    """ Mutates pdb anf outputs pdb """
    
    
    
    
    
    # Name of the pdb file ( without path and .pdb)
    pdbName = pdbIn.split('/')[-1].split('.')[0]
    
    # Load the structure
    pymol.cmd.load(pdbIn, pdbName)
    pymol.cmd.wizard('mutagenesis')
    
    for iRes, mutRes, chain in zip(arr_iRes,arr_mutRes,arr_chain):
        print(iRes, mutRes, chain)
        
        # Minty fresh
        pymol.cmd.do('refresh_wizard')

        # Set target residue type
        pymol.cmd.get_wizard().set_mode(str(mutRes))

        # Select residue
        selection = '/{}//{}/{}'.format(pdbName,chain,iRes)

        # Notify the wizard about the selection
        pymol.cmd.get_wizard().do_select(selection)

        # Select frame 1
        pymol.cmd.frame('1')

        # Apply the mutation
        pymol.cmd.get_wizard().apply()
    

    # Create
    pymol.cmd.save(pdbOut)
    
    # and destroy
    pymol.cmd.delete('all')


pdb = '../Add_missing_residues_1/AMPA/5weo_cleaned_filled_pymol_aligned_glutamate.pdb'
ABCD = lh.seq_in_chains(pdb,chains='ABCD')




i0 = 1
arr_iRes,arr_mutRes,arr_chain =[], [], []

for chain in 'ABCD':
    hih = lh.diff_seq(ABCD[chain],lh.seq_ampa)
    
    for iRes, _, mutRes in hih:
        
        arr_iRes.append(iRes+i0)
        arr_mutRes.append(lh.one2three[mutRes])
        arr_chain.append(chain)
    
    i0 += len(ABCD[chain])

pdb_out = '5weo_cleaned_filled_pymol_aligned_glutamate_mutated.pdb'


mutate_protein_pdb(pdb,pdb_out,arr_iRes,arr_mutRes,arr_chain)