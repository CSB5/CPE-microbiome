from pymol import cmd
import pandas as pd
from tqdm import tqdm
import numpy as np
from distancetoatom_script import distancetoatom

def get_atom_from_selection(sele):
    """
    sele is selection string
    """
    # Lok through selection using either cmd.iterate vs cmd.get_model
    myspace = {'atomID': []}
    cmd.iterate(sele, 'atomID.append([ID,name,resn,resi,chain])', space=myspace)
    atomDF = pd.DataFrame(myspace['atomID'],columns=['ID','name','resn','resi','chain'])
    atomDF = atomDF.astype({'resi':'int32','ID':'int32'})
    return atomDF

def find_dist(query, target, distance=10):
    """
    query and target are selection criteria
    E.g. find_dist('chain A', 'chain B', 20)
    """
    cmd.select('all_atoms','all')
    cmd.select("query", query)
    cmd.select("target", target)
    cmd.select("query2target", f'query near_to {distance} of target')
    try:
        atom_count = cmd.select("target2query", f'target near_to {distance} of query')
        if atom_count==0:
            raise Exception
    except:
        atom_count = cmd.select("target2query", 'target')
    atomDF = get_atom_from_selection('all_atoms')
    atomQueryPrefilterDF = get_atom_from_selection('query2target')
    dist2target = []
    for id in tqdm(atomDF['ID']):
        # check if atom has coordinate associated with it
        try:
            if cmd.count_atoms(f"id {id}") == 1:
                # atom coordinates
                cmd.get_atom_coords(f"id {id}")
        except:
            dist2target.append(np.nan)
            continue
        # check if atom in query near target. If not, nan
        if (atomQueryPrefilterDF['ID']==id).sum()==0:
            dist2target.append(np.nan)
            continue
        # Calculate min distance from closest atom in dimer to closest atom in DNA
        dst = distancetoatom(
            origin=f"id {id}",
            cutoff=f'{distance}',
            selection='target2query',
        )
        dstDf = pd.DataFrame(dst[2:],columns=dst[0])
        dist2target.append(dstDf['distance_to_origin'].min())
    return dist2target
