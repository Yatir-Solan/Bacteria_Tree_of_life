# modules import:
import io
import re
import os
import sys
import random
import pickle
import itertools
import collections
import pandas as pd
import subprocess as sp
from Bio import SeqIO,AlignIO
from ete3 import NCBITaxa,Tree
from Bio.Align.Applications import MafftCommandline
# python scripts import:
sys.path.append('/davidb/yatirsolan/scripts/python/bio_utilities')
sys.path.append('/davidb/yatirsolan/scripts/python/cluster_utilities')
import computational_tools
import phylogenetics
from phylogenetics import name_to_taxid, taxid_to_name
import fasta_import
import databases
import HMM
#
########################################################

def regs():
    return {'annot':re.compile(r'(?<=\s).+(?= OS=)'),
           'minim_desc':re.compile(r'OS=.+(?= GN=)'),
           'tax':re.compile(r'(?<=OS=).+(?= OX=)'),
           'id':re.compile(r'(?<=OX=)\d+\S')}

########################################################

def hmmsearch_filter(hmmsdf):
    hmmsdf.sort_values(by=['E-value'], ascending=[True], inplace=True) # all values are now sorted by their E-values. significant -> less-significant. 
    # all the drop_duplicate lines to come (in this function) are based on this sort.
    hmmsdf['queryName'] = hmmsdf['queryName'].apply(lambda x:x.split('.')[0]) # K02874.1 -> K02874 (all Hmms targets of the same ribosomal protein, are named the same).
    hmmsdf.drop_duplicates(subset=['queryName','targetName'], keep='first', inplace=True) # this line here to keep only line number 798 from the example given bellow. 

    # it keeps the first apperance because the file originally is sorted by evalue. "first" is the default and it is here only for emphasising its importance.
    # 
    #       targetName                       queryName E-value             tmp_targetName
    # 798   tr|A0A4V3RXN9|A0A4V3RXN9_9PROT    K02874  6.900000e-66             tr
    # 7994  tr|A0A4V3RXN9|A0A4V3RXN9_9PROT    K02874  6.900000e-13             tr
    # 8072  tr|A0A4V3RXN9|A0A4V3RXN9_9PROT    K02874  1.600000e-13             tr
    #

    hmmsdf['OX'] = hmmsdf['targetDescription'].apply(lambda x:regs().get('id').search(x).group()) # creates a column for the species taxid of the matched protein.
    hmmsdf['OS'] = hmmsdf['targetDescription'].apply(lambda x:regs().get('tax').search(x).group()) # creates a column for the species (in text) of the matched protein.
    hmmsdf.drop_duplicates(subset=['targetName'], keep='first', inplace=True) # deletes double occurances of targetNames that were mapped to different ribosomal proteins.
    hmmsdf.drop_duplicates(subset=['queryName','OX'], keep='first', inplace=True) # deletes double occurances of different orf that correspond to the same hmm. 

    hmm_hits_per_taxid = {taxid:sub_df.shape[0] for taxid, sub_df in hmmsdf.groupby(by='OX')} # the groupby as an iterator returns a tuple -> (sort_by_variable, sub_dataframe)
    taxids_to_filter = dict(filter(lambda elem: elem[1] < 15, hmm_hits_per_taxid.items())) # find all the taxids of which are abscent of 15 mathces. (15 is the number or ribosomal proteins that in search in first place...)
    hmmsdf = hmmsdf[~hmmsdf.OX.isin(taxids_to_filter.keys())] # filtering the main df - excluding all taxids in taxids_to_filter (which are missing some ribosomal proteins matches)

    return hmmsdf

########################################################

def hmmsdf_rank_annotation(hmmsdf, rank, taxid_map):
    hmmsdf[rank] = hmmsdf.loc[:,'OX'].apply(lambda rnk_txid:taxid_map.get(rnk_txid,dict()).get(rank)) # dict() is here, because if rnk_txid is not found in taxid_map a empty dictionary will be...
    # returned so None will be finally established, ready to be filtered in the necxt line...
    hmmsdf = hmmsdf[~(hmmsdf[rank].isnull())] # filters all family/order None values - in case they were not found within taxid_map.
    # hmmsdf = hmmsdf[hmmsdf[rank]!='None']
    return hmmsdf

########################################################

def candidatus_filter(hmmsdf, rank, ncbi, taxid_map):
    hmmsdf['phylum_name'] = hmmsdf[rank].apply(lambda rnk_txd:taxid_map.get(rnk_txd).get('phylum'))
    hmmsdf = hmmsdf[~(hmmsdf['phylum_name'].isnull())] # filters all phylums None values - in case they were not found within taxid_map.
    hmmsdf['phylum_name'] = hmmsdf['phylum_name'].apply(lambda phylum_txd:taxid_to_name(phylum_txd, ncbi))
    hmmsdf[f'{rank}_name'] = hmmsdf[rank].apply(lambda rnk_txd: taxid_to_name(rnk_txd, ncbi))
    hmmsdf = hmmsdf[~(hmmsdf['phylum_name'].isnull())] # filters all phylums None values - in case they were not found within taxid_map.
    hmmsdf = hmmsdf[~(hmmsdf[f'{rank}_name'].isnull())] # filters all phylums None values - in case they were not found within taxid_map.
    hmmsdf = hmmsdf[~(hmmsdf[f'phylum_name'].str.contains('Candidatus'))] # filters Candidtuses in the phylum rank dgree.
    hmmsdf = hmmsdf[~(hmmsdf[f'{rank}_name'].str.contains('Candidatus'))] # filters Candidtuses in the family/order rank dgree.
    hmmsdf = hmmsdf[~(hmmsdf.OS.str.contains('Candidatus'))] # filters Candidtuses in the rank (species/strain - I believe...) dgree.
    hmmsdf.drop(columns=['phylum_name', f'{rank}_name'], axis=1, inplace=True)

    return hmmsdf

########################################################

def get_headers(hmmsdf, rank, ncbi, offset_leaves=None):
    headers_coordinates = {rnk:(0,15) for rnk in set(hmmsdf.loc[:,rank])}
    if offset_leaves:
        for leaf in offset_leaves:
            next_start_idx = 15
            next_end_idx = 30
            leaf_txid = name_to_taxid(leaf, ncbi) # all txids are handled as strings.
            if hmmsdf[hmmsdf.loc[:,rank] == leaf_txid].sort_values(by=['OX'], ascending=[True]).iloc[next_start_idx:next_end_idx:].empty: 
                # inside this if statement -> means that the dataframe is empty, and that means that, as far as my analyze is correct...
                # ... there is no taxonomical valid representation of the the family/order in refseq dataset. 
                hmmsdf = hmmsdf[hmmsdf.loc[:,rank]!=leaf_txid]
                headers_coordinates[leaf] = 'not-valid'
            else:
                headers_coordinates[leaf] = (next_start_idx, next_end_idx)
    
    rib_hdrs_dic = {rib:{rnk:None for rnk in set(hmmsdf.loc[:,rank])} for rib in set(hmmsdf.queryName)} # {K02874:{12358:None, 5491:None,...}, K02342:{12358:None, 5491:None,...}}

    for rnk, rnk_df in hmmsdf.groupby(by=rank): # the groupby as an iterator returns a tuple -> (sort_by_variable, sub_dataframe)
        start_idx = headers_coordinates.get(rnk)[0]
        end_idx = headers_coordinates.get(rnk)[1]
        for rib,hdr in tuple(rnk_df.sort_values(by=['OX'], ascending=[True]).iloc[start_idx:end_idx:,].loc[:,['queryName','targetName']].values):
            rib_hdrs_dic.get(rib)[rnk] = hdr # e.g -> {K02874:{12358:tr|A0A5R9QKH3|A0A5R9QKH3_9BACT, 5491:sp|Q2YRA5|RL14_BRUA2,...}, K02342:{12358:sp|C3MAZ0|RL14_SINFN, 5491:tr|A0A316J933|A0A316J933_9RHIZ,...}}

    rib_hdrs_dic = {rib:dict(collections.OrderedDict(sorted(hdrs.items()))) for rib,hdrs in rib_hdrs_dic.items()} # sotring the hdrs by rank taxids
    rib_hdrs_dic = {rib:list(hdrs.values()) for rib,hdrs in rib_hdrs_dic.items()} 
    # after the sotring is done, excluding the ranks, leaving only the hdrs as list in the vales of the dic : e.g -> {K02874:[tr|A0A5R9QKH3|A0A5R9QKH3_9BACT, sp|Q2YRA5|RL14_BRUA2, ...], K02342:[tr|..., sp|..., ]}
    return rib_hdrs_dic

########################################################

def differential_alignments(rib_fasta_dic, directory):
    # the functions returnes a list of MSA objects.
    print("differential_alignments")
    MSAs_lst = list()
    for rib,rcrds in rib_fasta_dic.items():
        fasta_file = os.path.join(directory,f'{rib}.faa')
        SeqIO.write(rcrds, fasta_file, 'fasta')
        mafft_cline = MafftCommandline(cmd=computational_tools.mafft_paths(algorithm='linsi'), input=fasta_file)
        stdout,stderr = mafft_cline()
        alignment = AlignIO.read(io.StringIO(stdout), 'fasta')
        MSAs_lst.append(alignment)
        os.remove(fasta_file)
    return MSAs_lst

########################################################

def alignments_concatenation(MSAs_lst, rank, taxid_map, ncbi):
    # the functions returnes an MSA objects (the single concatenated one - built out of all the MSAs together).
    print("alignments_concatenation")
    for MSA in MSAs_lst:
        try:
            cncatntd_MSA += MSA
        except NameError:
            cncatntd_MSA = MSA

    header_names = list()
    for rec in MSAs_lst[0]:
        desc = rec.description
        header_names.append({'taxid':taxid_map.get(regs().get('id').search(desc).group(),dict()).get(rank),
                             'family_name':taxid_to_name(taxid_map.get(regs().get('id').search(desc).group(),dict()).get(rank),ncbi),
                             'description':regs().get('minim_desc').search(desc).group()})

    for hdr,names in zip(cncatntd_MSA, header_names):
        hdr.id = names.get('family_name').replace(' ','_') # if there are spaces there, the fasta file 'thinks' it's where the description starts...
        hdr.name = hdr.id
        hdr.description = names.get('description')
    return cncatntd_MSA

########################################################

def phyl_tree(alignment_file, directory):
    tree_file = os.path.join(directory, f'{os.path.splitext(alignment_file)[0]}.nw')
    sp.run(f"{computational_tools.tree_path(algorithm='fast_tree')} {alignment_file} > {tree_file}",shell=True,capture_output=True,text=True) # bash command to run directly within the lab cluster
    return tree_file

########################################################

def tree_construction_pipeline(hmmsdf, rank, target, directory, taxid_map, ncbi, filter_leaves, offset_leaves=None):
    rib_hdrs_dic = get_headers(hmmsdf, rank, ncbi, offset_leaves)
    rib_fasta_dic = {rib:fasta_import.pullseq(headers=hdrs, target=target) for rib,hdrs in rib_hdrs_dic.items()}
    MSAs_lst = differential_alignments(rib_fasta_dic, directory)
    cncatntd_MSA = alignments_concatenation(MSAs_lst, rank, taxid_map, ncbi)
    msa_file = os.path.join(directory,f"review_{rank}{'_filtered' if filter_leaves else str()}.aln")
    AlignIO.write(cncatntd_MSA, msa_file, 'fasta')
    tree_file = phyl_tree(msa_file, directory)
    return tree_file,cncatntd_MSA

########################################################

def leaf_verifier(tree_file, rank, taxid_map, ncbi, directory):
    # The function aimes to target false taxonomical annotations within the RefSeq data.
    def error_clades_finder(tree, taxid_map, ncbi):
        error_clades_lst = list()
        for leaf in Tree(tree, format=1):
            upper_node = leaf.up
            leaves = [leaf.name for leaf in upper_node.get_leaves()] # e.g. -> ['Sphaerobacteraceae', 'Herpetosiphonaceae', 'Roseiflexaceae', 'Oscillochloridaceae', 'Chloroflexaceae']
            leaves_taxids = [name_to_taxid(leaf.replace('_',' '), ncbi) for leaf in leaves] # e.g. -> ['85002', '189773', '1508635', '104174', '1106']
            phylums = [taxid_to_name(taxid_map.get(v).get('phylum'), ncbi) for v in leaves_taxids] # e.g. -> ['Chloroflexi', 'Chloroflexi', 'Chloroflexi', 'Chloroflexi', 'Chloroflexi']
            if len(set(phylums)) > 1: # when leafs in a clade are members of different phyla, it is concluded into the error clades list.
                error_clades_lst.append(dict(zip(leaves, phylums))) # e.g. -> [{'Brachyspiraceae': 'Spirochaetes', 'Gloeobacteraceae': 'Cyanobacteria', 'Gloeomargaritaceae': 'Cyanobacteria'}, {...}, ...]
        error_clades_lst = sorted(error_clades_lst, key=len, reverse=True) # this sort is probabley not crutial for the script... but it is more elegant.
        return error_clades_lst

    def subclades_filter(error_clades_lst):
        sublist_to_remove = list()
        for clade in error_clades_lst:
            for subset_clade in error_clades_lst:
                if clade == subset_clade:
                    continue
                if subset_clade.items() <= clade.items(): # the operator '<=' checks if the left vairable is a subset of the right one.
                    sublist_to_remove.append(subset_clade)
        return filter(lambda x:x not in sublist_to_remove, error_clades_lst)

    error_ranks = list()
    error_clades_lst = error_clades_finder(tree_file,taxid_map,ncbi) # 
    error_clades_lst = subclades_filter(error_clades_lst)
    for error_clade in error_clades_lst:
        phyla_cnt_ranks_dic = {phyla:{'count':None,'ranks':list()} for phyla in error_clade.values()}
        for rnk,phyla in error_clade.items():
            phyla_cnt_ranks_dic.get(phyla).get('ranks').append(rnk)
        for phyla in phyla_cnt_ranks_dic.keys():
            phyla_cnt_ranks_dic.get(phyla)['count'] = len(phyla_cnt_ranks_dic.get(phyla).get('ranks'))
        rank_conuts = [cnt_rnklst.get('count') for cnt_rnklst in phyla_cnt_ranks_dic.values()]
        mn, mx = min(rank_conuts), max(rank_conuts)
        ranks = [cnt_rnklst.get('ranks') for cnt_rnklst in filter(lambda cnt_rnklst:cnt_rnklst.get('count') == mn or cnt_rnklst.get('count') / mx < 0.3, [cnt_rnklst for cnt_rnklst in phyla_cnt_ranks_dic.values()])]
        error_ranks.append(list(itertools.chain.from_iterable(ranks)))
    error_ranks = list(itertools.chain.from_iterable(error_ranks))
    error_ranks = [str(rnk) for rnk in error_ranks] # added after figuering that the output is actualy integers.
    with open(os.path.join(directory,f'review_{rank}_to_offset.pkl'), 'wb') as f:
        pickle.dump(error_ranks,f)

    return error_ranks

########################################################

def metadata(conc_algn, taxid_map, ncbi):
    header_phylum_dic = {aligned_seq.name.split(' ',1)[0]:aligned_seq.description.split('OX=')[-1] for aligned_seq in conc_algn} # e.g. -> {'Oscillochloridaceae':'765420', 'Gallionellaceae': 1188319',...}
    header_phylum_dic = {header_name:taxid_map.get(rank_taxid).get('phylum') for header_name,rank_taxid in header_phylum_dic.items()} # the taxids here are of the phylum ranks.
    header_phylum_dic = {header_name:taxid_to_name(rank_phylum, ncbi) for header_name,rank_phylum in header_phylum_dic.items()} # values here altered from phylums taxids to their actual names.
    dic_for_df = {'rnk_txn':list(),'rnk_txd':list(),'phylum':list()}
    for hdr,phylum in header_phylum_dic.items():
        dic_for_df.get('rnk_txn').append(hdr)
        dic_for_df.get('rnk_txd').append(name_to_taxid(hdr.replace('_',' '), ncbi))
        dic_for_df.get('phylum').append(phylum)
    mtdta_df = pd.DataFrame.from_dict(dic_for_df)
    
    no_of_colors = len(set(header_phylum_dic.values())) # number of colors as the number of phylums present
    colors = ["#"+''.join([random.choice('0123456789ABCDEF') for i in range(6)]) for j in range(no_of_colors)] # list of colors as the nunber of phylums...
    colors = dict(zip(set(header_phylum_dic.values()),colors)) # mapping all the phylums into color
    mtdta_df['color'] = mtdta_df.loc[:,'phylum'].apply(lambda phyl:colors.get(phyl))
    
    # adds another column with phylum/class annotation - 'phyla_class'
    # while all phylums are being annotated as phylum within the new column, 
    # Proteobacteria only, are annotated as class. that is because the proteobcateria is so large clade, to the point 
    # of it being uninformative enough, so the lower rank (class) was chosen regarding members of it.

    mtdta_df['phyla_class'] = mtdta_df.apply(lambda row: name_to_taxid(row['rnk_txn'].replace('_',' '), ncbi)
                                                         if row['phylum'] == 'Proteobacteria'
                                                         else None,
                                                         axis=1) # phyla_class columns is now None for every non-proteobacteria leaf, and family/order taxid for every protebacteria leaf.
        
    mtdta_df['phyla_class'] = mtdta_df['phyla_class'].apply(lambda x: taxid_map.get(x) if x else None) # non-proteobacteria leaf -> None. proteobcateria -> the taxid_map header (dict) in respect to the corresponding family/order.
    mtdta_df['phyla_class'] = mtdta_df['phyla_class'].apply(lambda x: x.get('class') if x else None) # non-proteobacteria leaf -> None. proteobcateria -> taxid for the class rank of every leaf.
    mtdta_df['phyla_class'] = mtdta_df['phyla_class'].apply(lambda x: taxid_to_name(x, ncbi) if x else None) # non-proteobacteria leaf -> None. proteobcateria -> name (taxid->name) of the class rank of every leaf.
    mtdta_df['phyla_class'] = mtdta_df.apply(lambda x: x['phylum'] if not x['phyla_class'] else x['phyla_class'], axis=1) # non-proteobacteria leaf -> name of the phylum rank of every leaf. None. proteobcateria -> no change.

    return mtdta_df

########################################################

def naming_inner_nodes(tree_file, mtdta_df):
    # with the help of ete3 get_common_ancsestor method, inner nodes who represent phylum/class monophylentic groups are annotated within the
    # tree object. later, the tree object will be exported to ggtree as newick foramt, and ggtree will use the inner
    # nodes markers to annotated the figure. when the inner node is the root, which means that some phylogenetic problem is found within
    # the tree, nothing is marked, but the missing_inner_nodes list is updated to be examined later on.
    
    tree = Tree(tree_file)
    missing_inner_nodes = list()
    for phyla_name, sub_df in mtdta_df.groupby('phyla_class'):
        node = tree.get_common_ancestor(list(sub_df.loc[:,'rnk_txn']))    
        if node.is_root():
            missing_inner_nodes.append(phyla_name)
            continue
        node.name = phyla_name
    return tree, missing_inner_nodes

########################################################

def main(hmm = phylogenetics.ribosomal_HMMs('HMMs'),
         hmms = phylogenetics.ribosomal_HMMs('hmm_search_tsv'),
         target = databases.ref_seq_file_paths('bacteria_fasta'),
         rank = 'family',
         filter_candidatus = True,
         filter_leaves = False):

    taxid_map = phylogenetics.taxid_mapping()
    ncbi = NCBITaxa()

    directory = os.path.join(r'/davidb/yatirsolan/review_tree/',rank)
    if not os.path.isdir(directory):
        os.mkdir(directory)

    if not os.path.isfile(hmms): # hmms is first being taken from the main function default 'hmms' value
        hmms = HMM.hmmsearch(hmm,target,evalue=1e-12)
        hmms = HMM.hmmsearch_to_tsv(hmms)

    hmmsdf = pd.read_table(hmms)
    hmmsdf = hmmsearch_filter(hmmsdf)
    hmmsdf = hmmsdf_rank_annotation(hmmsdf, rank, taxid_map)

    if filter_candidatus:
        hmmsdf = candidatus_filter(hmmsdf, rank, ncbi, taxid_map)

    pickle_file = databases.ref_seq_file_paths('bacteria_SeqIO_dic_pickle') # pickle containing all RefSeq data 
    if os.path.isfile(pickle_file):
        target = pickle.load(open(pickle_file,'rb')) # if the 'if' statement is False the target is being brought from the function 'target' default value.

    if filter_leaves:
        tree_file,cncatntd_MSA = tree_construction_pipeline(hmmsdf, rank, target, directory, taxid_map, ncbi, filter_leaves)
        offset_leaves = leaf_verifier(tree_file, rank, taxid_map, ncbi, directory)
        if offset_leaves: # empty list is considered False in python...
            tree_file,cncatntd_MSA = tree_construction_pipeline(hmmsdf, rank, target, directory, taxid_map, ncbi, filter_leaves, offset_leaves)
    else:
        tree_file,cncatntd_MSA = tree_construction_pipeline(hmmsdf, rank, target, directory, taxid_map, ncbi, filter_leaves)

    mtdta_df = metadata(cncatntd_MSA, taxid_map, ncbi)
    tree, missing_inner_nodes = naming_inner_nodes(tree_file, mtdta_df)
    tree.write(outfile=tree_file, format=1) # overides the original newick file. (format 1 is nesceserry for the label name to be inside the newick output).
    mtdta_df.to_csv(os.path.join(directory, f"review_{rank}_phylogenetic_mtdta{'_filtered' if filter_leaves else str()}.tsv"), index=None, sep='\t')
    with open(os.path.join(directory, f"review_{rank}_missing_inner_nodes{'_filtered' if filter_leaves else str()}.pkl"), 'wb') as f:
        pickle.dump(missing_inner_nodes, f)

### main ### main ### main ### main ### main ### main ###

if __name__ == '__main__':
    main(rank='family')
