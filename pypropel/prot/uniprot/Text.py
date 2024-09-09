__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import json
import numpy as np
import pandas as pd
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.ComputLib import ComputLib as computlib


class Text:

    def __init__(
            self,
            text_fpn: str,
            sv_json_fpn: str,
            sv_df_fpn: str,
    ):
        self.text_fpn = text_fpn
        self.sv_json_fpn = sv_json_fpn
        self.sv_df_fpn = sv_df_fpn
        self.pfwriter = pfwriter()

    def parse(self, ) -> pd.DataFrame:

        f = open(self.text_fpn)
        cands = {}
        # print(f.readlines())
        ids = []
        pac = []
        gns = []
        descriptions = []
        ensembl_G_ids = []
        ensembl_T_ids = []
        seqs = []
        tms_nums = []
        binding_nums = []
        binding_sites = []
        pdb_nums = []
        pdb_ids = []
        pdb_chains = []
        pdb_rezs = []
        count = 0
        pac_tmp = []
        gn_tmp = []
        de_tmp = []
        hpa_tmp = []
        tm_tmp = []
        binding_tmp = []
        pdb_tmp = []
        region_value_tmp = []
        region_tmp = []
        seq_tmp = []
        inExtracting = False
        inRecordingSeqMode = False
        for line in f.read().splitlines():
            # if line.startswith('//') != True:
                # print(1)
            acronym_line = line[:2]
            # print(acronym_line)
            if acronym_line != '//':
                if acronym_line == 'ID':
                    # print(line.split()[1])
                    id = line.split()[1]
                    ids.append(id)
                    count += 1
                if acronym_line == 'DE':
                    de_tmp.append(line.split())
                if acronym_line == 'AC':
                    pac_tmp.append(line.split())
                if acronym_line == 'GN':
                    gn_tmp.append(line.split()[1])
                    # print(gn_tmp)
                if acronym_line == 'DR':
                    ft_op = line.split()
                    # print(ft_op)
                    if ft_op[1] == 'Ensembl;':
                        hpa_tmp.append(ft_op)
                    if ft_op[1] == 'PDB;':
                        pdb_tmp.append(ft_op)
                        # print(pdb_tmp)
                if acronym_line == 'FT':
                    ft_op = line.split()
                    # print(ft_op)
                    if ft_op[1] == 'TRANSMEM':
                        tm_tmp.append(ft_op)
                        # print(tm_tmp)
                    if line.startswith('FT   TOPO_DOM'):
                        region_value_tmp.append(ft_op)
                        inExtracting = True
                        continue
                    if inExtracting == True:
                        region_tmp.append(ft_op)
                        inExtracting = False
                    if ft_op[1] == 'BINDING':
                        # print(ft_op)
                        binding_tmp.append(ft_op)
                        # print(binding_tmp)
                if not inRecordingSeqMode:
                    if line.startswith('SQ'):
                        inRecordingSeqMode = True
                        # print(line.split())
                elif line.startswith('ID'):
                    inRecordingSeqMode = False
                else:
                    seq_tmp.append(line.split())
            else:
                descriptions.append(' '.join([' '.join(part[1:]) for part in de_tmp]))
                seqs.append(''.join([''.join(part) for part in seq_tmp]))
                primary_ac = pac_tmp[0][1].split(';')[0]
                # print(primary_ac)
                pac.append(primary_ac)
                # print(gn_tmp)
                if gn_tmp != []:
                    gene_name = gn_tmp[0].split('=')[1].split(';')[0]
                else:
                    gene_name = 'NaN'
                # print(gene_name)
                gns.append(gene_name)
                # print(hpa_tmp)
                if hpa_tmp != []:
                    ensembl_id = hpa_tmp[0][4].split('.')[0]
                else:
                    ensembl_id = 'NaN'
                ensembl_G_ids.append(ensembl_id)
                ensembl_T_ids.append(';'.join([dr_i[2].split(';')[0] for dr_i in hpa_tmp]))
                # print(transcript_id)
                # print(ensembl_id)
                binding_sites.append(';'.join([bind_i[2].split(';')[0] for bind_i in binding_tmp]))
                pdb_ids.append(';'.join([pdb_i[2].split(';')[0] for pdb_i in pdb_tmp]))
                pdb_chains.append(';'.join([pdb_i[5].split('.')[0] if pdb_i[4] == '-;' else pdb_i[6].split('.')[0] for pdb_i in pdb_tmp]))
                pdb_rezs.append(';'.join(['0' if pdb_i[4] == '-;' else pdb_i[4].split(';')[0] for pdb_i in pdb_tmp]))

                if len(tm_tmp):
                    cands[ids[-1]] = {}
                    cands[ids[-1]]['pac'] = primary_ac
                    cands[ids[-1]]['gn'] = gene_name
                    cands[ids[-1]]['ensembl_id'] = ensembl_id
                    # print(len(tm_tmp))
                    tms_nums.append(len(tm_tmp))
                    binding_nums.append(len(binding_tmp))
                    pdb_nums.append(len(pdb_tmp))
                    cands[ids[-1]]['binding_num'] = len(binding_tmp)
                    cands[ids[-1]]['pdb_num'] = len(pdb_tmp)
                    x = np.array(tm_tmp)[:, 2].tolist()
                    x_ = [str(t).split('..') for t in x]
                    # tm_intervals.append(tm_tmp[:, 3])
                    tm_lower = []
                    tm_upper = []
                    cyto_lower = []
                    cyto_upper = []
                    extra_lower = []
                    extra_upper = []
                    for i in x_:
                        try:
                            tm_lower.append(int(i[0]))
                            tm_upper.append(int(i[1]))
                        except:
                            tm_lower.append(int(i[0].split(':')[1]))
                            tm_upper.append(int(i[1]))

                    # print(tm_lower)
                    # print(tm_upper)
                    # print(len(seqs[-1]))
                    nontm_lower, nontm_upper= computlib().reverseInterval(tm_lower, tm_upper, len(seqs[-1]))
                    # print(nontm_lower)
                    # print(nontm_upper)

                    # ### /*** block for region cyto and extra ***/
                    if len(region_value_tmp):
                        for rk, rv in enumerate(region_tmp):
                            if rv[1].split('"')[1] == 'Cytoplasmic':
                                cc = region_value_tmp[rk][2].split('..')
                                if len(cc) < 2:
                                    # print(cc)
                                    cyto_lower.append(int(cc[0]))
                                    cyto_upper.append(int(cc[0]))
                                else:
                                    if cc[0] == '?':
                                        continue
                                    elif cc[1] == '?':
                                        continue
                                    else:
                                        try:
                                            cyto_lower.append(int(cc[0]))
                                            cyto_upper.append(int(cc[1]))
                                        except:
                                            cyto_lower.append(int(cc[0].split(':')[1]))
                                            cyto_upper.append(int(cc[1]))

                            if rv[1].split('"')[1] == 'Extracellular':
                                cc = region_value_tmp[rk][2].split('..')
                                if len(cc) < 2:
                                    # print(cc)
                                    extra_lower.append(int(cc[0]))
                                    extra_upper.append(int(cc[0]))
                                else:
                                    if cc[0] == '?':
                                        continue
                                    elif cc[1] == '?':
                                        continue
                                    else:
                                        try:
                                            extra_lower.append(int(cc[0]))
                                            extra_upper.append(int(cc[1]))
                                        except:
                                            extra_lower.append(int(cc[0].split(':')[1]))
                                            extra_upper.append(int(cc[1]))
                        # print(cyto_lower)
                        # print(cyto_upper)
                        # print(extra_lower)
                        # print(extra_upper)
                    cands[ids[-1]]['cyto_lower'] = cyto_lower
                    cands[ids[-1]]['cyto_upper'] = cyto_upper
                    cands[ids[-1]]['extra_lower'] = extra_lower
                    cands[ids[-1]]['extra_upper'] = extra_upper
                    cands[ids[-1]]['tmh_lower'] = tm_lower
                    cands[ids[-1]]['tmh_upper'] = tm_upper
                    cands[ids[-1]]['nontmh_lower'] = nontm_lower
                    cands[ids[-1]]['nontmh_upper'] = nontm_upper
                    cands[ids[-1]]['seq'] = seqs[-1]
                else:
                    tms_nums.append(0)
                    binding_nums.append(len(binding_tmp))
                    pdb_nums.append(len(pdb_tmp))
                de_tmp = []
                hpa_tmp = []
                tm_tmp = []
                binding_tmp = []
                pdb_tmp = []
                seq_tmp = []
                pac_tmp = []
                gn_tmp = []
                region_value_tmp = []
                region_tmp = []
        # print(pac)
        # print(len(pac))
        # print(cands)
        # print(len(cands))
        with open(self.sv_json_fpn, 'w') as fp:
            json.dump(cands, fp)
        df = pd.DataFrame({
            'ID': ids,
            'AC': pac,
            'DE': descriptions,
            'GN': gns,
            'Ensembl_G_id': ensembl_G_ids,
            'Ensembl_T_id': ensembl_T_ids,
            'SQ': seqs,
            'tms_nums': tms_nums,
            'binding_nums': binding_nums,
            'binding_sites': binding_sites,
            'pdb_nums': pdb_nums,
            'pdb_ids': pdb_ids,
            'pdb_chains': pdb_chains,
            'pdb_rezs': pdb_rezs,
        })
        self.pfwriter.generic(df, sv_fpn=self.sv_df_fpn, header=True)
        return df


if __name__ == "__main__":
    from pypropel.path import to
    p = Text(
        # text_fpn=to('data/uniprot/text/uniprot-human-filtered-organism.txt'),
        text_fpn=to('data/uniprot/text/uniprotkb_Human_AND_reviewed_true_AND_m_2024_07_27.txt'),
        sv_json_fpn=to('data/uniprot/text/human.json'),
        sv_df_fpn=to('data/uniprot/text/human.txt'),
    )
    print(p.parse())
    df = p.parse()
    print(df.loc[df['tms_nums'] != 0])