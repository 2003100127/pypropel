__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np
import pandas as pd
from pypropel.util.Reader import Reader as pfreader
from pypropel.util.Writer import Writer as pfwriter
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.ComputLib import ComputLib as computlib
from pypropel.util.Console import Console


class Reader:

    def __init__(
            self,
            verbose: bool = True,
    ):
        self.__sort_ = -1
        self.pfreader = pfreader()
        self.pfwriter = pfwriter()
        self.computlib = computlib()
        self.console = Console()
        self.console.verbose = verbose

    @property
    def sort_(self):
        return self.__sort_

    @sort_.setter
    def sort_(
            self, value):
        print('Please note that you are attempting externally.')
        if value > 7 or value < 0:
            raise ValueError(
                '`sort_` has yet to reach there.',
                '| 1: return results for entire-chain residue contacts.',
                '| 2: return results of residue contacts by given pairs of interest.',
                '| 3: return sorted results by `score`',
                '| 4: return sorted results by `contact_id_1` and `contact_id_2`',
                '| 5: return dict results of a predictor',
                '| 6: return results of a residue of a predictor',
                '| 7: return cumulative dict results of a predictor',
                '| else: return raw results of a predictor',
                '| beyond: you need to choose one of opts above.',
            )
        else:
            self.__sort_ = value

    def sort_1(
            self, recombine, dp, df_dist):
        """
        ..  @description
            ------------

        :param recombine:
        :param dp: dist ob
        :param df_dist: a df of distances of residue pairs
        :return:
        """
        recombine_ = recombine
        # print(recombine_)
        recombine_dist = dp.segment(df_dist, recombine_)
        return recombine_, recombine_dist

    def sort_2(
            self, recombine, by='score'):
        """
        :param recombine:
        :param by: 'interact_id', 'score'
        :return:
        """
        recombine_ = recombine
        recombine_ = recombine_.sort_values(
            by=[by],
            ascending=True
        ).reset_index(drop=True)
        return recombine_

    def sort_3(
            self, recombine, is_uniform=False, uniform_df=None):
        """
        ..  @description:
            -------------
            select data by specifying seq_sep_inferior and seq_sep_superior.
            The select data can be sorted by two ways:
            1.  'score'
            2.  'contact_id_1' and 'contact_id_2'

        :param recombine: results of a predictor
        :param is_sort: False
        :return:
        """
        recombine_ = self.sort_2(recombine, by='interact_id')
        # # /*** block 1 ***/
        if is_uniform:
            # print(recombine_)
            predict_dict = self.computlib.tactic6(recombine_.values.tolist())
            # print(predict_dict)
            uniform_df[1] = 0
            # print(uniform_df)
            uniform = uniform_df.values.tolist()
            uniform_shape = len(uniform)
            for i in range(uniform_shape):
                id = uniform[i][0]
                try:
                    uniform[i][1] = predict_dict[id]
                except KeyError:
                    continue
            recombine_ = pd.DataFrame(uniform)
            # print(recombine_)
            recombine_.columns = [
                'interact_id',
                'score'
            ]
        return recombine_

    def mbpredcombined(
            self,
            mbp_cyto_path,
            mbp_tmh_path,
            mbp_extra_path,
            file_name,
            file_chain,
            phobius_path,
            xml_path,
            df_dist=None,
    ):
        seq_chain = chainname().seqchain(file_chain)
        mbp_cyto = self.mbpred(
            mbp_path=mbp_cyto_path,
            file_name=file_name,
            file_chain=file_chain,
            sort_=2,
            by='interact_id'
        )
        mbp_tmh = self.mbpred(
            mbp_path=mbp_tmh_path,
            file_name=file_name,
            file_chain=file_chain,
            sort_=2,
            by='interact_id'
        )
        mbp_extra = self.mbpred(
            mbp_path=mbp_extra_path,
            file_name=file_name,
            file_chain=file_chain,
            sort_=2,
            by='interact_id'
        )
        pdbtm_seg, phobius_seg = pdbnontmh().judge(
            df_dist=df_dist,
            pred_path=phobius_path,
            prot_name=file_name,
            xml_path=xml_path,
            seq_chain=seq_chain,
            kind='interaction_site'
        )
        # print('sjf')
        # print(pdbtm_seg)
        # print(sscenario().toSingle(phobius_seg['cyto_lower'], phobius_seg['cyto_upper']))
        pdbtm_pos_cyto = np.squeeze(sscenario().toSingle(pdbtm_seg['cyto_lower'], pdbtm_seg['cyto_upper']), axis=1) if pdbtm_seg['cyto_upper'] != [] else []
        pdbtm_pos_tmh = np.squeeze(sscenario().toSingle(pdbtm_seg['tmh_lower'], pdbtm_seg['tmh_upper']), axis=1) if pdbtm_seg['tmh_upper'] != [] else []
        pdbtm_pos_extra = np.squeeze(sscenario().toSingle(pdbtm_seg['extra_lower'], pdbtm_seg['extra_upper']), axis=1) if pdbtm_seg['extra_upper'] != [] else []
        phobius_pos_cyto = np.squeeze(sscenario().toSingle(phobius_seg['cyto_lower'], phobius_seg['cyto_upper']), axis=1) if phobius_seg['cyto_upper'] != [] else []
        phobius_pos_tmh = np.squeeze(sscenario().toSingle(phobius_seg['tmh_lower'], phobius_seg['tmh_upper']), axis=1) if phobius_seg['tmh_upper'] != [] else []
        phobius_pos_extra = np.squeeze(sscenario().toSingle(phobius_seg['extra_lower'], phobius_seg['extra_upper']), axis=1) if phobius_seg['extra_upper'] != [] else []

        mbp_pdbtm_cyto = mbp_cyto[mbp_cyto['interact_id'].isin(pdbtm_pos_cyto)]
        mbp_pdbtm_tmh = mbp_tmh[mbp_tmh['interact_id'].isin(pdbtm_pos_tmh)]
        mbp_pdbtm_extra = mbp_extra[mbp_extra['interact_id'].isin(pdbtm_pos_extra)]
        mbp_phobius_cyto = mbp_cyto[mbp_cyto['interact_id'].isin(phobius_pos_cyto)]
        mbp_phobius_tmh = mbp_tmh[mbp_tmh['interact_id'].isin(phobius_pos_tmh)]
        mbp_phobius_extra = mbp_extra[mbp_extra['interact_id'].isin(phobius_pos_extra)]
        # print(phobius_pos_cyto)
        # print(mbp_phobius_cyto)
        # print(mbp_phobius_tmh)
        # print(mbp_phobius_extra)
        pdbtm = pd.concat([mbp_pdbtm_cyto, mbp_pdbtm_tmh, mbp_pdbtm_extra]).reset_index(drop=True)
        phobius = pd.concat([mbp_phobius_cyto, mbp_phobius_tmh, mbp_phobius_extra]).reset_index(drop=True)
        placeholders_pdbtm = pd.DataFrame(np.zeros(len(pdbtm)))
        placeholders_phobius = pd.DataFrame(np.zeros(len(phobius)))
        # print(placeholders_phobius)
        pdbtm = pd.concat([placeholders_pdbtm, placeholders_pdbtm, pdbtm], axis=1).reset_index(drop=True)
        phobius = pd.concat([placeholders_phobius, placeholders_phobius, phobius], axis=1).reset_index(drop=True)
        pdbtm.columns = [
            '',
            'Res',
            'AA',
            'Score'
        ]
        phobius.columns = [
            '',
            'Res',
            'AA',
            'Score'
        ]
        return pdbtm, phobius

    def generate(
            self,
            list_fpn,
            mbp_cyto_path,
            mbp_tmh_path,
            mbp_extra_path,
            dist_path,
            phobius_path,
            xml_path,
            sv_pdbtm_fpn=None,
            sv_phobius_fpn=None,
    ):
        prot_df = self.pfreader.generic(list_fpn)
        for i, prot_name in enumerate(prot_df[0]):
            prot_chain = prot_df[1][i]
            file_chain = chainname().chain(prot_chain)
            dp = dlable(
                dist_path=dist_path,
                prot_name=prot_name,
                file_chain=file_chain,
            )
            df_dist = dp.attach()
            pdbtm, phobius = self.mbpredcombined(
                mbp_cyto_path=mbp_cyto_path,
                mbp_tmh_path=mbp_tmh_path,
                mbp_extra_path=mbp_extra_path,
                file_name=prot_name,
                file_chain=file_chain,
                phobius_path=phobius_path,
                xml_path=xml_path,
                df_dist=df_dist,
            )
            # pdbtm.columns = []
            self.pfwriter.generic(pdbtm, sv_fpn=sv_pdbtm_fpn + prot_name + file_chain + '.mbpred', df_sep=',', header=True)
            self.pfwriter.generic(phobius, sv_fpn=sv_phobius_fpn + prot_name + file_chain + '.mbpred', df_sep=',', header=True)
        return 0

    def mbpred(
            self, mbp_path, file_name, file_chain, dp=None, by=None, df_dist=None, pos_single_list=None, sort_=0, is_sort=False, id=0, L=50, len_seq=50):
        self.__sort_ = sort_
        results = self.pfreader.generic(
            mbp_path + file_name + file_chain + '.mbpred',
            df_sep=',',
            header=0,
            is_utf8=True
        )
        results.columns = [
            'index',
            'aa',
            'interact_id',
            'score'
        ]
        # print(results)
        results['aa'] = results['aa'].astype(str)
        recombine = results[[
            'interact_id',
            'score'
        ]]
        # print(recombine)
        # print(recombine.dtypes)
        if self.__sort_ == 1:
            pair_df = pd.DataFrame(pos_single_list)
            # print(pair_df)
            recombine = self.sort_3(
                recombine,
                is_uniform=True,
                uniform_df=pair_df
            )
            # print(recombine)
            dist_mbpred, dist_true = self.sort_1(recombine, dp, df_dist)
            return dist_mbpred, dist_true
        elif self.__sort_ == 2:
            return self.sort_2(recombine, by=by)
        elif self.__sort_ == 3:
            pass
        elif self.__sort_ == 4:
            pass
        elif self.__sort_ == 5:
            pass
        else:
            return recombine

    def delphi(
            self,
            delphi_path,
            file_name,
            file_chain,
            dp=None,
            by=None,
            df_dist=None,
            pos_single_list=None,
            sort_=0,
            is_sort=False,
            id=0,
            L=50,
            len_seq=50,
    ):
        self.__sort_ = sort_
        delphi_fpn = delphi_path + file_name + file_chain + '.txt'
        with open(delphi_fpn) as file:
            cues = []
            for line in file:
                if line.split()[0] == '#':
                    continue
                else:
                    cues.append(line.split())
            # print(cues)
            results = pd.DataFrame(cues)
            # print(results)
        results.columns = [
            'interact_id',
            'aa',
            'score'
        ]
        # print(results)
        results['aa'] = results['aa'].astype(str)
        results['interact_id'] = results['interact_id'].astype(int)
        results['score'] = results['score'].astype(float)
        recombine = results[[
            'interact_id',
            'score'
        ]]
        # print(recombine)
        # print(recombine.dtypes)
        if self.__sort_ == 1:
            pair_df = pd.DataFrame(pos_single_list)
            print(pair_df)
            recombine = self.sort_3(
                recombine,
                is_uniform=True,
                uniform_df=pair_df
            )
            # print(recombine)
            dist_delphi, dist_true = self.sort_1(recombine, dp, df_dist)
            return dist_delphi, dist_true
        elif self.__sort_ == 2:
            return self.sort_2(recombine, by=by)
        elif self.__sort_ == 3:
            pass
        elif self.__sort_ == 4:
            pass
        elif self.__sort_ == 5:
            pass
        else:
            return recombine

    def tma300(
            self,
            tma300_path,
            file_name,
            file_chain,
            dp=None,
            by=None,
            df_dist=None,
            pos_single_list=None,
            sort_=0,
            is_sort=False,
            id=0,
            L=50,
            len_seq=50,
    ):
        self.__sort_ = sort_
        results = self.pfreader.generic(
            tma300_path + file_name + file_chain + '.tma300',
            df_sep='\t',
            header=None,
        )
        # print(results)
        results.columns = [
            'interact_id',
            'aa',
            'score'
        ]
        results['aa'] = results['aa'].astype(str)
        recombine = results[[
            'interact_id',
            'score'
        ]]
        # print(recombine)
        # print(recombine.dtypes)
        if self.__sort_ == 1:
            pair_df = pd.DataFrame(pos_single_list)
            # print(pair_df)
            recombine = self.sort_3(
                recombine,
                is_uniform=True,
                uniform_df=pair_df
            )
            # print(recombine)
            dist_tma300, dist_true = self.sort_1(recombine, dp, df_dist)
            return dist_tma300, dist_true
        elif self.__sort_ == 2:
            return self.sort_2(recombine, by=by)
        elif self.__sort_ == 3:
            pass
        elif self.__sort_ == 4:
            pass
        elif self.__sort_ == 5:
            pass
        else:
            return recombine


if __name__ == "__main__":
    from pypropel.path import to
    from pypropel.prot.sequence.Fasta import Fasta as sfasta
    import tmkit as tmk
    from tmkit.seqnetrr.combo.Length import length as lscenario
    from tmkit.position.scenario.Segment import Segment as sscenario
    from pypropel.prot.structure.distance.isite.Label import Label as dlabel

    INIT = {
        'delphi_path': to('data/predictor/ppi/delphi/tm_alpha_n30/'),
        'mbp_cyto_path': to('data/predictor/ppi/mbpred/mbpredcyto/tm_alpha_n60/'),
        'mbp_tmh_path': to('data/predictor/ppi/mbpred/mbpredtm/tm_alpha_n60/'),
        'mbp_extra_path': to('data/predictor/ppi/mbpred/mbpredextra/tm_alpha_n60/'),
        # 'mbp_combined_path': to('data/predictor/ppi/mbpred/tm_alpha_n60/combined/phobius/'),
        'sv_pdbtm_fpn': to('data/predictor/ppi/mbpred/mbpredcombined/tm_alpha_n60/pdbtm/'),
        'sv_phobius_fpn': to('data/predictor/ppi/mbpred/mbpredcombined/tm_alpha_n60/phobius/'),
    }

    dp = dlabel(
        dist_path=to('data/pdb/complex/pdbtm/'),
        prot_name='1aij',
        file_chain='L',
    )

    df_dist = dp.attach()
    print(df_dist)

    sequence = sfasta().get(
        fasta_fpn=to("data/fasta/1aijL.fasta")
    )
    print(sequence)
    len_seq = len(sequence)
    length_pos_list = lscenario().tosgl(len_seq)
    # print(length_pos_list)

    pdbtm_seg, pred_seg = tmk.topo.cepdbtm(
        pdb_fp=to('data/pdb/pdbtm/'),
        prot_name='1aij',
        seq_chain='L',
        file_chain='L',
        topo_fp=to('data/phobius/1aijL.jphobius'),
        xml_fp=to('data/xml/'),
        fasta_fp=to('data/fasta/'),
    )
    print(pdbtm_seg)
    print(pred_seg)

    pdbtm_pos_cyto = sscenario().to_single(pdbtm_seg['cyto_lower'], pdbtm_seg['cyto_upper'])
    pdbtm_pos_tmh = sscenario().to_single(pdbtm_seg['tmh_lower'], pdbtm_seg['tmh_upper'])
    pdbtm_pos_extra = sscenario().to_single(pdbtm_seg['extra_lower'], pdbtm_seg['extra_upper'])
    phobius_pos_cyto = sscenario().to_single(pred_seg['cyto_lower'], pred_seg['cyto_upper'])
    phobius_pos_tmh = sscenario().to_single(pred_seg['tmh_lower'], pred_seg['tmh_upper'])
    phobius_pos_extra = sscenario().to_single(pred_seg['extra_lower'], pred_seg['extra_upper'])

    phobius_list = phobius_pos_cyto + phobius_pos_tmh + phobius_pos_extra
    print(phobius_list)
    print(pdbtm_pos_tmh)
    p = Reader()

    # print(p.delphi(
    #     delphi_path=to('data/isite/delphi/'),
    #     file_name='3jcu',
    #     file_chain='H',
    #     pos_single_list=pdbtm_pos_tmh,
    #     # dp=dp,
    #     # df_dist=df_dist,
    #     sort_=1
    # ))

    print(p.tma300(
        tma300_path=to('data/isite/deeptminter/'),
        file_name='1aij',
        file_chain='L',
        pos_single_list=pdbtm_pos_tmh,
        dp=dp,
        df_dist=df_dist,
        sort_=1
    ))

    # print(p.mbpred(
    #     mbp_path=INIT['mbp_combined_path'],
    #     file_name=INIT['prot_name'],
    #     file_chain=INIT['file_chain'],
    #     pos_single_list=phobius_list,
    #     dp=dp,
    #     df_dist=df_dist,
    #     sort_=1
    # ))

    # print(p.mbpredcombined(
    #     mbp_cyto_path=INIT['mbp_cyto_path'],
    #     mbp_tmh_path=INIT['mbp_tmh_path'],
    #     mbp_extra_path=INIT['mbp_extra_path'],
    #     file_name=INIT['prot_name'],
    #     file_chain=INIT['file_chain'],
    #     phobius_path=INIT['phobius_path'],
    #     xml_path=INIT['xml_path'],
    #     df_dist=df_dist,
    # ))

    # print(p.generate(
    #     list_fpn=INIT['list_fpn'],
    #     mbp_cyto_path=INIT['mbp_cyto_path'],
    #     mbp_tmh_path=INIT['mbp_tmh_path'],
    #     mbp_extra_path=INIT['mbp_extra_path'],
    #     phobius_path=INIT['phobius_path'],
    #     xml_path=INIT['xml_path'],
    #     dist_path=INIT['dist_path'],
    #     sv_pdbtm_fpn=INIT['sv_pdbtm_fpn'],
    #     sv_phobius_fpn=INIT['sv_phobius_fpn'],
    # ))