__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import re
from protein.precision.ppi.Dispatch import dispatch
from protein.file.Control import control


class controller(object):

    def __init__(self, ):
        pass

    def scroll(self, ):
        p = dispatch()
        # meet = '/uniprot2016.02/aaanalysis/w3_9/'
        # meet = '/uniclust2020.01/aaanalysis/w3_9/filtered/reduced/n301/'
        # meet = '/uniclust2020.01/aaanalysis/w3_9/filtered/reduced/n301/i4/'
        meet = '/uniclust2020.01/aaanalysis/w3_9/filtered/reduced/n301_i0/'
        # meet = '/uniclust2020.01/aaanalysis/w3_9/filtered/reduced/cv/'
        # meet = '/uniclust2020.01/aaanalysis/w3_9/filtered/reduced/cv/cv4/'
        # category = 'tma203'
        category = 'tma301'
        tools = [
            'tma300',
            # 'delphi',
            # 'mbpred/mbpredcyto',
            # 'mbpred/mbpredtm',
            # 'mbpred/mbpredextra',
            # 'mbpred/mbpredall',
            # 'mbpred/mbpredcombined',
        ]
        models = [
            # random
            # '301-cv0-599',
            # '301-cv0-430',
            # '301-cv0-150',
            '301-cv0-300',

            # 'stacking/head',
            # 'stacking/tail',
            # 'stacking',

            # ### /* aaanalysis + filtered n301 */
            # '301-cv0-60',
            # '301-cv0-59',
            # '301-cv0-58',
            # '301-cv0-63',
            # '301-cv0-67',

            # ### /* aaanalysis + filtered cross validation */
            # '301-cv0-65',
            # '301-cv0-68',
            # '301-cv0-60',
            # '301-cv0-65',
            # '301-cv0-66',

            # ### /* n301_i0 */
            # '301-cv0-57',
            # '301-cv0-58',
            # '301-cv0-59',
            # '301-cv0-60',
            # '301-cv0-61',
            # '301-cv0-150',
            # '301-cv0-270',
            # '301-cv0-300',
            # '301-cv0-370',
            # '301-cv0-430',
            # '301-cv0-599',
            # '301-cv0-611',
            # '301-cv0-618',
            # '301-cv0-628',
            # '301-cv0-637',
            # 'ave',

            # ### /* Non-AAanalysis for tma300 */
            # '/300-cv0-401/',
            # '/300-cv0-532/',
            # '/300-cv0-536/',
            # '/300-cv0-538/',
            # '/300-cv0-539/',

            # AAanalysis for tma300
            # '/300-cv0-407/',
            # '/300-cv0-411/',
            # '/300-cv0-534/',
            # '/300-cv0-591/',
            # '/300-cv0-666/',

            # AAanalysis for tma203
            # '203-cv0-204',
            # '203-cv0-298',
            # '203-cv0-307',
            # '203-cv0-347',
            # '203-cv0-366',
        ]
        topologies = [
            'pdbtm',
            # 'phobius',
        ]
        regions = {
            'membrane': 'tmh',
            'cytoplasmic': 'cyto',
            'extracellular': 'extra',
            'combined': 'combined',
        }
        datasets = {
            # 'tm_alpha_n31': 'prot_n31.txt',
            'tm_alpha_n36': 'prot_n36.txt',
            # 'tm_alpha_n36': 'prot_delphi.txt',
            # 'tm_alpha_n36': 'prot_mbpredall.txt',
            # 'tm_alpha_n101': 'prot_n101.txt',
            # 'tm_alpha_n101': 'prot_delphi.txt',
            # 'tm_alpha_n101': 'prot_mbpredall.txt',
            # 'tm_alpha_n60': 'prot_n60.txt',
            # 'tm_alpha_n30': 'prot_mbpredall.txt',
            # 'tm_alpha_n30': 'prot_n30.txt',
            # 'tm_alpha_n30': 'prot_delphi.txt',
            # 'tm_alpha_n91': 'prot_n91.txt',
        }
        thresholds = {
            4: 'BordInter',
            # 5.5: 'FuchInter',
            # 6: 'RostInter',
        }
        for tool in tools:
            for dataset, prot_list in datasets.items():
                for topology in topologies:
                    for region, region_abbr in regions.items():
                        for thre, definition in thresholds.items():
                            if tool == 'tma300':
                                working_path = to('data/tool/ppisite/') + category + '/' + meet + '/' + dataset + '/' + models[0] + '/'
                                sv_fp = to('data/al/prediction/ppisite/') + category + '/' + meet + '/' + dataset + '/' + topology + '/' + region + '/' + definition + '/' + models[0] + '/'
                                control().create(
                                    DIRECTORY=sv_fp,
                                    mode='dir'
                                )
                            elif tool == 'delphi':
                                working_path = to('data/tool/ppisite/') + tool + '/' + dataset + '/'
                                sv_fp = to('data/al/prediction/ppisite/') + tool + '/' + dataset + '/' + topology + '/' + region + '/' + definition + '/'
                                control().create(
                                    DIRECTORY=sv_fp,
                                    mode='dir'
                                )
                            else:
                                working_path = to('data/tool/ppisite/') + tool + '/' + dataset + '/'
                                # working_path = to('data/tool/ppisite/') + tool + '/' + dataset + '/' + topology + '/'
                                sv_fp = to('data/al/prediction/ppisite/') + tool + '/' + dataset + '/' + topology + '/' + region + '/' + definition + '/'
                                control().create(
                                    DIRECTORY=sv_fp,
                                    mode='dir'
                                )
                            PARAMS = {
                                'sort': 1,
                                'mode': topology + '_' + region_abbr,
                                'list_fpn': to('data/protein/pdb/') + dataset + '/' + prot_list,
                                'fasta_path': to('data/protein/fasta/') + dataset + '/',
                                'xml_path': to('data/protein/xml/') + dataset + '/',
                                'dist_path': to('data/protein/dist/') + dataset + '/',
                                'phobius_path': to('data/tool/topology/phobius/') + dataset + '/',
                            }
                            print('{}: {} - {} - {} - {}'.format(tool, dataset, topology, region, definition))
                            print(p.segment(
                                list_fpn=PARAMS['list_fpn'],
                                xml_path=PARAMS['xml_path'],
                                mode=PARAMS['mode'],
                                tool_file_path=working_path,
                                tool=re.sub(r'/.*$', '', tool),
                                phobius_path=PARAMS['phobius_path'],
                                dist_path=PARAMS['dist_path'],
                                dist_limit=thre,
                                sort=PARAMS['sort'],
                                sv_fp=sv_fp,
                            ))
        return 0

    def rsa(self, ):
        p = dispatch()
        # meet = '/uniclust2020.01/aaanalysis/w3_9/filtered/reduced/n301/'
        # meet = '/uniclust2020.01/aaanalysis/w3_9/filtered/reduced/n301/i4/'
        meet = '/uniclust2020.01/aaanalysis/w3_9/filtered/reduced/n301_i0/'
        # meet = '/uniclust2020.01/aaanalysis/w3_9/filtered/reduced/cv/'
        # meet = '/uniclust2020.01/aaanalysis/w3_9/filtered/reduced/cv/cv4/'
        # category = 'tma203'
        category = 'tma301'
        tools = [
            'tma300',
            # 'delphi',
            # 'mbpred/mbpredcyto',
            # 'mbpred/mbpredtm',
            # 'mbpred/mbpredextra',
            # 'mbpred/mbpredall',
            # 'mbpred/mbpredcombined',
        ]
        models = [
            # random
            # '301-cv0-599',
            # '301-cv0-430',
            '301-cv0-300',
            # '301-cv0-150',

            # 'stacking/head',
            # 'stacking/tail',
            # 'stacking',

            # ### /* aaanalysis + filtered n301 */
            # '301-cv0-60',
            # '301-cv0-59',
            # '301-cv0-58',
            # '301-cv0-63',
            # '301-cv0-67',

            # ### /* aaanalysis + filtered cross validation */
            # '301-cv0-65',
            # '301-cv0-68',
            # '301-cv0-60',
            # '301-cv0-65',
            # '301-cv0-66',

            # ### /* n301_i0 */
            # '301-cv0-57',
            # '301-cv0-58',
            # '301-cv0-59',
            # '301-cv0-60',
            # '301-cv0-61',
            # '301-cv0-150',
            # '301-cv0-270',
            # '301-cv0-300',
            # '301-cv0-370',
            # '301-cv0-430',
            # '301-cv0-599',
            # '301-cv0-611',
            # '301-cv0-618',
            # '301-cv0-628',
            # '301-cv0-637',
            # 'ave',

            # ### /* Non-AAanalysis for tma300 */
            # '/300-cv0-401/',
            # '/300-cv0-532/',
            # '/300-cv0-536/',
            # '/300-cv0-538/',
            # '/300-cv0-539/',

            # AAanalysis for tma300
            # '/300-cv0-407/',
            # '/300-cv0-411/',
            # '/300-cv0-534/',
            # '/300-cv0-591/',
            # '/300-cv0-666/',

            # AAanalysis for tma203
            # '203-cv0-204',
            # '203-cv0-298',
            # '203-cv0-307',
            # '203-cv0-347',
            # '203-cv0-366',
        ]
        topologies = [
            'pdbtm',
            # 'phobius',
        ]
        regions = {
            # 'membrane': 'tmh',
            # 'cytoplasmic': 'cyto',
            # 'extracellular': 'extra',
            'combined': 'combined',
        }
        datasets = {
            # 'tm_alpha_n31': 'prot_n31.txt',
            'tm_alpha_n36': 'prot_n36.txt',
            # 'tm_alpha_n36': 'prot_delphi.txt',
            # 'tm_alpha_n36': 'prot_mbpredall.txt',
            # 'tm_alpha_n101': 'prot_n101.txt',
            # 'tm_alpha_n101': 'prot_delphi.txt',
            # 'tm_alpha_n101': 'prot_mbpredall.txt',
            # 'tm_alpha_n60': 'prot_n60.txt',
            # 'tm_alpha_n30': 'prot_mbpredall.txt',
            # 'tm_alpha_n30': 'prot_n30.txt',
            # 'tm_alpha_n30': 'prot_delphi.txt',
            # 'tm_alpha_n91': 'prot_n91.txt',
        }
        thresholds = {
            # 4: 'BordInter',
            # 5.5: 'FuchInter',
            6: 'RostInter',
        }
        for tool in tools:
            for dataset, prot_list in datasets.items():
                for topology in topologies:
                    for region, region_abbr in regions.items():
                        for thre, definition in thresholds.items():
                            if tool == 'tma300':
                                working_path = to('data/tool/ppisite/') + category + '/' + meet + '/' + dataset + '/' + models[0] + '/'
                                sv_fp = to('data/al/prediction/ppisite/') + category + '/' + meet + '/' + dataset + '/' + topology + '/' + region + '/' + definition + '/' + models[0] + '/rsa/'
                                control().create(
                                    DIRECTORY=sv_fp,
                                    mode='dir'
                                )
                            elif tool == 'delphi':
                                working_path = to('data/tool/ppisite/') + tool + '/' + dataset + '/'
                                sv_fp = to('data/al/prediction/ppisite/') + tool + '/' + dataset + '/' + topology + '/' + region + '/' + definition + '/rsa/'
                                control().create(
                                    DIRECTORY=sv_fp,
                                    mode='dir'
                                )
                            else:
                                working_path = to('data/tool/ppisite/') + tool + '/' + dataset + '/'
                                # working_path = to('data/tool/ppisite/') + tool + '/' + dataset + '/' + topology + '/'
                                sv_fp = to('data/al/prediction/ppisite/') + tool + '/' + dataset + '/' + topology + '/' + region + '/' + definition + '/rsa/'
                                control().create(
                                    DIRECTORY=sv_fp,
                                    mode='dir'
                                )
                            PARAMS = {
                                'sort': 1,
                                'rsa_path': to('data/tool/dssp/rsa/') + dataset + '/',
                                'rsa_thres': [0.1, 0.15, 0.2, 0.25],
                                'mode': topology + '_' + region_abbr,
                                'list_fpn': to('data/protein/pdb/') + dataset + '/' + prot_list,
                                'fasta_path': to('data/protein/fasta/') + dataset + '/',
                                'xml_path': to('data/protein/xml/') + dataset + '/',
                                'dist_path': to('data/protein/dist/') + dataset + '/',
                                'phobius_path': to('data/tool/topology/phobius/') + dataset + '/',
                            }
                            print('{}: {} - {} - {} - {}'.format(tool, dataset, topology, region, definition))
                            print(p.rsa(
                                list_fpn=PARAMS['list_fpn'],
                                xml_path=PARAMS['xml_path'],
                                mode=PARAMS['mode'],
                                tool_file_path=working_path,
                                tool=re.sub(r'/.*$', '', tool),
                                phobius_path=PARAMS['phobius_path'],
                                dist_path=PARAMS['dist_path'],
                                dist_limit=thre,
                                sort=PARAMS['sort'],
                                rsa_path=PARAMS['rsa_path'],
                                rsa_thres=PARAMS['rsa_thres'],
                                sv_fp=sv_fp,
                            ))
        return 0


if __name__ == "__main__":
    p = controller()
    print(p.scroll())
    # print(p.rsa())