__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import xml.etree.ElementTree as ET


class XML:

    def get(
            self,
            xml_path,
            xml_name,
            seq_chain,
    ):
        xml_fpn = xml_path + xml_name + '.xml'
        tree = ET.parse(xml_fpn)
        parser_pdb = tree.getroot()
        for chains in parser_pdb:
            if chains.tag == '{http://pdbtm.enzim.hu}CHAIN':
                for seqs in chains.iter('{http://pdbtm.enzim.hu}SEQ'):
                    if chains.get('CHAINID') == seq_chain:
                        # print(seqs.text)
                        fasta_seq = ''.join(seqs.text.split())
                        return fasta_seq