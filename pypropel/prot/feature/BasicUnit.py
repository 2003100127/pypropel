__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import re
import pandas as pd
from pypropel.prot.feature.alignment.Conservation import Conservation as conser


class BasicUnit:

    def __init__(self, ):
        self.masp1_unit = 'GGAGQGGQGGYGRGGYGQGGAGQGGAGAAAAAAAAAR'
        self.masp2_unit = 'GGAGPGRQQAYGPGGSGAAAAAAAAAR'
        self.masp1 = ''
        for i in range(12):
            if i != 11:
                self.masp1 = self.masp1 + self.masp1_unit
            else:
                self.masp1 = self.masp1 + 'GGAGQGGQGGYGRGGYGQGGAGQGGAGAAAAAAAAAS'
        # print(len(self.masp1))
        self.masp2 = ''
        for i in range(32):
            if i != 31:
                self.masp2 = self.masp2 + self.masp2_unit
            else:
                self.masp2 = self.masp2 + 'GGAGPGRQQAYGPGGSGAAAAAAAAAS'
        # print(len(self.masp2))

    def helix(self, sequence):
        start = []
        end = []
        for m in re.finditer(r"GG[A-Z]", sequence):
            # print('%d-%d: %s' % (m.start()+1, m.end(), m.group(0)))
            start.append(m.start() + 1)
            end.append(m.end())
        return start, end

    def sheet(self, sequence):
        start = []
        end = []
        for m in re.finditer(r"AAAAAAAA", sequence):
            # print('%d-%d: %s' % (m.start()+1, m.end(), m.group(0)))
            start.append(m.start() + 1)
            end.append(m.end())
        return start, end

    def turn(self, sequence):
        start = []
        end = []
        for m in re.finditer(r"GPG", sequence):
            # print('%d-%d: %s' % (m.start()+1, m.end(), m.group(0)))
            start.append(m.start() + 1)
            end.append(m.end())
        for m in re.finditer(r"QQ", sequence):
            # print('%d-%d: %s' % (m.start()+1, m.end(), m.group(0)))
            start.append(m.start() + 1)
            end.append(m.end())
        return start, end

    def localize(self, helix_start, helix_end, sheet_start, sheet_end, turn_start, turn_end, conser_fpn):
        masp1_jsd = conser().jsd(jsd_fpn=conser_fpn)
        # print(masp1_jsd)
        masp1_jsd['ss_cls'] = 'Other'
        for i in range(masp1_jsd.shape[0]):
            mark = masp1_jsd.at[i, 0] + 1
            for u, x in enumerate(helix_start):
                helix_interv = pd.Interval(x, helix_end[u], closed='both')
                if mark in helix_interv:
                    # print(mark, helix_interv)
                    masp1_jsd.at[i, 'ss_cls'] = 'H'
            for v, y in enumerate(sheet_start):
                sheet_interv = pd.Interval(y, sheet_end[v], closed='both')
                if mark in sheet_interv:
                    # print(mark, sheet_interv)
                    masp1_jsd.at[i, 'ss_cls'] = 'B'
            for w, z in enumerate(turn_start):
                turn_interv = pd.Interval(z, turn_end[w], closed='both')
                if mark in turn_interv:
                    # print(mark, turn_interv)
                    masp1_jsd.at[i, 'ss_cls'] = 'T'
        return masp1_jsd

    def localize2(self, helix_start, helix_end, sheet_start, sheet_end, turn_start, turn_end, consurf_fpn):
        masp1_consurf = conser().consurf(consurf_fpn=consurf_fpn)
        masp1_consurf['ss_cls'] = 'Other'
        for i in range(masp1_consurf.shape[0]):
            mark = masp1_consurf.at[i, 0] + 1
            for u, x in enumerate(helix_start):
                helix_interv = pd.Interval(x, helix_end[u], closed='both')
                if mark in helix_interv:
                    # print(mark, helix_interv)
                    masp1_consurf.at[i, 'ss_cls'] = 'H'
            for v, y in enumerate(sheet_start):
                sheet_interv = pd.Interval(y, sheet_end[v], closed='both')
                if mark in sheet_interv:
                    # print(mark, sheet_interv)
                    masp1_consurf.at[i, 'ss_cls'] = 'B'
            for w, z in enumerate(turn_start):
                turn_interv = pd.Interval(z, turn_end[w], closed='both')
                if mark in turn_interv:
                    # print(mark, turn_interv)
                    masp1_consurf.at[i, 'ss_cls'] = 'T'
        return masp1_consurf


if __name__ == "__main__":
    p = basicUnit()
    print(p.masp1)
    print(p.masp2)
    # print(p.helix(p.masp1))
    # print(p.sheet(p.masp1))
    # print(p.turn(p.masp2))

    conser_fpn = to('data/predictor/conservation/jsdivergence/masp/MaSp1.jsd')
    # conser_fpn = to('data/predictor/conservation/jsdivergence/masp/MaSp2.jsd')
    consurf_fpn = to('data/predictor/conservation/consurf/masp/MaSp1.consurf')

    helix_start = p.helix(p.masp1)[0]
    helix_end = p.helix(p.masp1)[1]
    sheet_start = p.sheet(p.masp1)[0]
    sheet_end = p.sheet(p.masp1)[1]
    turn_start = p.turn(p.masp1)[0]
    turn_end = p.turn(p.masp1)[1]
    # print(sheet_start)
    # print(sheet_end)
    # print(p.localize(helix_start, helix_end, sheet_start, sheet_end, turn_start, turn_end, conser_fpn))
    print(p.localize2(helix_start, helix_end, sheet_start, sheet_end, turn_start, turn_end, consurf_fpn))