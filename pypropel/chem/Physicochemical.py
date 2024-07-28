__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"

import numpy as np
from rdkit.Chem import Crippen
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors


class Physicochemical:

    def __init__(self):
        descriptors_interest = [
            'MaxEStateIndex', 'MinEStateIndex', 'MaxAbsEStateIndex',
            'MinAbsEStateIndex', 'qed', 'MolWt', 'HeavyAtomMolWt',
            'ExactMolWt', 'NumValenceElectrons', 'NumRadicalElectrons',
            'FpDensityMorgan1', 'FpDensityMorgan2',
            'FpDensityMorgan3', 'BalabanJ', 'BertzCT', 'Chi0', 'Chi0n',
            'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n',
            'Chi3v', 'Chi4n', 'Chi4v', 'HallKierAlpha', 'Kappa1',
            'Kappa2', 'Kappa3', 'LabuteASA', 'PEOE_VSA1', 'PEOE_VSA10',
            'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14',
            'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6',
            'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'SMR_VSA1', 'SMR_VSA10',
            'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6',
            'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10',
            'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3',
            'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7',
            'SlogP_VSA8', 'SlogP_VSA9', 'TPSA', 'EState_VSA1', 'EState_VSA10',
            'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4',
            'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8',
            'EState_VSA9', 'VSA_EState1', 'VSA_EState10', 'VSA_EState2',
            'VSA_EState3', 'VSA_EState4', 'VSA_EState5', 'VSA_EState6',
            'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'FractionCSP3',
            'HeavyAtomCount', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles',
            'NumAliphaticHeterocycles', 'NumAliphaticRings',
            'NumAromaticCarbocycles', 'NumAromaticHeterocycles',
            'NumAromaticRings', 'NumHAcceptors', 'NumHDonors', 'NumHeteroatoms',
            'NumRotatableBonds', 'NumSaturatedCarbocycles',
            'NumSaturatedHeterocycles', 'NumSaturatedRings', 'RingCount',
            'MolLogP', 'MolMR', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert',
            'fr_ArN', 'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO',
            'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S', 'fr_HOCCN',
            'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O',
            'fr_Ndealkylation1', 'fr_Ndealkylation2', 'fr_Nhpyrrole',
            'fr_SH', 'fr_aldehyde', 'fr_alkyl_carbamate', 'fr_alkyl_halide',
            'fr_allylic_oxid', 'fr_amide', 'fr_amidine', 'fr_aniline',
            'fr_aryl_methyl', 'fr_azide', 'fr_azo', 'fr_barbitur',
            'fr_benzene', 'fr_benzodiazepine', 'fr_bicyclic', 'fr_diazo',
            'fr_dihydropyridine', 'fr_epoxide', 'fr_ester', 'fr_ether',
            'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine',
            'fr_hdrzone', 'fr_imidazole', 'fr_imide', 'fr_isocyan',
            'fr_isothiocyan', 'fr_ketone', 'fr_ketone_Topliss', 'fr_lactam',
            'fr_lactone', 'fr_methoxy', 'fr_morpholine', 'fr_nitrile',
            'fr_nitro', 'fr_nitro_arom', 'fr_nitro_arom_nonortho',
            'fr_nitroso', 'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation',
            'fr_phenol', 'fr_phenol_noOrthoHbond', 'fr_phos_acid',
            'fr_phos_ester', 'fr_piperdine', 'fr_piperzine', 'fr_priamide',
            'fr_prisulfonamd', 'fr_pyridine', 'fr_quatN', 'fr_sulfide',
            'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole',
            'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_unbrch_alkane',
            'fr_urea'
        ]
        self.allowed_descriptors = set(descriptors_interest)

    def real(self, mol):
        return mol.GetPropsAsDict()

    def dgbid(self, mol):
        return mol.GetProp('DATABASE_ID')

    def descriptors1(self, mol, is_dict=False):
        """
        ..  @description:
            -------------
            It is essentially same with the way of descriptors2().
            The difference is from the module of "MoleculeDescriptors".

        ..  @example:
            ---------
            1. mol_Bivalirudin = Chem.MolFromSmiles('CC[C@H](C)[C@H](NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](CC(O)=O)NC(=O)CNC(=O)[C@H](CC(N)=O)NC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)[C@@H]1CCCN1C(=O)[C@H](CCCNC(N)=N)NC(=O)[C@@H]1CCCN1C(=O)[C@H](N)CC1=CC=CC=C1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC1=CC=C(O)C=C1)C(=O)N[C@@H](CC(C)C)C(O)=O')
            or
            2. mol_Bivalirudin = mols[0]

        :param mols:
        :return: a list (200)
        """
        names = [d[0] for d in Descriptors._descList]
        dpt = MoleculeDescriptors.MolecularDescriptorCalculator(names)
        descriptor = dpt.CalcDescriptors(mol)
        print(names)
        if is_dict:
            descriptors = dict()
            for i, v in enumerate(names):
                descriptors[v] = descriptor[i]
        else:
            descriptors = np.array(descriptor)
        return descriptors

    def descriptors2(self, mol, is_dict=False):
        """
        ..  @description
            ------------
            1. Descriptors._descList is a list consisting of elements each
                of which is a 2d array.
            2. It is essentially same with the way of descriptors1().

        :param mol: should be single molecule.
        :return: a dict or an array
        """
        descriptors_dict = {}
        _descList = []
        descriptors_list = []
        # print(len(Descriptors._descList))
        # print(len(self.allowedDescriptors))
        for descriptor, function in Descriptors._descList:
            if descriptor in self.allowed_descriptors:
                if is_dict:
                    descriptors_dict[descriptor] = function(mol)
                else:
                    descriptors_list.append(function(mol))
                    _descList.append((descriptor, function))
        if is_dict:
            descriptors = descriptors_dict
        else:
            descriptors = np.array(descriptors_list)
        return descriptors

    def crippen(self, mol, is_dict=False):
        LogP = Crippen.MolLogP(mol)
        molar_refractivity = Crippen.MolMR(mol)
        descriptor = [LogP, molar_refractivity]
        if is_dict:
            descriptors = {'LogP': LogP, 'MR': molar_refractivity}
        else:
            descriptors = np.array(descriptor)
        return descriptors


if __name__ == "__main__":
    from pypropel.path import to
    from pypropel.chem.SDF import SDF

    p1 = Physicochemical()
    mp = SDF(sdf_fpn=to('data/chem/approved.sdf'))
    mp.kind = 'all_dict'
    mols = mp.extract()

    # print(p1.dgbid(mols['DB09517']))

    # print(p1.real(mols['DB09517']))


    # print(p1.descriptors1(mols['DB09517'], is_dict=True))
    # print(p1.descriptors2(mols['DB01169'], is_dict=True))

    # print(p1.crippen(mols['DB09517'], is_dict=True))