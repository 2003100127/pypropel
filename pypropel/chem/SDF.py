__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"

from rdkit import Chem
from pypropel.util.Console import Console


class SDF:

    def __init__(
            self,
            sdf_fpn,
            kind='all',
            verbose=False,
    ):
        self.console = Console()
        self.verbose = verbose
        self.console.verbose = self.verbose
        self._kind = kind
        self.sdf_fpn = sdf_fpn
        self.all_mols = self.__all()
        self.console.print('===>Fetching mol info...')

    @property
    def kind(self):
        return self._kind

    @kind.setter
    def kind(self, value):
        self.console.print('======>Please note that this is from an external attempt')
        if value not in ['all', 'approved', 'expt', 'all_2d', 'approved_2d', 'expt_2d', 'all_dict', 'approved_dict',
                         'expt_dict']:
            raise ValueError(
                'kind has not yet reached there.',
                '| all: ',
                '| approved: ',
                '| expt: ',
                '| all_2d: ',
                '| approved_2d: ',
                '| expt_2d: ',
                '| all_dict: ',
                '| approved_dict: ',
                '| expt_dict: ',
            )
        else:
            self._kind = value

    def __all(self, ):
        return Chem.SDMolSupplier(self.sdf_fpn)

    def extract(self, ):
        if self._kind == 'all':
            self.console.print('=========>You are attempting to get all molecules (all).')
            return [mol for mol in self.all_mols if mol]
        elif self._kind == 'approved':
            self.console.print('=========>You are attempting to get approved molecules (approved).')
            return [mol for mol in self.all_mols if mol and 'approved' in mol.GetProp('DRUG_GROUPS')]
        elif self._kind == 'expt':
            self.console.print('=========>You are attempting to get approved molecules (expt).')
            return [mol for mol in self.all_mols if mol and 'experimental' in mol.GetProp('DRUG_GROUPS')]
        elif self._kind == 'all_2d':
            self.console.print('=========>You are attempting to get all name<->molecules (all_2d).')
            return [[str(mol.GetProp('DATABASE_ID')), mol] for mol in self.all_mols if mol]
        elif self._kind == 'approved_2d':
            self.console.print('=========>You are attempting to get approved name<->molecules (approved_2d).')
            return [[str(mol.GetProp('DATABASE_ID')), mol] for mol in self.all_mols if
                    mol and 'approved' in mol.GetProp('DRUG_GROUPS')]
        elif self._kind == 'expt_2d':
            self.console.print('=========>You are attempting to get experimental name<->molecules (expt_2d).')
            return [[str(mol.GetProp('DATABASE_ID')), mol] for mol in self.all_mols if
                    mol and 'experimental' in mol.GetProp('DRUG_GROUPS')]
        elif self._kind == 'all_dict':
            self.console.print('=========>You are attempting to get all molecules dictionary (all_dict).')
            mols = dict()
            for mol in self.all_mols:
                if mol:
                    # self.console.print(mol.GetProp('DRUG_GROUPS'))
                    mols[mol.GetProp('DATABASE_ID')] = mol
            return mols
        elif self._kind == 'approved_dict':
            self.console.print('=========>You are attempting to get approved molecules dictionary (approved_dict).')
            mols = dict()
            for mol in self.all_mols:
                if mol and 'approved' in mol.GetProp('DRUG_GROUPS'):
                    # self.console.print(mol.GetProp('DRUG_GROUPS'))
                    mols[mol.GetProp('DATABASE_ID')] = mol
            return mols
        elif self._kind == 'expt_dict':
            self.console.print('=========>You are attempting to get approved molecules dictionary (approved_dict).')
            mols = dict()
            for mol in self.all_mols:
                if mol and 'experimental' in mol.GetProp('DRUG_GROUPS'):
                    # print(mol.GetProp('DRUG_GROUPS'))
                    mols[mol.GetProp('DATABASE_ID')] = mol
            return mols


if __name__ == "__main__":
    from pypropel.path import to

    p = SDF(sdf_fpn=to('data/chem/approved.sdf'))
    # p.kind = 'expt_dict' # expt approved all
    p.kind = 'approved'
    print(p.kind)
    # print(p.extract())
    print(len(p.extract()))