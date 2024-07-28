__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Avalon.pyAvalonTools import GetAvalonFP


class Fingerprint:

    def __init__(self, ):
        pass

    def get(
            self,
            mols,
            type='Morgan',
            is_to_vec=False,
            single=True,
    ):
        """

        Parameters
        ----------
        mols
        type
        self.dgb_mols_all or self.dgb_mols_approved
            1. MACCSKeys: 167bit,
            2. Morgan: changeable,
            3. Hashed_Torsion,
            4. Torsion,
            5. Atom_Pair,
            6. RDK,
            7. Avalon: 512 bit.
        is_to_vec
        single

        Returns
        -------

        """
        # print('===>calculating {} fingerprints'.format(type))
        if type == 'MACCSKeys':
            if single:
                fps = AllChem.GetMACCSKeysFingerprint(mols)
            else:
                fps = [[mol[0], AllChem.GetMACCSKeysFingerprint(mol[1])] for mol in mols]
        elif type == 'Morgan':
            if single:
                fps = AllChem.GetMorganFingerprintAsBitVect(mols, radius=2, nBits=1024)
            else:
                fps = [[mol[0], AllChem.GetMorganFingerprintAsBitVect(mol[1], radius=2, nBits=1024)] for mol in mols]
        elif type == 'Hashed_Torsion':
            if single:
                fps = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mols)
            else:
                fps = [[mol[0], rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol[1])] for mol in mols]
        elif type == 'Torsion':
            if single:
                fps = Torsions.GetTopologicalTorsionFingerprintAsIntVect(mols)
            else:
                fps = [[mol[0], Torsions.GetTopologicalTorsionFingerprintAsIntVect(mol[1])] for mol in mols]
        elif type == 'Atom_Pair':
            if single:
                fps = Pairs.GetAtomPairFingerprintAsBitVect(mols)
            else:
                fps = [[mol[0], Pairs.GetAtomPairFingerprintAsBitVect(mol[1])] for mol in mols]
        elif type == 'RDK':
            if single:
                fps = Chem.RDKFingerprint(mols, maxPath=5)
            else:
                fps = [[mol[0], Chem.RDKFingerprint(mol[1], maxPath=5)] for mol in mols]
        elif type == 'Avalon':
            if single:
                fps = GetAvalonFP(mols)
            else:
                fps = [[mol[0], GetAvalonFP(mol[1])] for mol in mols]
        else:
            if single:
                fps = Chem.RDKFingerprint(mols, maxPath=5)
            else:
                fps = [[mol[0], Chem.RDKFingerprint(mol[1], maxPath=5)] for mol in mols]
        if is_to_vec:
            if single:
                fps = fps.ToBitString()
            else:
                fps = [[fp[0], fp[1].ToBitString()] for i, fp in enumerate(fps)]
        return fps


if __name__ == "__main__":
    from pypropel.path import to
    from pypropel.chem.SDF import SDF

    p = Fingerprint()
    mp = SDF(sdf_fpn=to('data/chem/approved.sdf'))
    mp.kind = 'all_dict' # approved
    mols = mp.extract()
    # print(mols)
    fps = p.get(
        mols['DB17308'],
        is_to_vec=True,
        type='Morgan',
        single=True
    )
    print(fps)

    # print(len(fps))

    print(p.get(mols, is_to_vec=True, single=False))