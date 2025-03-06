import os
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import DataStructs
#from clip import helpers
from multiprocessing import Pool
from rdkit.Chem import rdFingerprintGenerator 



def morgan_from_smiles(smiles, radius=3, nbits=1024, chiral=True):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return print(f"Could not generate mol object from SMILES: {smiles}")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=nbits, includeChirality=True)
    fp = generator.GetFingerprint(mol)


    #fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=nbits, useChirality=chiral)
    arr = np.zeros((0,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp,arr)
    return arr


if __name__ == '__main__':
    indir = "/home/vera/vsc5/data/AI/masterthesis/"
    index = "source_1_named_extract.csv"
    index = os.path.join(indir, index)

    outdir = "/home/vera/vsc5/data/AI/masterthesis/"
    outfile_hdf = "morgan_chiral_fps.hdf5"
    outfile_hdf = os.path.join(outdir, outfile_hdf)

    #n_cpus = 60

    csv = pd.read_csv(index)

    #csv["ID"] = csv.apply(lambda row: "-".join([str(row["PLATE_ID"]), str(row["WELL_POSITION"]),  str(row["SITE"])]), axis=1)

    ids = csv["Metadata_Sample_ID"]
    smiles = csv["SMILES"]


    fps = map(morgan_from_smiles, smiles)
    fps = list(fps)


    #fps = helpers.parallelize(morgan_from_smiles, smiles, n_cpus)

    columns = [str(i) for i in range(fps[0].shape[0])]

    df = pd.DataFrame(fps, index=ids, columns=columns)
    df.to_hdf(outfile_hdf, key="df", mode="w")

