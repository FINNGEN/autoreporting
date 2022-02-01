from .db import PhenoData, PhenoInfoDAO
import pandas as pd

class FGPhenoInfo(PhenoInfoDAO):
    def __init__(self, fname: str):
        #load data
        try:
            self.data = pd.read_csv(fname,sep="\t")
        except:
            print("Exception during phenotype info file loading")
            raise

    def get_pheno_info(self, phenotype: str) -> PhenoData:
        try:
            phenorow = self.data[self.data["phenocode"] == phenotype].iloc[0]
            data = PhenoData(
                phenotype,
                phenorow["name"],
                phenorow["category"],
                int(phenorow["num_cases"]),
                int(phenorow["num_controls"]),
            )
        #row not found
        except IndexError:
            print(f"PhenoInfo: no phenotype named {phenotype} found in phenotype data")
            data = PhenoData(phenotype,"NA","NA",0,0)
        except:
            print("Exception during phenotype info loading")
            raise
        return data

class UKBBPhenoInfo(PhenoInfoDAO):
    """Phenotype information DAO for Pan-UKBB data
    """
    def __init__(self, fname: str):
        try:
            data = pd.read_csv(fname, sep="\t")
            self.data = data[[
                "phenocode",
                "description",
                "category",
                "n_cases_EUR",
                "n_controls_EUR"
            ]]
        except:
            print("Exception during phenotype info file loading")
            raise

    def get_pheno_info(self, phenotype: str) -> PhenoData:
        try:
            phenorow = self.data[self.data["phenocode"] == phenotype].iloc[0]
            data = PhenoData(
                phenotype,
                phenorow["description"],
                phenorow["category"],
                int(phenorow["n_cases_EUR"]),
                int(phenorow["n_controls_EUR"]),
            )
        except:
            print("Exception during phenotype info loading")
            raise
        return data