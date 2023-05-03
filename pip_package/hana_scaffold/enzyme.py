# -*- coding: utf-8 -*-
from typing import List


class _EnzymeDict(dict):
    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            for key_item in key:
                super().__setitem__(key_item.upper(), value)
        else:
            super().__setitem__(key.upper(), value)


_HANA_ENZYME_MAP = _EnzymeDict()
_HANA_ENZYME_MAP["Arima"] = ["GATC", "GAAT", "GACT", "GAGT", "GATT"]
_HANA_ENZYME_MAP["TaqI", "Taq1"] = "TCGA"
_HANA_ENZYME_MAP["HaeIII", "Hae3"] = "GGCC"
_HANA_ENZYME_MAP["DpnI", "Dpn1", "DpnII", "Dpn2", "MboI", "Mbo1", "Sau3AI"] = "GATC"
_HANA_ENZYME_MAP["AluI", "Alu1"] = "AGCT"
_HANA_ENZYME_MAP["NlaIII", "Nla3", "FaeI", "Fae1", "FatI", "Fat1", "Hin1II", "Hsp92II"] = "CATG"
_HANA_ENZYME_MAP["HpaII", "Hpa2"] = "CCGG"
_HANA_ENZYME_MAP["FokI", "Fok1"] = "GGATG"
_HANA_ENZYME_MAP["AaaI", "Aaa1"] = "CGGCG"
_HANA_ENZYME_MAP["HgaI", "Hga1"] = "GACGC"
_HANA_ENZYME_MAP["BglII", "Bgl2"] = "AGATCT"
_HANA_ENZYME_MAP["EcoRV", "EcoR5"] = "GATATC"
_HANA_ENZYME_MAP["EcoRI", "EcoR1"] = "GAATTC"
_HANA_ENZYME_MAP["BamHI", "BamH1"] = "GGATCC"
_HANA_ENZYME_MAP["HindIII", "Hind3"] = "AAGCTT"
_HANA_ENZYME_MAP["KpnI", "Kpn1"] = "GGTACC"
_HANA_ENZYME_MAP["XbaI", "Xba1"] = "TCTAGA"
_HANA_ENZYME_MAP["XhoI", "Xho1"] = "CTCGAG"
_HANA_ENZYME_MAP["SacI", "Sac1"] = "GAGCTC"
_HANA_ENZYME_MAP["PstI", "Pst1"] = "CTGCAG"
_HANA_ENZYME_MAP["SmaI", "Sma1"] = "CCCGGG"
_HANA_ENZYME_MAP["PvuII", "Pvu2"] = "CAGCTG"
_HANA_ENZYME_MAP["SalI", "Sal1"] = "GTCGAC"
_HANA_ENZYME_MAP["ScaI", "Sca1"] = "AGTACT"
_HANA_ENZYME_MAP["SpeI", "Spe1"] = "ACTAGT"
_HANA_ENZYME_MAP["SphI", "Sph1"] = "GCATGC"
_HANA_ENZYME_MAP["StuI", "Stu1"] = "AGGCCT"
_HANA_ENZYME_MAP["NdeI", "Nde1"] = "CATATG"
_HANA_ENZYME_MAP["NotI", "Not1"] = "GCGGCCGC"


def enzyme_validator(enzyme_arg: any) -> List[str]:
    def __upper_map_convert(enzyme: str):
        enzyme = enzyme.upper()
        # Check map exists.
        if enzyme in _HANA_ENZYME_MAP:
            enzyme = _HANA_ENZYME_MAP[enzyme]
        if isinstance(enzyme, str):
            return [enzyme]
        return enzyme

    # Check the input argument type.
    if isinstance(enzyme_arg, str):
        enzyme_arg = __upper_map_convert(enzyme_arg)
    if isinstance(enzyme_arg, list):
        # Go through each element, convert them to uppercase.
        enzyme_result = []
        for x in enzyme_arg:
            enzyme_result += __upper_map_convert(x)
        enzyme_arg = enzyme_result

    def __enzyme_valid(enzyme: str):
        if not isinstance(enzyme, str):
            raise Exception('{} is not a valid restriction sites.'.format(enzyme))
        # Force enzyme to be upper.
        for x in enzyme:
            if x not in ('A', 'C', 'T', 'G'):
                raise Exception('Restriction site "{}" contains invalid base pair "{}"'.format(enzyme, x))

    for x in enzyme_arg:
        __enzyme_valid(x)
    return enzyme_arg
