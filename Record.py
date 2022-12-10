import copy


class PatientRecord(object):
    def __init__(self, patient_id=None, correlation=None, recist_score=None):
        self.patient_id = patient_id
        self.recist_score = recist_score
        if correlation is None:
            self.correlation = []
        else:
            self.correlation = copy.copy(correlation)

        self.drug_synonyms_list = None


class CellineData(object):
    def __init__(self, index=0, drug_name='', sampleid='',sensitivity=None):
        self.index = index
        self.drug_name = drug_name
        self.sensitivity = sensitivity
        self.sampleid = sampleid


class FinalRecord(object):
    def __init__(self, patient_id: str, drug_name: str, celline: str, correlation: float,
                 sensitivity: dict, recist_score=None):
        self.patient_id = patient_id
        self.drug_name = drug_name
        self.celline = celline
        self.correlation = correlation

        self.ic50 = sensitivity['ic50_published']
        self.auc = 1-sensitivity['aac_published']
        self.amax = sensitivity['amax_published']
        self.HS = sensitivity['HS']
        self.E_inf = sensitivity['E_inf']
        self.EC50 = sensitivity['EC50']

        self.recist_score = recist_score

