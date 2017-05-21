#########################
# py.test test functions
#########################


import pytest
from Qpyl.core.qdyninp import QDynInput, QDynInputError

class TestQDynInput: 

    def test_parse_input(self):
        # parse and write back, see if it's the same
        qdis = open("data/qdyn.inp.1").read()
        assert QDynInput(qdis).get_string() == qdis

    def test_update_parms(self):
        # update the params, see what happens
        qdis = open("data/qdyn.inp.1").read()
        qdi = QDynInput(qdis, parameters = {"md": {"temperature":100}})
        assert qdi.parameters["md"]["temperature"] == 100
        qdi.update(input_string=qdis)
        assert int(qdi.parameters["md"]["temperature"]) == 300
        qdi.update(parameters = {"md": {"temperature":100}})
        assert qdi.parameters["md"]["temperature"] == 100

    def test_parse_input_fail(self):
        # fail on typos
        qdis = open("data/qdyn.inp.1").read()
        with pytest.raises(QDynInputError):
            QDynInput(qdis.replace("temperature", "tmperature"))

        with pytest.raises(QDynInputError):
            QDynInput(qdis.replace("[lambdas]\n", ""))

        with pytest.raises(QDynInputError):
            QDynInput(qdis.replace("off", "whatisthis?"))

        with pytest.raises(QDynInputError):
            QDynInput(qdis.replace("300", "300_"))

        # FEP and lambda must exist at the same time
        with pytest.raises(QDynInputError):
            qdi = QDynInput(qdis)
            del qdi.parameters["lambdas"]
            qdi.get_string()

        # random_seed and initial_vel must exist at the same time
        with pytest.raises(QDynInputError):
            qdi = QDynInput(qdis)
            del qdi.parameters["md"]["random_seed"]
            qdi.get_string()

        # energy file and interval must be both specified
        with pytest.raises(QDynInputError):
            qdi = QDynInput(qdis)
            del qdi.parameters["intervals"]["energy"]
            qdi.get_string()




