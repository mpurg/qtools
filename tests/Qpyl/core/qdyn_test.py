#########################
# py.test test functions
#########################


from __future__ import absolute_import
from __future__ import print_function
import pytest
from Qpyl.core.qdyn import QDynInput, QDynInputError
from Qpyl.core.qdyn import QDynOutput, QDynOutput
import six
from six.moves import range
from six.moves import zip

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


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



class TestQDynOutput:

    @pytest.fixture(scope='session')
    def qdo1(self):
        return QDynOutput("data/qdyn5.log")

    @pytest.fixture(scope='session')
    def qdo2(self):
        return QDynOutput("data/qdyn6.log")

    def test_header(self, qdo1):
        # no point in having multiple tests for this
        assert qdo1.header.qdyn_version == "5.10.8"
        assert is_close(qdo1.header.stepsize, 1.0)
        assert qdo1.header.offdiagonals[1] == ('1', '9')
        assert qdo1.header.top_file == "wat.top"
        assert qdo1.header.fep_file == "wat.fep"
        assert qdo1.header.nstates == 2
        assert qdo1.header.nsteps == 100
        assert is_close(qdo1.time_end, 0.1)
        assert qdo1.time_unit == "ps"

    def test_energies(self, qdo1, qdo2):
        # at step 16
        refs = {\
        "solute": [0.016, 0.00, 0.00, 6.73, 17.30, -1.50, 0.05],
        "solvent": [0.016, -1602.97, -831.03, 0.0, 0.0, 0.0, 0.0],
        "solute_solvent": [0.016, -0.09, 0.34],
        "LRF": [0.016, -3.11],
        "Q_atom": [0.016, -191.78, 83.74, -337.75, 31.05, 1.14, 0.0],
        "restraints": [0.016, -542.1, 0.0, -584.91, 42.78, 0.0, 0.03],
        "SUM": [0.016, -3319.31, -3369.97, 50.66]
        }

        for k, ref in six.iteritems(refs):
            print("Testing ", k)
            val = getattr(qdo1, "data_E_"+k).get_rows()[15]
            val2 = getattr(qdo1, "data_E_"+k).get_rows()[15]
            val3 = getattr(qdo1, "data_E_"+k).get_rows()[16]
            assert str(val) == str(ref)
            assert str(val2) == str(ref)
            assert str(val3) != str(ref)

    def test_qenergies(self, qdo1, qdo2):
        # at step 16
        refs = {\
        "Q": [[0.016, 0.5, -95.89, 51.9], [0.016, 0.5, -165.76, 87.73]],
        "prot": [[0.016, 0.5, -3.89, 0.45], [0.016, 0.5, -0.87, 0.46]],
        "wat": [[0.016, 0.5, -53.07, 12.57], [0.016, 0.5, -64.09, 14.38]],
        "surr": [[0.016, 0.5, -56.96, 13.02], [0.016, 0.5, -64.96, 14.84]],
        "SUM": [[0.016, 0.5, -446.89, 0.0], [0.016, 0.5, -380.3, 0.0]],
        "any": [[0.016, 0.5, -152.85, 64.92, -374.88, 15.65, 0.27, 0.0],
                [0.016, 0.5, -230.72, 102.56, -300.61, 46.45, 2.02, 0.0]]
        }

        for k, ref in six.iteritems(refs):
            for i in range(2): # two states
                print("Testing ", k, ", state ", i)
                val = getattr(qdo1, "data_EQ_"+k)[i].get_rows()[15]
                val2 = getattr(qdo1, "data_EQ_"+k)[i].get_rows()[15]
                val3 = getattr(qdo1, "data_EQ_"+k)[i].get_rows()[21]
                assert str(val) == str(ref[i])
                assert str(val2) == str(ref[i])
                assert str(val3) != str(ref[i])

    def test_temps(self, qdo1, qdo2):
        # note that zeroth step was removed from qdyn5.log manually
        val1 = qdo1.data_temp.get_rows()[15]
        val2 = qdo2.data_temp.get_rows()[15]
        ref = [0.016, 4.0, 4.0, 4.0, 4.0]
        for v1, v2, r in zip(val1, val2, ref):
            assert is_close(v1, r)
            assert is_close(v2, r)

