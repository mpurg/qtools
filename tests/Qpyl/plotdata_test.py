#########################
# py.test test functions
#########################

from __future__ import absolute_import
import re
import pytest

from Qpyl.plotdata import PlotData, PlotDataError
from Qpyl.plotdata import PlotDataJSONEncoder, PlotDataJSONDecoder

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def round_json(inps):
    # round all numbers to 10 digits to float errors
    pattern = r'(-?\d+\.\d+)'
    return re.sub(pattern,
                  lambda m: "{!r}".format(round(float(m.group(1)), 10)), inps)

class TestPlotData:
    def test_encoder(self):
        pl = PlotData("random_title", "bar", "xlab", "ylab", "zlab")
        pl.add_subplot("1.subplot", [1,2,3,4], [2,4,8,16],
                       [2,4,8,16], [1,2,3,4])
        pl.add_subplot("2.subplot", [1,2,3,4], [3,4,5,6],
                       [2,3,4,5], [1,1,1,1])
        jsonenc = PlotDataJSONEncoder(indent=2, separators=(',', ': '))
        ref_data = open("data/PlotData.json").read()
        assert ref_data == round_json(jsonenc.encode(pl))


    def test_decoder(self):
        # decode json to PlotData
        qaf_data = open("data/qaf.PlotData.json").read()
        jsondec = PlotDataJSONDecoder()
        plots = jsondec.decode(qaf_data)
        assert len(plots) == 44

        dgde = plots["dgde"]
        assert dgde.xlabel == "E1-E2  [kcal/mol]"
        assert dgde.ylabel == "Free energy  [kcal/mol]"
        assert dgde.zlabel == None
        assert dgde.title == "Free-energy profile (bin-averaged, norm.)"
        assert dgde.plot_type == "line"

        sub1 = dgde.subplots["data/qfep.out.1"]
        assert is_close(sub1["ydata"][0], 5.41)
        assert is_close(sub1["xdata"][-1], 215.47)


    def test_exporting(self):
        # export PlotData object to xmgrace output
        qaf_data = open("data/qaf.PlotData.json").read()
        jsondec = PlotDataJSONDecoder()
        plots = jsondec.decode(qaf_data)
        sampling_data = plots["pts_egap"]

        ref_values = open("data/pts_egap.agr").read()
        assert sampling_data.export_grace() == ref_values

