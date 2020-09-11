#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# MIT License
# 
# Copyright (c) 2018  Miha Purg <miha.purg@gmail.com>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
#

"""
This module contains the class for custom data structure PlotData
used for storing various analysis data.
It also contains custom JSON encoder and decoder classes 
which allow saving the objects in JSON format.
"""

from __future__ import absolute_import, division, unicode_literals
import six
from six.moves import zip

import json
import sys
from collections import OrderedDict as ODict


class PlotDataJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, PlotData):
            return ODict([("__type__",  "PlotData"),
                          ("title",     obj.title),
                          ("plot_type", obj.plot_type),
                          ("xlabel",    obj.xlabel),
                          ("ylabel",    obj.ylabel),
                          ("subplots",  obj.subplots)])
        else:
            return json.JSONEncoder.default(self, obj)


class PlotDataJSONDecoder(json.JSONDecoder):
    def __init__(self):
        if sys.version_info < (2, 7):
            # object_pairs_hook is supported only in version 2.7
            raise PlotDataError("You need python 2.7 or later to run this"
                                "script, sorry (it's json's fault)!")
        super(PlotDataJSONDecoder,
                self).__init__(object_pairs_hook=self.decode_plotdata)

    def decode_plotdata(self, d):
        d = ODict(d)
        if "__type__" not in d:
            return d
        t = d["__type__"]
        if t == "PlotData":
            pd = PlotData(d["title"], d["plot_type"], d["xlabel"], d["ylabel"])
            pd.subplots = d["subplots"]
            return pd
        else:
            return d

class PlotDataError(Exception):
    pass

class PlotData(object):
    def  __init__(self, title, plot_type="line", 
                  xlabel=None, ylabel=None, zlabel=None):
        self.title = title
        plot_types = ["line", "bar", "scatter", "wireframe"]
        if plot_type not in plot_types:
            raise ValueError("'plot_type' {} not supported. "
                             "Try one of these instead: {}"
                             "".format(plot_type, ",".join(plot_types)))
        self.plot_type = plot_type
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel
        self.subplots = ODict()

    def add_subplot(self, label, xdata, ydata, zdata=None, yerror=None):
        self.subplots[label] = ODict([("xdata", xdata),
                                      ("ydata", ydata),
                                      ("zdata", zdata),
                                      ("yerror", yerror)])


    # TODO: clean this code, looks almost as bad as Javascript
    def export_grace(self):
        if self.plot_type == "line":
            typ = "xy"
        elif self.plot_type == "bar":
            typ = "bar"
        elif self.plot_type == "scatter":
            typ = "xy"
        elif self.plot_type == "wireframe":
            raise PlotDataError("Cannot export wireframe data to grace...")

        set_config = ""
        sets = ""
        for i, (label, sp) in enumerate(six.iteritems(self.subplots)):
            # create this:
            # @s0 legend "rep_000"
            # @s1 legend "rep_001" ...
            set_config += "@s{} legend \"{}\" \n".format(i, label) 
            if typ == "bar":
                # don't show the line in bar plots
                set_config += "@s{} line type 0 \n".format(i)

            # add the data
            if not sp["yerror"] or len(sp["yerror"]) != len(sp["xdata"]):
                yerror=["" for x in sp["xdata"]]
            else:
                yerror = sp["yerror"]
                typ = typ + "dy"
            for x, y, dy in zip(sp["xdata"], sp["ydata"], yerror):
                sets += "{} {} {}\n".format(x, y, dy)
            sets += "&\n"
        
        return """#
@type {typ}
@title "{title}"
@xaxis label "{xlabel}"
@yaxis label "{ylabel}"
{set_config}
{sets}

""".format(typ=typ, title=self.title, xlabel=self.xlabel, ylabel=self.ylabel, set_config=set_config, sets=sets)


