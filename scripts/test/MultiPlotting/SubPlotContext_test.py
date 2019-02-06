# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest

from MultiPlotting.subplot.subplot_context import subplotContext
from mantid import plots

try:
    from unittest import mock
except ImportError:
    import mock


class line(object):

    def __init__(self):
        self.label = "test"

    def get_label(self):
        return self.label

    def get_color(self):
        return "red"

    def get_marker(self):
        return "star"

    def remove(self):
        return


class label(object):

    def __init__(self, name, protected):
        self.text = name
        self.protected = protected


def errors():
    return tuple([line(), [line()], [line()]])


class SubPlotContextTest(unittest.TestCase):

    def setUp(self):
        name = "test"
        self.subplot = mock.MagicMock()
        self.context = subplotContext(name, self.subplot)

    def test_add_line_no_erros(self):
        ws = mock.MagicMock()
        with mock.patch("mantid.plots.plotfunctions.plot") as patch:
            patch.return_value = tuple([line()])
            self.context.addLine(ws, 3)
            self.assertEquals(patch.call_count, 1)
            patch.assert_called_with(self.subplot, ws, specNum=3)

    def test_add_line_errors(self):
        ws = mock.MagicMock()
        self.context.change_errors(True)
        lines = line()
        with mock.patch("mantid.plots.plotfunctions.plot") as plot:
            plot.return_value = tuple([lines])
            with mock.patch("mantid.plots.plotfunctions.errorbar") as patch:
                patch.return_value = errors()
                self.context.addLine(ws, 3)
                self.assertEquals(plot.call_count, 1)
                self.assertEquals(patch.call_count, 1)
                patch.assert_called_with(
                    self.subplot,
                    ws,
                    specNum=3,
                    label=lines.get_label())

    def test_redraw_errors(self):
        ws = mock.MagicMock()
        self.context.change_errors(True)
        lines = line()
        # add first line
        with mock.patch("mantid.plots.plotfunctions.plot") as plot:
            plot.return_value = tuple([lines])
            with mock.patch("mantid.plots.plotfunctions.errorbar") as patch:
                patch.return_value = errors()
                self.context.addLine(ws, 3)
                self.assertEquals(plot.call_count, 1)
                self.assertEquals(patch.call_count, 1)
                # redraw
                self.context.redraw(lines.get_label())
                self.assertEquals(patch.call_count, 2)
                patch.assert_called_with(
                    self.subplot,
                    ws,
                    specNum=3,
                    color=lines.get_color(),
                    marker=lines.get_marker(),
                    label=lines.get_label())

    def test_redraw_no_errors(self):
        ws = mock.MagicMock()
        with mock.patch("mantid.plots.plotfunctions.plot") as patch:
            lines = line()
            patch.return_value = tuple([lines])
            self.context.addLine(ws, 3)
            self.assertEquals(patch.call_count, 1)
            patch.assert_called_with(self.subplot, ws, specNum=3)
            # redraw
            self.context.redraw(lines.get_label())
            self.assertEquals(patch.call_count, 2)
            patch.assert_called_with(
                self.subplot,
                ws,
                specNum=3,
                color=lines.get_color(),
                marker=lines.get_marker(),
                label=lines.get_label())

    def test_change_errors(self):
        self.context._lines = {"one": 1, "two": 2, "three": 3}
        self.context.redraw = mock.MagicMock()
        self.context.change_errors(True)
        self.assertEquals(self.context.redraw.call_count, 3)

    def test_change_auto(self):
        self.context._lines = {"one": 1, "two": 2, "three": 3}
        self.context._subplot.autoscale = mock.MagicMock()
        self.context.change_auto(True)
        self.assertEquals(self.context._subplot.autoscale.call_count, 3)

    def test_vlines(self):
        self.context._labelObjects = {
            "one": label("one", True), "two": label("two", False), "three": label("three", False)}
        self.context._vLines = {
            "two": mock.MagicMock(), "four": mock.MagicMock()}
        result = self.context.vlines
        self.assertEquals(["three", "two", "four"], result)

if __name__ == "__main__":
    unittest.main()
