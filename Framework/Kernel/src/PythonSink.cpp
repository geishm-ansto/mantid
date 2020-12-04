// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidKernel/PythonSink.h"

#include <boost/iostreams/categories.hpp> // sink_tag
#include <iosfwd>                         // streamsize
#include <boost/format.hpp>

#include <Python.h>

// test_ostream  std::cout

std::streamsize pysys_stdout_sink::write(const char *s, std::streamsize n) {
  // PySys_WriteStdout truncates to 1000 chars
  static const std::streamsize MAXSIZE = 1000;

  std::streamsize written = std::min(n, MAXSIZE);
  PySys_WriteStdout( (boost::format("%%.%1%s") % written).str().c_str(), s );
  // std::cout << s;

  return written;
}
