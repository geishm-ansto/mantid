# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

# pylint: disable=too-few-public-methods, invalid-name

from __future__ import (absolute_import, division, print_function)
import os
import re
import h5py as h5
import tarfile
import tempfile
from mantid.kernel import DateAndTime
from sans.common.enums import (SANSInstrument, SANSFacility, FileType, SampleShape)
from sans.common.file_information import (SANSFileInformation)

def is_ansto_file_name(file_name):
    """
    Confirms that the filename matches the ANSTO format. 

    :param file_name: the full file name to an existing file.
    :return: bool
    """
    patterns = [r'[A-Z]{3}[0-9]{7}\.tar',
               r'[A-Z]{3}[0-9]{7}\.nx\.hdf']
    base = os.path.basename(file_name)
    for pattern in patterns:
        if re.match(pattern, base):
            return True
    return False

def get_run_number_for_ansto_file(file_name, prefix, extn):
    """
    Extracts the 7 digit run number from the file name. If the file_name
    does not match the expected format it returns -1 

    :param file_name: the full file name to an existing raw file.
    :return: run
    """
    pattern = prefix +  r'([0-9]{7}?)' + extn
    base = os.path.basename(file_name)
    m = re.match(pattern, base)
    return int(m.group(1)) if m else -1

class BILBYFileInformation(SANSFileInformation):
    def __init__(self, sample_file, transmission_file):
        super(BILBYFileInformation, self).__init__(sample_file)
        
        self._number_of_periods = 1
        self._run_number = get_run_number_for_ansto_file(sample_file, prefix='BBY', extn='.tar')

        # get the parameters from the sample file, for each key pass the 
        # hdf tag and a default value, if there is no default value an 
        # exception will be raised if it is missing
        tags = {'start_time': ('/entry1/start_time', None),
                'frame_source': ('/entry1/instrument/detector/frame_source', ['external']),
                'aperture': ('/entry1/sample/sample_aperture', None),
                'duration': ('/entry1/instrument/detector/time', None),
                }
        values = self._extract_hdf_params(sample_file, tags)

        self._instrument = SANSInstrument.BILBY
        self._facility = SANSFacility.ANSTO
        self._date = DateAndTime(values['start_time'][0])
        frame_source = values['frame_source'][0]
        self._is_event_mode = frame_source.lower() == 'external'
        self._aperture = values['aperture'][0]
        self._duration = values['duration'][0]

        # now collect the transmission file
        tags = {'att_pos': ('/entry1/instrument/att_pos', None)}
        values = self._extract_hdf_params(transmission_file, tags)
        self._att_pos = int(round(values['att_pos'][0]))

        # set to default values
        self._thickness = 1.
        self._height = 1.
        self._width = 1.
        self._shape = SampleShape.Disc

    def _extract_hdf_params(self, file_path, tags):

        # extract hdf from tar file and complete loading 
        try:
            temp_dir = tempfile.mkdtemp()
            hdf_path = None            
            with tarfile.open(file_path, mode='r') as tf:
                for name in tf.getnames():
                    if name.endswith('.nx.hdf'):
                        tf.extract(name, temp_dir)
                        hdf_path = os.path.normpath(os.path.join(temp_dir, name))
                        break
            if hdf_path is None:
                raise RuntimeError('{} tar file is missing hdf file.'.format(file_path))

            # Setup the parameters from the hdf file, if it is missing it returns an
            # empty list
            values = {}
            try:
                with h5.File(hdf_path, 'r') as fp:
                    for k, (hdf_tag, def_val) in tags.items():            
                        try:
                            values[k] = fp[hdf_tag].value
                        except KeyError:
                            if def_val is None:
                                raise RuntimeError('{} - not included in the hdf file'.format(hdf_tag))
                            else:
                                values[k] = def_val
            except IOError:
                raise RuntimeError('{} - contained hdf parameter error missing hdf file.'.format(file_path))

        finally:
            os.remove(hdf_path)
            os.rmdir(temp_dir)

        return values

    def get_file_name(self):
        return self._full_file_name

    def get_instrument(self):
        return self._instrument

    def get_facility(self):
        return self._facility

    def get_date(self):
        return self._date

    def get_number_of_periods(self):
        return self._number_of_periods

    def get_run_number(self):
        return self._run_number

    def get_type(self):
        return FileType.ANSTOTar

    def is_event_mode(self):
        return self._is_event_mode

    def is_added_data(self):
        return False

    def get_height(self):
        raise RuntimeError('Height is not available from BILBY files')

    def get_width(self):
        raise RuntimeError('Width is not available from BILBY files')

    def get_thickness(self):
        raise RuntimeError('Thickness is not available from BILBY files')

    def get_shape(self):
        raise RuntimeError('Shape is not available from BILBY files')

    def get_duration(self):
        return self._duration

    def get_att_pos(self):
        return self._att_pos
