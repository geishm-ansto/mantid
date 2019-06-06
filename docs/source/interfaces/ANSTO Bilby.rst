.. _ANSTO_Bilby-ref:

ANSTO BILBY
============

.. image::  ../images/sans_isis_v2_whole_gui.png
   :align: right
   :width: 800px

.. contents:: Table of Contents
  :local:

Interface Overview
------------------

This interface is used to reduce ANSTO data for BILBY. 
The interface can be accessed from the main menu of MantidPlot, in *Interfaces → ANSTO → BILBY*.
This interface takes it appearance and behaviour from the ISIS SANS interface but deviates from 
the ISIS implementation because of the specific processing chain implemented for BILBY. 

Runs
----

.. _Runs:

The *Runs* tab is the entry-point. It allows the user to specify the user file and
batch file. Alternatively the user can set the data sets for reduction manually on the data table.
In addition it allows the user to specify her preferred save settings (see more below). The actual
parameters which are loaded from the user file are accessible from the sub-tabs on the Settings_ tab.

Data Table
^^^^^^^^^^

.. _RunsDataTable:

.. image::  ../images/sans_isis_v2_run_tab_data_table.png
   :align: center
   :width: 800px

:Insert:
	Adds a row after the currently selected row. 
:Delete:
	Deletes a selected row.   
:Copy:
	Creates a copy of the selected rows. 
:Cut:
	Cuts the selected rows.  
:Paste:
	Pastes rows from the clipboard.    
:Clear:
	Clears the entries from selected rows. It however does not the delete the rows themselves.
:Export Table:
	Export the table to a csv file that can be reloaded later.  
:Output Folder:
	A sub-folder in the save directory used to save results.  
:Process Selected:
	Process the selected rows from the table. The scheduled rows are highlighted. 
:Process All:
	Process all row regardless of selection.  
:Cancel Processing:
	A sub-folder in the save directory used to save results.
	

Columns
^^^^^^^
:Sample:
	Mandatory scatter data file to use.
:T Empty Beam:
	Mandatory empty transmission data. Mandatory field.
:Trans. Sample:
	Mandatory transmission data file to use. 
:Thickness:
	Sample thickness in cm.
:Blocked Beam:
	Optional blocked beam data file to use.
:Start Time:
	Process event pulses after start of test in seconds. Skip check if blank.
:End Time:
	Process event pulses before end time from start of in seconds. Skip check if blank. 
:Sample Mask:
	Mask file applied to sample data. If blank uses the mask file listed in the settings tab if available.
:Trans. Mask:
	Mask file applied to transmission data. If blank uses the mask file listed in the settings tab if available.
:Suffix:
	Suffix and description are appended to the output file name.
:Description:
	Suffix and description are appended to the output file name.


Settings
--------

.. _Settings:

.. image::  ../images/sans_isis_v2_general_tab_whole.png
   :align: right
   :width: 800px

The Settings tab and its sub-tabs allow for manipulating and inspecting the reduction parameters which were
initially set through loading a user file.  Currently there are three sub-tabs:

:Reduction parameters:
   Settings for the wavelength and momentum transfer binning.
   
      AS to add ...

:Transmission parameters: 
   Transmission fit options and binning parameters.
   
      AS to add ...

:Advanced parameters: 
   Contains a mix of processing options to account for gravity, solid angle, wide angle and blocked beam corrections along with the default mask files and radius and wave cut parameters.
   
      AS to add ...


.. categories:: Interfaces SANS
