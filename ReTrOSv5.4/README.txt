ReTrOS: Reconstructing Transcription Open Software (2017),
Giorgos Minas; Hiroshi Momiji; Dafyd Jenkins; Maria J. Costa; David A. Rand; Barbel Finkenstadt.
Copyright (C) 2017.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


*******


CONTENTS

<  m files  >
ConfEnvRatesKR.m	ksrlin.m		trapeziumHM.m
backcalculationKR.m	optimalBandWidth.m
calcdegrateHM.m		runConfEnvRatesKR.m
drawDegRates.m		runReTrOS.m (TOP file)

<  mat files  >
data17.mat
data27.mat

<  txt files  >
README.txt (this file)	license.txt
testdata.txt (sample data - three time series in one file)

<  sample input/output files in the "testData" folder  >

mockBCdata.txt (3 mock time series)

In the "output_examples" subfolder:

mockBCdata.png (runtime figure, mockBCdata.txt)
mockBCdata_outC.txt (output, mockBCdata.txt)
mockBCdata_outD.txt (output, mockBCdata.txt)
mockBCdata_outM.mat (output, mockBCdata.txt)


*******


HOW TO CREATE AN INPUT DATA FILE

The files need to be in the "Windows Formatted Text" format. The sample file (mockBCdata.txt) was created and saved in Excel.

Read with "mockBCdata.txt" -- The first five columns are reserved for annotations, while numerical data are stored in column 6 and onwards, where the first row contain the time, and the following rows contain the measurements.

There are further rules to follow: (1) Columns 1--5 should contain only texts; (2) With the header "TEMP" on the first row, the second column contains measurement temperatures in such a text format as "12C" (see testdata.txt). This column needs to be filled (with dummy temperatures, if necessary, e.g.  when a user provides the information about the degradation rates e.g. as "gLUCDEG = [...]". (3) Column 6 should contain sample names.


*******


HOW TO USE

(0)  Edit the top file "runReTrOS.m", if necessary.  It includes four switches.  They are described in the file, with possible choices and recommendations.

(1)  In Matlab's command window, do

>> runReTrOS;

Then, a window pops up to ask you to choose an input (multiple) time series file.

(2)  When you chose an input file, calculation starts, with a result figure window popping up (see BCtest.png).  For each time series, the figure window is updated for each time series, and  you will see something like the following in the command window:

>> runReTrOS
  1/  1: WT:LHY(XX) in calc...
1.500000	21278265873.880745
...

opt. band width = 1.500000
   9: done
  19: done
  ...
  89: done
  99: done

The first half is monitoring the bandwidth optimisation, while the second half is monitoring boot-strapping.

(3)  After the completion of all calculations, you are asked to give output file name.  Type it (XXXXXXX) WITHOUT file extension.  You will then have three output files: XXXXXXX_outC.txt; XXXXXXX_outD.txt; XXXXXXX_outM.mat.  Text files contain the median transcription activity, in the continuous (C) and in the discrete (D) time, while more information including confidence envelope is stored in the mat file.

(4)  The average and the variance (NOT std dev!!) of native mRNA degradation rate can be given on L22 of "runConfEnvRatesKR.m".

(5)  You may also increase the number of boot-strapping from "99" (default) say to 999, on L24 of "ConfEnvRatesKR.m".

(6)  The user of an earlier version of Matlab may need to change L42 of "optimalBandWidth.m"

from  [~,jhmin] = min(mse);
to    [dummy,jhmin] = min(mse);
