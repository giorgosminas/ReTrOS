function [] = ReTrOSGUI

% ReTrOS, Reconstructing Transcriptional Activity from Gene and Protein Expression Data
% by 
% Copyright (C) 2016 Giorgos Minas, Hiroshi Momiji, Dafyd Jenkins, Maria
% Costa, David Rand, Barbel Finkenstadt
% and University of Warwick
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


currWarning = warning('off','all');
%SET UP DEFAULT VALUES FIRST:
%General:
%degradation rates
%
%Smooth specific:
%bootstrap samples
%kernel type
%polynomial degree
%bandwidth (?)
%time resolution
%
%Switch specific:
%MCMC iterations
%burn-in
%posterior samples
%plot samples
%baseline factor

detrendNone_str = java.lang.String('none');
detrendLinear_str = java.lang.String('linear');
detrendLogLinear_str = java.lang.String('log-linear');
detrendArray = javaArray('java.lang.String',3);
detrendArray(1) = detrendNone_str;
detrendArray(2) = detrendLinear_str;
detrendArray(3) = detrendLogLinear_str;

mRNA_degRate1 = 0;
mRNA_degRate2 = 0;

protein_degRate1 = 0;
protein_degRate2 = 0;

degRateDistArray = javaArray('java.lang.String',2);
degRateDistArray(1) = java.lang.String('Normal');
degRateDistArray(2) = java.lang.String('Gamma');
degRateDistIdx = 1; %Gamma

degRateDefaultArray = javaArray('java.lang.String',3);
degRateDefaultArray(1) = java.lang.String('user input');
degRateDefaultArray(2) = java.lang.String('Arabidopsis'); %Gamma,mRNA: a=1.4828, b=0.1823 (Plant Cell. 2007 Nov;19(11):3418-36), protein: a=3.3765, b=0.0822 (Mol Cell Proteomics. 2012 June; 11(6): M111.010025)
degRateDefaultArray(3) = java.lang.String('luc reporter');%Normal,ROBUST
degRateActualDefaultIdx = 1; %Arabidopsis
degRateReporterDefaultIdx = 2; %luc reporter

mRNA_str = java.lang.String('mRNA');
protein_reporter_str = java.lang.String('protein reporter');

%DEFAULT DEG RATES
ARAB_MRNA_ALPHA = 1.4828;
ARAB_MRNA_BETA = 0.1823;
ARAB_PROTEIN_ALPHA = 3.3765;
ARAB_PROTEIN_BETA = 0.0822;

LUC_MRNA_MEAN = 2.2823; % calculated at 17oC
LUC_MRNA_SDEV = 0.4569;
LUC_PROTEIN_MEAN = 0.1295;
LUC_PROTEIN_SDEV = 0.0104;

smoothBootstrapSamplesDefault = 99;
smoothBandwidthMinDefault = 1;
smoothBandwidthMaxDefault = 5;
smoothBandwidthIncDefault = 0.1;
smoothTimeResolutionDefault = 0.1;
smoothKernelDegreeArray = javaArray('java.lang.String',3);
smoothKernelDegreeArray(1) = java.lang.String('1');
smoothKernelDegreeArray(2) = java.lang.String('2');
smoothKernelDegreeArray(3) = java.lang.String('3');
smoothKernelDegreeDefaultIdx = 1; %2
smoothKernelTypeArray = javaArray('java.lang.String',3);
smoothKernelTypeArray(1) = java.lang.String('Gaussian');
smoothKernelTypeArray(2) = java.lang.String('Epanechnikov');
smoothKernelTypeArray(3) = java.lang.String('Triweight');
smoothKernelTypeDefaultIdx = 0; %Gaussian

previewFigH = [];

switchIterationsDefault = 100000;
switchRegressionMethodArray = javaArray('java.lang.String',2);
switchRegressionMethodArray(1) = java.lang.String('least squares');
switchRegressionMethodArray(2) = java.lang.String('weighted least squares');
switchRegressionMethodDefaultIdx = 0; %least squares
switchBurnInDefault = 0.25;
switchPosteriorSampleNumDefault = switchIterationsDefault;
switchBaselineFactorDefault = 1;
switchPlotSampleNumDefault = 5000;
switchMinimumSwStrDefault = 0.2;
switchMaxSwDeviationFactorDefault = 0.5;

randomSeedDefault = 0;

acceptedInputs = false;
hidFig = figure('visible','off');
frame = javaObjectEDT('javax.swing.JFrame');
hFrame = handle(frame,'CallbackProperties');
set(hFrame,'WindowClosingCallback',@closeFrame);
% javaMethodEDT('setDefaultCloseOperation',frame,javax.swing.JFrame.DISPOSE_ON_CLOSE);
javaMethodEDT('setDefaultCloseOperation',frame,javax.swing.JFrame.DO_NOTHING_ON_CLOSE);
javaMethodEDT('setResizable',frame,false);

RETROS_RUNNING = false;
RETROS_CANCELLED = false;

dataStruct = [];
dataStruct.hasData = false;
dataStruct.hasAlgorithmParameters = false;
dataSelectUI(frame);

% %TEMP
% if ishghandle(hidFig)
%     close(hidFig);
% end

waitfor(hidFig);
warning(currWarning);

% if acceptedInputs
%     for i = 1:10
%         value = 1;
%         for j = 1:100000000
%             value = value + j;
%         end
%         value;
%     end
% end

    function [] = closeFrame(varargin)
        javaMethodEDT('dispose',frame);
        if ishghandle(previewFigH)
            close(previewFigH);
        end
        clear global RETROS_RUNNING;
        close(hidFig);
    end

    function [] = dataSelectUI(frame)
        
        javaMethodEDT('setTitle',frame,'ReTrOS: Select Data');
        contentPane = javaMethodEDT('getContentPane',frame);
        javaMethodEDT('removeAll',contentPane);
        
        % dataStruct = [];
        
        inputPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.BorderLayout'));
        
        selectedFile = 'No file selected';
        isValidFile = false;
        
        %set up file select components (NORTH)
        filePanel = javaObjectEDT('javax.swing.JPanel');
        filePanelBorder = javaMethodEDT('createTitledBorder','javax.swing.BorderFactory','Data file selection');
        javaMethodEDT('setBorder',filePanel,filePanelBorder);
        fileSelectButton = javaObjectEDT('javax.swing.JButton','Select data file');
        fileSelectButtonH = handle(fileSelectButton,'CallbackProperties');
        set(fileSelectButtonH,'ActionPerformedCallback',@fileSelectButtonEvent);
        p = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridBagLayout'));
        javaMethodEDT('add',p,fileSelectButton);
        javaMethodEDT('add',filePanel,p,java.awt.BorderLayout.WEST);
        fileSelectField = javaObjectEDT('javax.swing.JTextField',selectedFile,50);
        javaMethodEDT('add',filePanel,fileSelectField,java.awt.BorderLayout.CENTER);
        
        %set up data preview components (CENTER)
        previewPanel = javaObjectEDT('javax.swing.JPanel');
        previewBorder = javaMethodEDT('createTitledBorder','javax.swing.BorderFactory','Data preview');
        javaMethodEDT('setBorder',previewPanel,previewBorder);
        previewTable = javaObjectEDT('javax.swing.JTable',1,1);
        javaMethodEDT('setPreferredScrollableViewportSize',previewTable,javaObjectEDT('java.awt.Dimension',850,300));
        javaMethodEDT('setAutoResizeMode',previewTable,javax.swing.JTable.AUTO_RESIZE_OFF);
        javaMethodEDT('setFillsViewportHeight',previewTable,true);
        javaMethodEDT('setShowGrid',previewTable,true);
        javaMethodEDT('setReorderingAllowed',javaMethodEDT('getTableHeader',previewTable),false);
        tableScrollPane = javaObjectEDT('javax.swing.JScrollPane',previewTable);
        javaMethodEDT('add',previewPanel,tableScrollPane);
        % javaMethodEDT('add',previewPanel,previewTable);
        
        
        %set up properties components (SOUTH)
        propertiesPanel = javaObjectEDT('javax.swing.JPanel');
        propertiesBorder = javaMethodEDT('createTitledBorder','javax.swing.BorderFactory','Data properties');
        javaMethodEDT('setBorder',propertiesPanel,propertiesBorder);
        
        % data properties:
        %combineRows
        %numcolumns
        %numRows
        %timepoints
        %timescale
        %replicates
        %protein / mRNA
        
        
        expTypeArray = javaArray('java.lang.String',2);
        expTypeArray(1) = mRNA_str;
        expTypeArray(2) = protein_reporter_str;
        
        timescaleLabel = javaObjectEDT('javax.swing.JLabel','Timescale:');
        timescaleField = javaObjectEDT('javax.swing.JTextField','',60);
        timepointsLabel = javaObjectEDT('javax.swing.JLabel','# timepoints:');
        timepointsField = javaObjectEDT('javax.swing.JTextField','',4);
        replicatesLabel = javaObjectEDT('javax.swing.JLabel','# replicates:');
        replicatesField = javaObjectEDT('javax.swing.JTextField','',4);
        %         datapointsLabel = javaObjectEDT('javax.swing.JLabel','# datapoints per profile:');
        %         datapointsField = javaObjectEDT('javax.swing.JTextField','',3);
        profilesLabel = javaObjectEDT('javax.swing.JLabel','# unique profiles:');
        profilesField = javaObjectEDT('javax.swing.JTextField',4);
        %         combineReplicatesBox = javaObjectEDT('javax.swing.JCheckBox','Combine non-unique names?',false);
        %         combineRepsH = handle(combineReplicatesBox,'CallbackProperties');
        %         set(combineRepsH,'ActionPerformedCallback',@combineReplicatesEvent);
        profileNameLabel = javaObjectEDT('javax.swing.JLabel','Profile names:');
        profileNameBox = javaObjectEDT('javax.swing.JComboBox');
        profileNameH = handle(profileNameBox,'CallbackProperties');
        set(profileNameH,'ActionPerformedCallback',@profileNameSelectEvent);
        profileTypeLabel = javaObjectEDT('javax.swing.JLabel','Profile type:');
        profileTypeBox = javaObjectEDT('javax.swing.JComboBox',expTypeArray);
        profileTypeH = handle(profileTypeBox,'CallbackProperties');
        set(profileTypeH,'ActionPerformedCallback',@profileTypeSelectEvent);
        
        detrendLabel = javaObjectEDT('javax.swing.JLabel','Data detrend:');
        detrendBox = javaObjectEDT('javax.swing.JComboBox',detrendArray);
        previewButton = javaObjectEDT('javax.swing.JButton','Preview');
        javaMethodEDT('setEnabled',previewButton,false);
        previewButtonH = handle(previewButton,'CallbackProperties');
        set(previewButtonH,'ActionPerformedCallback',@showProfilePreview);
        
        javaMethodEDT('setEditable',timescaleField,false);
        javaMethodEDT('setEditable',timepointsField,false);
        javaMethodEDT('setEditable',replicatesField,false);
        %         javaMethodEDT('setEditable',datapointsField,false);
        javaMethodEDT('setEditable',profilesField,false);
        javaMethodEDT('setEditable',fileSelectField,false);
        %         javaMethodEDT('setEnabled',combineReplicatesBox,false);
        
        gb = javaObjectEDT('java.awt.GridBagLayout');
        javaMethodEDT('setLayout',propertiesPanel,gb);
        c = javaObjectEDT('java.awt.GridBagConstraints');
        c.fill = java.awt.GridBagConstraints.BOTH;
        c.weightx = 0.5;
        javaMethodEDT('setConstraints',gb,timescaleLabel,c);
        javaMethodEDT('add',propertiesPanel,timescaleLabel);
        c.gridwidth=java.awt.GridBagConstraints.REMAINDER;
        c.gridx = 1;
        javaMethodEDT('setConstraints',gb,timescaleField,c);
        javaMethodEDT('add',propertiesPanel,timescaleField);
        c.gridx=0;
        c.gridy=1;
        c.gridwidth=1;
        javaMethodEDT('setConstraints',gb,timepointsLabel,c);
        javaMethodEDT('add',propertiesPanel,timepointsLabel);
        c.gridx=1;
        c.gridwidth = 1;
        javaMethodEDT('setConstraints',gb,timepointsField,c);
        javaMethodEDT('add',propertiesPanel,timepointsField);
        c.gridx=2;
        javaMethodEDT('setConstraints',gb,replicatesLabel,c);
        javaMethodEDT('add',propertiesPanel,replicatesLabel);
        c.gridx=3;
        javaMethodEDT('setConstraints',gb,replicatesField,c);
        javaMethodEDT('add',propertiesPanel,replicatesField);
        c.gridx=4;
        javaMethodEDT('setConstraints',gb,detrendLabel,c);
        javaMethodEDT('add',propertiesPanel,detrendLabel);
        c.gridx=5;
        javaMethodEDT('setConstraints',gb,detrendBox,c);
        javaMethodEDT('add',propertiesPanel,detrendBox);
        c.gridx = 6;
        c.gridheight=2;
        javaMethodEDT('setConstraints',gb,previewButton,c);
        javaMethodEDT('add',propertiesPanel,previewButton);
        c.gridheight=1;
        %         javaMethodEDT('setConstraints',gb,datapointsLabel,c);
        %         javaMethodEDT('add',propertiesPanel,datapointsLabel);
        %         c.gridx=5;
        %         javaMethodEDT('setConstraints',gb,datapointsField,c);
        %         javaMethodEDT('add',propertiesPanel,datapointsField);
        c.gridx=0;
        c.gridy=2;
        javaMethodEDT('setConstraints',gb,profilesLabel,c);
        javaMethodEDT('add',propertiesPanel,profilesLabel);
        c.gridx=1;
        javaMethodEDT('setConstraints',gb,profilesField,c);
        javaMethodEDT('add',propertiesPanel,profilesField);
        c.gridx=2;
        javaMethodEDT('setConstraints',gb,profileNameLabel,c);
        javaMethodEDT('add',propertiesPanel,profileNameLabel);
        c.gridx=3;
        javaMethodEDT('setConstraints',gb,profileNameBox,c);
        javaMethodEDT('add',propertiesPanel,profileNameBox);
        %         c.gridwidth=2;
        %         javaMethodEDT('setConstraints',gb,combineReplicatesBox,c);
        %         javaMethodEDT('add',propertiesPanel,combineReplicatesBox);
        c.gridx=4;
        c.gridwidth=1;
        javaMethodEDT('setConstraints',gb,profileTypeLabel,c);
        javaMethodEDT('add',propertiesPanel,profileTypeLabel);
        c.gridx=5;
        javaMethodEDT('setConstraints',gb,profileTypeBox,c);
        javaMethodEDT('add',propertiesPanel,profileTypeBox);
        
        javaMethodEDT('add',inputPanel,filePanel,java.awt.BorderLayout.NORTH);
        javaMethodEDT('add',inputPanel,previewPanel,java.awt.BorderLayout.CENTER);
        javaMethodEDT('add',inputPanel,propertiesPanel,java.awt.BorderLayout.SOUTH);
        
        %Set up next/back/cancel buttons
        
        buttonsPanel = javaObjectEDT('javax.swing.JPanel');
        cancelButton = javaObjectEDT('javax.swing.JButton','Cancel');
        nextButton = javaObjectEDT('javax.swing.JButton','Next');
        javaMethodEDT('setEnabled',nextButton,false);
        backButton = javaObjectEDT('javax.swing.JButton','Back');
        javaMethodEDT('setEnabled',backButton,false);
        
        cancelButtonH = handle(cancelButton,'CallbackProperties');
        set(cancelButtonH,'ActionPerformedCallback',@cancelButtonEvent);
        nextButtonH = handle(nextButton,'CallbackProperties');
        set(nextButtonH,'ActionPerformedCallback',@nextButtonEvent);
        backButtonH = handle(backButton,'CallbackProperties');
        set(backButtonH,'ActionPerformedCallback',@fileSelectButtonEvent);
        
        
        javaMethodEDT('setLayout',buttonsPanel,javaObjectEDT('java.awt.FlowLayout'));
        javaMethodEDT('add',buttonsPanel,backButton);
        javaMethodEDT('add',buttonsPanel,cancelButton);
        javaMethodEDT('add',buttonsPanel,nextButton);
        
        
        %add to content pane
        javaMethodEDT('add',contentPane,inputPanel,java.awt.BorderLayout.CENTER);
        javaMethodEDT('add',contentPane,buttonsPanel,java.awt.BorderLayout.SOUTH);
        
        %display frame
        javaMethodEDT('pack',frame);
        
        %workaround to get horizontal scroll bar to appear after having empty table
        javaMethodEDT('setModel',previewTable,javaObjectEDT('javax.swing.table.DefaultTableModel'));
        
        javaMethodEDT('setVisible',frame,true);
        
        if dataStruct.hasData
            updateUI;
        end
        
        
        
        
        
        function [y, trend, beta] = detrendLinear(data,timescale,replicates)
            
            t = timescale ./ timescale(end);
            t=t';
            data=data';
            
            %least squares regression
            X  = [repmat(t,1,replicates) ones(length(t) * replicates,1)];
            beta = (X' * X) \ X' * data;
            trend = (t.*beta(1)) + beta(2);
            y = data - trend;
            
        end
        
        function [] = showProfilePreview(varargin)
            if isempty(previewFigH) || ~ishghandle(previewFigH)
                previewFigH = figure;
            end
            figure(previewFigH);
            
            timescale = dataStruct.timescale;
            %get first profile to plot
            gene = dataStruct.rawText(2,dataStruct.nameColumnSelected);
            allNames = dataStruct.rawText(2:end,dataStruct.nameColumnSelected);
            if sum(strcmpi(gene,allNames)) > 1
                %                 disp('combining non-unique rows');
                idx = find(strcmpi(allNames,gene));
                d = dataStruct.rawData(2:end,:);
                data = d(idx,:);
                
                reps = length(idx) * dataStruct.replicates;
                
                d = data';
                data = d(:)';
                numProfiles = length(idx);
            else
                reps = dataStruct.replicates;
                data = dataStruct.rawData(2,:);
                numProfiles = 1;
            end
            detrendName = char(detrendArray(javaMethodEDT('getSelectedIndex',detrendBox)+1));
            set(previewFigH,'Name',['Preview data - ' char(gene) ' (' int2str(numProfiles) ' profiles)  Data detrending: ' detrendName],'NumberTitle','off');
            
            %switch javaMethodEDT('getSelectedIndex',detrendBox)
            switch javaMethodEDT('getSelectedItem',detrendBox)
                case char(detrendLogLinear_str) %log linear
                    newData = nan(1,length(data));
                    trendData = nan(1,length(data));
                    betas = nan(reps,2);
                    timepoints = length(timescale);
                    for z = 1:reps
                        
                        [detrendData, trend, beta] = detrendLinear(log(data( ((z-1)*timepoints+1):(z*timepoints) )),timescale-timescale(1),1);
                        newData(((z-1)*timepoints+1):(z*timepoints)) = exp(detrendData)';
                        trendData(((z-1)*timepoints+1):(z*timepoints)) = exp(trend)';
                        %newData(((z-1)*timepoints+1):(z*timepoints)) = newData(((z-1)*timepoints+1):(z*timepoints)) +  trendData((z-1)*timepoints+1);
                        betas(z,:) = beta;
                    end
                    data = newData;
                case char(detrendLinear_str) %linear
                    newData = nan(1,length(data));
                    trendData = nan(1,length(data));
                    betas = nan(reps,2);
                    timepoints = length(timescale);
                    for z = 1:reps
                        
                        [detrendData, trend, beta] = detrendLinear(data( ((z-1)*timepoints+1):(z*timepoints) ),timescale-timescale(1),1);
                        newData(((z-1)*timepoints+1):(z*timepoints)) = detrendData';
                        trendData(((z-1)*timepoints+1):(z*timepoints)) = trend';
                        %newData(((z-1)*timepoints+1):(z*timepoints)) = newData(((z-1)*timepoints+1):(z*timepoints)) +  trendData((z-1)*timepoints+1);
                        betas(z,:) = beta;
                    end
                    data = newData;
            end
            
            plot(repmat(timescale,1,reps),data,'ob');
        end
        
        function [] = fileSelectButtonEvent(varargin)
            fileChooser = javaObjectEDT('javax.swing.JFileChooser',pwd);
            returnVal = javaMethodEDT('showOpenDialog',fileChooser,frame);
            if returnVal == javax.swing.JFileChooser.APPROVE_OPTION
                currentSelectedFile = javaMethodEDT('getAbsolutePath',javaMethodEDT('getSelectedFile',fileChooser));
                parseFile(char(currentSelectedFile));
            end
        end
        
        function [] = nextButtonEvent(varargin)
            if ishghandle(previewFigH)
                close(previewFigH);
            end
            dataStruct.dataType = profileTypeBox.getSelectedItem;
            dataStruct.dataDetrend = detrendBox.getSelectedItem;
            algorithmSelectUI(frame);
        end
        
        function [] = profileNameSelectEvent(parent,event)
            idx = get(parent,'SelectedIndex')+1;
            dataStruct.nameColumnSelected = idx;
            
            uniqueRowNames = length(unique(dataStruct.rawText(2:end,idx)));
            dataStruct.uniqueRowNames = uniqueRowNames;
            javaMethodEDT('setText',profilesField,java.lang.String([int2str(uniqueRowNames) '/' int2str(size(dataStruct.rawText,1)-1)]));
        end
        
        function [] = profileTypeSelectEvent(parent,event)
            profType =  get(parent,'SelectedItem');
            dataStruct.dataType = profType;
            
        end
        
        function [] = cancelButtonEvent(varargin)
            %         java.system.GC;
            dataStruct = [];
            acceptedInputs = false;
            javaMethodEDT('dispose',frame);
            
            if ishghandle(previewFigH)
                close(previewFigH);
            end
            
            close(hidFig);
        end
        
        %         function [] = combineReplicatesEvent(parent,event)
        %             dataStruct.combineReplicates = combineReplicatesBox.isSelected;
        %             updateProfiles;
        %         end
        
        function [] = updateProfiles
            uniqueRowNames = dataStruct.uniqueRowNames;
            textData = dataStruct.rawText;
            if dataStruct.combineReplicates
                javaMethodEDT('setText',profilesField,java.lang.String([int2str(uniqueRowNames) '/' int2str(uniqueRowNames)]));
            else
                javaMethodEDT('setText',profilesField,java.lang.String([int2str(uniqueRowNames) '/' int2str(size(textData,1)-1)]));
            end
        end
        
        function [] = parseFile(currentFile)
            
            %fix delim to '\t' for now
            delim = '\t';
            try
                d = importdata(currentFile,delim);
                
                %data needs to contain the following fields:
                % data
                % textdata
                % rowheaders
                
                fields = {'data','textdata'};
                hasFields = true;
                for x = 1:length(fields)
                    if ~isfield(d,fields{x})
                        hasFields = false;
                    end
                end
                
                if ~hasFields
                    error('Invalid file contents');
                end
                
                timescale = unique(d.data(1,:));
                replicates = 1;
                for x = 1:length(timescale)
                    if sum(timescale(x) == d.data(1,:)) > replicates
                        replicates = sum(timescale(x) == d.data(1,:));
                    end
                end
                
                timepoints = length(timescale);
                %                 dataStruct = struct;
                dataStruct.file = currentFile;
                dataStruct.rawData = d.data;
                dataStruct.rawText = d.textdata;
                dataStruct.timescale = timescale;
                dataStruct.timepoints = timepoints;
                dataStruct.replicates = replicates;
                dataStruct.nameColumnSelected = 1;
                dataStruct.totalProfileNum = size(dataStruct.rawData,1)-1;
                dataStruct.uniqueRowNames = length(unique(d.textdata(2:end,dataStruct.nameColumnSelected)));
                dataStruct.dataType = profileTypeBox.getSelectedItem;
                dataStruct.dataDetrend = detrendBox.getSelectedItem;
                
                %                 if dataStruct.uniqueRowNames < dataStruct.totalProfileNum
                %                     javaMethodEDT('setSelected',combineReplicatesBox,true);
                %                 else
                %                     javaMethodEDT('setSelected',combineReplicatesBox,false);
                %                 end
                %                 dataStruct.combineReplicates = combineReplicatesBox.isSelected;
                
                selectedFile = currentFile;
                dataStruct.hasData = true;
                
                updateUI;
            catch e
                [path,name,ext] = fileparts(currentFile);
                javaMethodEDT('showMessageDialog',javax.swing.JOptionPane,frame,['"' name ext '" is not a valid data file'],'Invalid data file',javax.swing.JOptionPane.ERROR_MESSAGE);
                
            end
        end
        
        function [] = updateUI
            data = dataStruct.rawData;
            textData = dataStruct.rawText;
            
            columnNames = javaArray('java.lang.String',size(textData,2)+size(data,2));
            for x = 1:size(textData,2)
                columnNames(x) = java.lang.String(char(textData{1,x}));
            end
            for x = 1:size(data,2)
                columnNames(x+size(textData,2)) = java.lang.String(num2str( data(1,x) ));
            end
            rowsToShow = min(size(textData,1)-1,100);
            tableData = javaArray('java.lang.Object',rowsToShow,columnNames.length);
            for x = 1:rowsToShow
                for y = 1:size(textData,2);
                    tableData(x,y) = java.lang.String(textData{x+1,y});
                end
            end
            for x = 1:rowsToShow
                for y = 1:size(data,2)
                    tableData(x,y+size(textData,2)) = java.lang.Double(data(x+1,y));
                end
            end
            tableModel = javaObjectEDT('javax.swing.table.DefaultTableModel',tableData,columnNames);
            javaMethodEDT('setModel',previewTable,tableModel);
            cm = javaMethodEDT('getColumnModel',previewTable);
            for x = 1:columnNames.length
                col = javaMethodEDT('getColumn',cm,x-1);
                javaMethodEDT('setPreferredWidth',col,100);
            end
            
            switch char(dataStruct.dataType)
                case char(mRNA_str)
                    javaMethodEDT('setSelectedIndex',profileTypeBox,0);
                case char(protein_reporter_str)
                    javaMethodEDT('setSelectedIndex',profileTypeBox,1);
                case char(protein_actual_str)
                    javaMethodEDT('setSelectedIndex',profileTypeBox,2);
            end
            timescale = dataStruct.timescale;
            timepoints = dataStruct.timepoints;
            replicates = dataStruct.replicates;
            uniqueRowNames = dataStruct.uniqueRowNames;
            
            javaMethodEDT('setText',fileSelectField,dataStruct.file);
            tsString = num2str(timescale(1));
            for x = 2:length(timescale)
                tsString = [tsString ', ' num2str(timescale(x))];
            end
            javaMethodEDT('setText',timescaleField,java.lang.String(tsString));
            javaMethodEDT('setText',timepointsField,java.lang.String(int2str(timepoints)));
            javaMethodEDT('setText',replicatesField,java.lang.String(int2str(replicates)));
            %             javaMethodEDT('setText',datapointsField,java.lang.String(int2str(size(data,2))));
            %             if dataStruct.combineReplicates
            %                 javaMethodEDT('setText',profilesField,java.lang.String([int2str(uniqueRowNames) '/' int2str(uniqueRowNames)]));
            %             else
            javaMethodEDT('setText',profilesField,java.lang.String([int2str(uniqueRowNames) '/' int2str(size(textData,1)-1)]));
            %             end
            javaMethodEDT('setEnabled',previewButton,true);
            
            javaMethodEDT('removeAllItems',profileNameBox);
            for x = 1:size(textData,2)
                javaMethodEDT('addItem',profileNameBox,java.lang.String(textData{1,x}));
            end
            
            javaMethodEDT('setSelectedIndex',profileNameBox,dataStruct.nameColumnSelected-1);
            javaMethodEDT('setSelectedItem',detrendBox,dataStruct.dataDetrend);
            
            javaMethodEDT('setEnabled',nextButton,true);
            
            %             dataStruct.combineReplicates = combineReplicatesBox.isSelected;
        end
        
    end

    function [] = algorithmSelectUI(frame)
        
        javaMethodEDT('setTitle',frame,'ReTrOS: Select Algorithm Parameters');
        contentPane = javaMethodEDT('getContentPane',frame);
        javaMethodEDT('removeAll',contentPane);
        
        %algorithm general details
        generalParamPanel = javaObject('javax.swing.JPanel',javaObjectEDT('java.awt.BorderLayout'));
        
        outputDir = [pwd filesep 'output'];
        [filePath, fileName, fileExt] = fileparts(dataStruct.file);
        outputLocationButton = javaObjectEDT('javax.swing.JButton','Location:');
        %not allowing changing output directory at the moment
        javaMethodEDT('setEnabled',outputLocationButton,false);
        outputLocationField = javaObjectEDT('javax.swing.JTextField',outputDir,50);
        outputNameLabel = javaObjectEDT('javax.swing.JLabel','Name:');
        outputNameField = javaObjectEDT('javax.swing.JTextField',fileName,50);
        
        gb = javaObjectEDT('java.awt.GridBagLayout');
        p = javaObjectEDT('javax.swing.JPanel',gb);
        outputPanelBorder = javaMethodEDT('createTitledBorder','javax.swing.BorderFactory','File output');
        javaMethodEDT('setBorder',p,outputPanelBorder);
        c = javaObjectEDT('java.awt.GridBagConstraints');
        c.fill = java.awt.GridBagConstraints.BOTH;
        javaMethodEDT('setConstraints',gb,outputLocationButton,c);
        javaMethodEDT('add',p,outputLocationButton);
        c.gridx = 1;
        c.gridwidth=java.awt.GridBagConstraints.REMAINDER;
        javaMethodEDT('setConstraints',gb,outputLocationField,c);
        javaMethodEDT('add',p,outputLocationField);
        c.gridx = 0;
        c.gridy = 1;
        javaMethodEDT('setConstraints',gb,outputNameLabel,c);
        javaMethodEDT('add',p,outputNameLabel);
        c.gridx = 1;
        javaMethodEDT('setConstraints',gb,outputNameField,c);
        javaMethodEDT('add',p,outputNameField);
        
        
        javaMethodEDT('add',generalParamPanel,p,java.awt.BorderLayout.SOUTH);
        
        degRateDistBox = javaObjectEDT('javax.swing.JComboBox',degRateDistArray);
        javaMethodEDT('setSelectedIndex',degRateDistBox,degRateDistIdx);
        degRateDistBoxH = handle(degRateDistBox,'CallbackProperties');
        set(degRateDistBoxH,'ActionPerformedCallback',@updateDegRateLabels);
        
        degRateDefaultBox = javaObjectEDT('javax.swing.JComboBox',degRateDefaultArray);
        if strfind(dataStruct.dataType,'reporter')
            javaMethodEDT('setSelectedIndex',degRateDefaultBox,degRateReporterDefaultIdx);
        else
            javaMethodEDT('setSelectedIndex',degRateDefaultBox,degRateActualDefaultIdx);
        end
        
        degRateDefaultBoxH = handle(degRateDefaultBox,'CallbackProperties');
        set(degRateDefaultBoxH,'ActionPerformedCallback',@updateDegRateValues);
        
        degRateDistLabel = javaObjectEDT('javax.swing.JLabel','Distribution type:');
        degRateDefaultLabel = javaObjectEDT('javax.swing.JLabel','Default values:');
        
        gb = javaObjectEDT('java.awt.GridBagLayout');
        dPanel = javaObjectEDT('javax.swing.JPanel',gb);
        c = javaObjectEDT('java.awt.GridBagConstraints');
        c.fill = java.awt.GridBagConstraints.BOTH;
        javaMethodEDT('setConstraints',gb,degRateDistLabel,c);
        javaMethodEDT('add',dPanel,degRateDistLabel);
        c.gridx = 1;
        c.gridwidth=java.awt.GridBagConstraints.REMAINDER;
        javaMethodEDT('setConstraints',gb,degRateDistBox,c);
        javaMethodEDT('add',dPanel,degRateDistBox);
        c.gridx = 0;
        c.gridy = 1;
        javaMethodEDT('setConstraints',gb,degRateDefaultLabel,c);
        javaMethodEDT('add',dPanel,degRateDefaultLabel);
        c.gridx = 1;
        javaMethodEDT('setConstraints',gb,degRateDefaultBox,c);
        javaMethodEDT('add',dPanel,degRateDefaultBox);
        
        
        mRNADegPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridLayout',1,4));
        %         mRNADegPanelBorder = javaMethodEDT('createTitledBorder','javax.swing.BorderFactory','mRNA (observed)');
        %         javaMethodEDT('setBorder',mRNADegPanel,mRNADegPanelBorder);
        mRNADegLabel1 = javaObjectEDT('javax.swing.JLabel','            ');
        mRNADegField1 = javaObjectEDT('javax.swing.JTextField',num2str(mRNA_degRate1),6);
        mRNADegLabel2 = javaObjectEDT('javax.swing.JLabel','            ');
        mRNADegField2 = javaObjectEDT('javax.swing.JTextField',num2str(mRNA_degRate2),6);
        javaMethodEDT('add',mRNADegPanel,mRNADegLabel1);
        javaMethodEDT('add',mRNADegPanel,mRNADegField1);
        javaMethodEDT('add',mRNADegPanel,mRNADegLabel2);
        javaMethodEDT('add',mRNADegPanel,mRNADegField2);
        
        javaMethodEDT('setEditable',mRNADegField1,false);
        javaMethodEDT('setEditable',mRNADegField2,false);
        
        degPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.FlowLayout'));
        javaMethodEDT('add',degPanel,dPanel);
        javaMethodEDT('add',degPanel,mRNADegPanel);
        
        if strfind(dataStruct.dataType,'protein')
            proteinDegPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridLayout',1,4));
            %             proteinDegPanelBorder = javaMethodEDT('createTitledBorder','javax.swing.BorderFactory','protein (observed)');
            %             javaMethodEDT('setBorder',proteinDegPanel,proteinDegPanelBorder);
            proteinDegLabel1 = javaObjectEDT('javax.swing.JLabel','      ');
            proteinDegField1 = javaObjectEDT('javax.swing.JTextField',num2str(protein_degRate1),6);
            proteinDegLabel2 = javaObjectEDT('javax.swing.JLabel','      ');
            proteinDegField2 = javaObjectEDT('javax.swing.JTextField',num2str(protein_degRate2),6);
            
            javaMethodEDT('add',proteinDegPanel,proteinDegLabel1);
            javaMethodEDT('add',proteinDegPanel,proteinDegField1);
            javaMethodEDT('add',proteinDegPanel,proteinDegLabel2);
            javaMethodEDT('add',proteinDegPanel,proteinDegField2);
            
            javaMethodEDT('setEditable',proteinDegField1,false);
            javaMethodEDT('setEditable',proteinDegField2,false);
            
            degPanel1 = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridLayout',2,1));
            javaMethodEDT('add',degPanel1,mRNADegPanel);
            javaMethodEDT('add',degPanel1,proteinDegPanel);
            
            javaMethodEDT('add',degPanel,degPanel1);
        end
        
        updateDegRateValues;
        
        degFilePanel1 = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.FlowLayout'));
        degFileCheckbox = javaObjectEDT('javax.swing.JCheckBox','Use rates from file?',false);
        degFileChooser = javaObjectEDT('javax.swing.JButton','Select file:');
        degFileChooserH = handle(degFileChooser,'CallbackProperties');
        set(degFileChooserH,'ActionPerformedCallback',@degRateFileChooseEvent);
        degFileField = javaObjectEDT('javax.swing.JTextField','No file selected',30);
        
        javaMethodEDT('setEnabled',degFileField,false);
        javaMethodEDT('add',degFilePanel1,degFileCheckbox);
        javaMethodEDT('add',degFilePanel1,degFileChooser);
        degFilePanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridLayout',2,1));
        javaMethodEDT('add',degFilePanel,degFilePanel1);
        javaMethodEDT('add',degFilePanel,degFileField);
        
        allDegPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.BorderLayout'));
        allDegPanelBorder = javaMethodEDT('createTitledBorder','javax.swing.BorderFactory','Degradation rates');
        javaMethodEDT('setBorder',allDegPanel,allDegPanelBorder);
        javaMethodEDT('add',allDegPanel,degPanel,java.awt.BorderLayout.CENTER);
        javaMethodEDT('add',allDegPanel,degFilePanel,java.awt.BorderLayout.SOUTH);
        
        %algorithm specific details:
        %smooth -
        %bootstrap samples
        %kernel
        %polynomial degree
        %time resolution
        %bandwidth range
        %
        %switch -
        %MCMC iterations
        %samples to plot
        %
        
        
        %general parameters select (these parameters required for both methods)
        paramSelectPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.BorderLayout'));
        paramPanelBorder = javaMethodEDT('createTitledBorder','javax.swing.BorderFactory','Algorithm Parameters');
        javaMethodEDT('setBorder',paramSelectPanel,paramPanelBorder);
        
        %         javaMethodEDT('add',paramSelectPanel,generalParamPanel,java.awt.BorderLayout.NORTH);
        javaMethodEDT('add',paramSelectPanel,allDegPanel,java.awt.BorderLayout.CENTER);
        
        %algorithm specific parameters
        tabbedPane = javaObjectEDT('javax.swing.JTabbedPane');
        
        %smooth options:
        %bootstrap samples
        %kernel type
        %degree polynomial
        %bandwidth range
        %time resolution
        
        smoothPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridLayout',8,2));
        %         paramLabel = javaObjectEDT('javax.swing.JLabel','Parameter');
        %         javaMethodEDT('add',smoothPanel,paramLabel);
        %             javaMethodEDT('getFont',paramLabel)
        %             javaMethodEDT('deriveFont',javaMethodEDT('getFont',paramLabel),java.awt.Font.BOLD)
        %             javaMethodEDT('setFont',paramLabel,javaMethodEDT('deriveFont',javaMethodEDT('getFont',paramLabel),java.awt.Font.BOLD));
        
        %         valueLabel = javaObjectEDT('javax.swing.JLabel','Value');
        %         javaMethodEDT('add',smoothPanel,valueLabel);
        
        
        smoothKernelTypeBox = javaObjectEDT('javax.swing.JComboBox',smoothKernelTypeArray);
        javaMethodEDT('setSelectedIndex',smoothKernelTypeBox,smoothKernelTypeDefaultIdx);
        javaMethodEDT('setToolTipText',smoothKernelTypeBox,'Selects which kernel to use to smooth input data');
        smoothKernelDegreeBox = javaObjectEDT('javax.swing.JComboBox',smoothKernelDegreeArray);
        javaMethodEDT('setSelectedIndex',smoothKernelDegreeBox,smoothKernelDegreeDefaultIdx);
        javaMethodEDT('setToolTipText',smoothKernelDegreeBox,'Selects the polynomial degree used in kernel smoothing regression');
        
        javaMethodEDT('add',smoothPanel,javaObjectEDT('javax.swing.JLabel','bootstrap samples'));
        smoothBootstrapField = javaObjectEDT('javax.swing.JTextField',num2str(smoothBootstrapSamplesDefault));
        javaMethodEDT('setToolTipText',smoothBootstrapField,'Selects the number of bootstrap samples to generate for each profile');
        javaMethodEDT('add',smoothPanel,smoothBootstrapField);
        javaMethodEDT('add',smoothPanel,javaObjectEDT('javax.swing.JLabel','combine profile names as replicates'));
        smoothCombineNamesField = javaObjectEDT('javax.swing.JCheckBox');
        javaMethodEDT('setToolTipText',smoothCombineNamesField,'<html>Selects whether to combine profiles with same name as replicates.<br>If false, then profiles with the same name will have _X appended to name</html>');
        if dataStruct.totalProfileNum > dataStruct.uniqueRowNames
            javaMethodEDT('setSelected',smoothCombineNamesField,true);
        else
            javaMethodEDT('setSelected',smoothCombineNamesField,false);
        end
        javaMethodEDT('add',smoothPanel,smoothCombineNamesField);
        javaMethodEDT('add',smoothPanel,javaObjectEDT('javax.swing.JLabel','kernel type'));
        javaMethodEDT('add',smoothPanel,smoothKernelTypeBox);
        javaMethodEDT('add',smoothPanel,javaObjectEDT('javax.swing.JLabel','kernel polynomial degree'));
        javaMethodEDT('add',smoothPanel,smoothKernelDegreeBox);
        javaMethodEDT('add',smoothPanel,javaObjectEDT('javax.swing.JLabel','Minimum bandwidth'));
        smoothMinBandwidthField = javaObjectEDT('javax.swing.JTextField',num2str(smoothBandwidthMinDefault));
        javaMethodEDT('setToolTipText',smoothMinBandwidthField,'Selects the minimum bandwidth (in hours) used in kernel smoothing');
        javaMethodEDT('add',smoothPanel,smoothMinBandwidthField);
        javaMethodEDT('add',smoothPanel,javaObjectEDT('javax.swing.JLabel','Maximum bandwidth'));
        smoothMaxBandwidthField = javaObjectEDT('javax.swing.JTextField',num2str(smoothBandwidthMaxDefault));
        javaMethodEDT('setToolTipText',smoothMaxBandwidthField,'Selects the maximum bandwidth (in hours) used in kernel smoothing');
        javaMethodEDT('add',smoothPanel,smoothMaxBandwidthField);
        javaMethodEDT('add',smoothPanel,javaObjectEDT('javax.swing.JLabel','Bandwidth increment'));
        smoothIncBandwidthField = javaObjectEDT('javax.swing.JTextField',num2str(smoothBandwidthIncDefault));
        javaMethodEDT('setToolTipText',smoothIncBandwidthField,'Selects the bandwidth increment (in hours) used in finding optimal bandwidth');
        javaMethodEDT('add',smoothPanel,smoothIncBandwidthField);
        javaMethodEDT('add',smoothPanel,javaObjectEDT('javax.swing.JLabel','Time resolution'));
        smoothTimeResField = javaObjectEDT('javax.swing.JTextField',num2str(smoothTimeResolutionDefault));
        javaMethodEDT('setToolTipText',smoothTimeResField,'Selects the time increment (in hours) used in the transcription back-calculation');
        javaMethodEDT('add',smoothPanel,smoothTimeResField);
        
        
        %MCMC iterations
        %use initial switch guess
        %switch strength
        %max switches
        %expected number of switches
        %switch delay
        %baseline factor
        %burn-in
        %number samples to plot
        %number samples to use for posteriors (thinning)
        
        switchPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridLayout',13,2));
        
        %         paramLabel = javaObjectEDT('javax.swing.JLabel','Parameter');
        %         javaMethodEDT('add',switchPanel,paramLabel);
        %         valueLabel = javaObjectEDT('javax.swing.JLabel','Value');
        %         javaMethodEDT('add',switchPanel,valueLabel);
        
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','MCMC iterations'));
        switchIterationsField = javaObjectEDT('javax.swing.JTextField',num2str(switchIterationsDefault));
        javaMethodEDT('setToolTipText',switchIterationsField,'Selects the number of MCMC iterations to run');
        javaMethodEDT('add',switchPanel,switchIterationsField);
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','combine profile names as replicates'));
        switchCombineNamesField = javaObjectEDT('javax.swing.JCheckBox');
        javaMethodEDT('setToolTipText',switchCombineNamesField,'<html>Selects whether to combine profiles with same name as replicates.<br>If false, then profiles with the same name will have _X appended to name</html>');
        if dataStruct.totalProfileNum > dataStruct.uniqueRowNames
            javaMethodEDT('setSelected',switchCombineNamesField,true);
        else
            javaMethodEDT('setSelected',switchCombineNamesField,false);
        end
        javaMethodEDT('add',switchPanel,switchCombineNamesField);
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','Regression method'));
        switchRegressionMethodBox = javaObjectEDT('javax.swing.JComboBox',switchRegressionMethodArray);
        javaMethodEDT('add',switchPanel,switchRegressionMethodBox);
        javaMethodEDT('setSelectedIndex',switchRegressionMethodBox,switchRegressionMethodDefaultIdx);
        javaMethodEDT('setToolTipText',switchRegressionMethodBox,'Selects the regression method used');
        
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','guess initial switch points'));
        switchGuessSwitchBox = javaObjectEDT('javax.swing.JCheckBox');
        javaMethodEDT('setToolTipText',switchGuessSwitchBox,'<html>Selects whether to estimate initial switch points from profile<br>If true, then switch points are estimated based on changes in gradient of profile</html>');
        javaMethodEDT('setSelected',switchGuessSwitchBox,true);
        javaMethodEDT('add',switchPanel,switchGuessSwitchBox);
        
        %estimate some switch parameters from data
        switchDelay = min(dataStruct.timescale(2:end) - dataStruct.timescale(1:(end-1)));
        maxSwitchesDelay = floor((dataStruct.timescale(end)-dataStruct.timescale(1)- 2*switchDelay) / switchDelay);
        maxSwitchesDelay = min(maxSwitchesDelay,length(dataStruct.timescale)-2);
        
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','maximum number of switch points'));
        switchMaxSwitchField = javaObjectEDT('javax.swing.JTextField',num2str(maxSwitchesDelay));
        javaMethodEDT('setToolTipText',switchMaxSwitchField,'Selects the maximum number of switch points in model');
        javaMethodEDT('add',switchPanel,switchMaxSwitchField);
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','expected number of switch points'));
        %switchExpectedSwitchField = javaObjectEDT('javax.swing.JTextField',num2str(ceil(maxSwitchesDelay/2)));
        switchExpectedSwitchField = javaObjectEDT('javax.swing.JTextField',num2str(ceil(dataStruct.timepoints/3)));
        javaMethodEDT('setToolTipText',switchExpectedSwitchField,'Selects the expected number of switch points in model');
        javaMethodEDT('add',switchPanel,switchExpectedSwitchField);
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','minimum time between switch points'));
        switchMinSwTimeField = javaObjectEDT('javax.swing.JTextField',num2str(switchDelay));
        javaMethodEDT('setToolTipText',switchMinSwTimeField,'Selects the minimum time (in hours) between switch points in model');
        javaMethodEDT('add',switchPanel,switchMinSwTimeField);
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','burn-in period'));
        switchBurnInField = javaObjectEDT('javax.swing.JTextField',num2str(switchBurnInDefault));
        javaMethodEDT('setToolTipText',switchBurnInField,'Selects the proportion of MCMC iterations to remove as "burn-in" from beginning of chains');
        javaMethodEDT('add',switchPanel,switchBurnInField);
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','number of samples to estimate posteriors'));
        switchSampleEstNumField = javaObjectEDT('javax.swing.JTextField',num2str(switchPosteriorSampleNumDefault));
        javaMethodEDT('setToolTipText',switchSampleEstNumField,'Selects the maximum number of MCMC samples to estimate posterior distributions from');
        javaMethodEDT('add',switchPanel,switchSampleEstNumField);
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','switch density baseline factor'));
        switchBaselineFactorField = javaObjectEDT('javax.swing.JTextField',num2str(switchBaselineFactorDefault));
        javaMethodEDT('setToolTipText',switchBaselineFactorField,'Selects the baseline scaling factor for switch point density function');
        javaMethodEDT('add',switchPanel,switchBaselineFactorField);
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','minimum switch strength'));
        switchMinSwStrField = javaObjectEDT('javax.swing.JTextField',num2str(switchMinimumSwStrDefault));
        javaMethodEDT('setToolTipText',switchMinSwStrField,'<html>Selects the minimum "strength" of each posterior switch distribution.<br>If a switch distribution is below this value it is removed from the posterior model</html>');
        javaMethodEDT('add',switchPanel,switchMinSwStrField);
        
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','maximum switch deviation factor'));
        switchMaxSwDevField = javaObjectEDT('javax.swing.JTextField',num2str(switchMaxSwDeviationFactorDefault));
        javaMethodEDT('setToolTipText',switchMaxSwDevField,'<html>Selects the maximum deviation of each posterior switch distribution as a factor of time point sampling.<br>If a switch deviation is above this factor times time point sampling it is removed from the posterior model</html>');
        javaMethodEDT('add',switchPanel,switchMaxSwDevField);
        
        javaMethodEDT('add',switchPanel,javaObjectEDT('javax.swing.JLabel','number of samples to plot'));
        switchSamplePlotNumField = javaObjectEDT('javax.swing.JTextField',num2str(switchPlotSampleNumDefault));
        javaMethodEDT('setToolTipText',switchSamplePlotNumField,'Selects the maximum number of MCMC samples to plot');
        javaMethodEDT('add',switchPanel,switchSamplePlotNumField);
        
        smoothScroll = javaObjectEDT('javax.swing.JScrollPane',smoothPanel);
        switchScroll = javaObjectEDT('javax.swing.JScrollPane',switchPanel);
        javaMethodEDT('setPreferredSize',smoothScroll,javaObjectEDT('java.awt.Dimension',200,215));
        javaMethodEDT('setPreferredSize',switchScroll,javaObjectEDT('java.awt.Dimension',200,215));
        
        javaMethodEDT('addTab',tabbedPane,'Smooth options',smoothScroll);
        javaMethodEDT('addTab',tabbedPane,'Switch options',switchScroll);
        
        if dataStruct.hasAlgorithmParameters
            updateParametersUI;
        end
        
        javaMethodEDT('add',paramSelectPanel,tabbedPane,java.awt.BorderLayout.SOUTH);
        
        %buttons (next/back/cancel)
        buttonsPanel = javaObjectEDT('javax.swing.JPanel');
        cancelButton = javaObjectEDT('javax.swing.JButton','Cancel');
        nextButton = javaObjectEDT('javax.swing.JButton','Next');
        % javaMethodEDT('setEnabled',nextButton,false);
        backButton = javaObjectEDT('javax.swing.JButton','Back');
        
        cancelButtonH = handle(cancelButton,'CallbackProperties');
        set(cancelButtonH,'ActionPerformedCallback',@cancelButtonEvent);
        nextButtonH = handle(nextButton,'CallbackProperties');
        set(nextButtonH,'ActionPerformedCallback',@nextButtonEvent);
        backButtonH = handle(backButton,'CallbackProperties');
        set(backButtonH,'ActionPerformedCallback',@backButtonEvent);
        
        javaMethodEDT('setLayout',buttonsPanel,javaObjectEDT('java.awt.FlowLayout'));
        javaMethodEDT('add',buttonsPanel,backButton);
        javaMethodEDT('add',buttonsPanel,cancelButton);
        javaMethodEDT('add',buttonsPanel,nextButton);
        
        %add to content pane
        javaMethodEDT('add',contentPane,generalParamPanel,java.awt.BorderLayout.CENTER);
        javaMethodEDT('add',contentPane,paramSelectPanel,java.awt.BorderLayout.NORTH);
        javaMethodEDT('add',contentPane,buttonsPanel,java.awt.BorderLayout.SOUTH);
        
        %display frame
        javaMethodEDT('pack',frame);
        
        
        
        
        
        
        
        
        
        %update functions
        function [] = updateDegRateValues(varargin)
            switch javaMethodEDT('getSelectedIndex',degRateDefaultBox)
                case 0 % USER
                    javaMethodEDT('setEnabled',degRateDistBox,true);
                    javaMethodEDT('setEditable',mRNADegField1,true);
                    javaMethodEDT('setEditable',mRNADegField2,true);
                    if strfind(dataStruct.dataType,'protein')
                        javaMethodEDT('setEditable',proteinDegField1,true);
                        javaMethodEDT('setEditable',proteinDegField2,true);
                    end
                case 1 %Arabidposis values
                    javaMethodEDT('setSelectedIndex',degRateDistBox,1);
                    javaMethodEDT('setEnabled',degRateDistBox,false);
                    mRNA_degRate1 = ARAB_MRNA_ALPHA;
                    mRNA_degRate2 = ARAB_MRNA_BETA;
                    javaMethodEDT('setEditable',mRNADegField1,false);
                    javaMethodEDT('setEditable',mRNADegField2,false);
                    
                    if strfind(dataStruct.dataType,'protein')
                        protein_degRate1 = ARAB_PROTEIN_ALPHA;
                        protein_degRate2 = ARAB_PROTEIN_BETA;
                        javaMethodEDT('setEditable',proteinDegField1,false);
                        javaMethodEDT('setEditable',proteinDegField2,false);
                    end
                    
                case 2 %luc reporter
                    javaMethodEDT('setSelectedIndex',degRateDistBox,0);
                    javaMethodEDT('setEnabled',degRateDistBox,false);
                    mRNA_degRate1 = LUC_MRNA_MEAN;
                    mRNA_degRate2 = LUC_MRNA_SDEV;
                    
                    javaMethodEDT('setEditable',mRNADegField1,false);
                    javaMethodEDT('setEditable',mRNADegField2,false);
                    
                    if strfind(dataStruct.dataType,'protein')
                        protein_degRate1 = LUC_PROTEIN_MEAN;
                        protein_degRate2 = LUC_PROTEIN_SDEV;
                        javaMethodEDT('setEditable',proteinDegField1,false);
                        javaMethodEDT('setEditable',proteinDegField2,false);
                    end
            end
            updateDegRateLabels;
        end
        
        function [] = updateDegRateLabels(varargin)
            if javaMethodEDT('getSelectedIndex',degRateDistBox) == 0 %Normal
                javaMethodEDT('setText',mRNADegLabel1,'mRNA mean:');
                javaMethodEDT('setText',mRNADegLabel2,'mRNA std dev:');
                if strfind(dataStruct.dataType,'protein')
                    javaMethodEDT('setText',proteinDegLabel1,'protein mean:');
                    javaMethodEDT('setText',proteinDegLabel2,'protein std dev:');
                end
            elseif javaMethodEDT('getSelectedIndex',degRateDistBox) == 1 %Gamma
                javaMethodEDT('setText',mRNADegLabel1,'mRNA alpha:');
                javaMethodEDT('setText',mRNADegLabel2,'mRNA beta:');
                if strfind(dataStruct.dataType,'protein')
                    javaMethodEDT('setText',proteinDegLabel1,'protein alpha:');
                    javaMethodEDT('setText',proteinDegLabel2,'protein beta:');
                end
            end
            
            javaMethodEDT('setText',mRNADegField1,num2str(mRNA_degRate1));
            javaMethodEDT('setText',mRNADegField2,num2str(mRNA_degRate2));
            if strfind(dataStruct.dataType,'protein')
                javaMethodEDT('setText',proteinDegField1,num2str(protein_degRate1));
                javaMethodEDT('setText',proteinDegField2,num2str(protein_degRate2));
            end
        end
        
        
        %other update function
        function [] = updateParametersUI
            %degradation rates
            javaMethodEDT('setSelectedItem',degRateDistBox,dataStruct.degDistTypeSelected);
            javaMethodEDT('setSelectedItem',degRateDefaultBox,dataStruct.degDefaultSelected);
            mRNA_degRate1 = dataStruct.degRates(1);
            mRNA_degRate2 = dataStruct.degRates(2);
            protein_degRate1 = dataStruct.degRates(3);
            protein_degRate2 = dataStruct.degRates(4);
            updateDegRateValues
            javaMethodEDT('setSelected',degFileCheckbox,dataStruct.useDegRateFile);
            if dataStruct.useDegRateFile
                javaMethodEDT('setText',degFileField,dataStruct.degRateFile);
            end
            %
            %smooth options:
            %bootstrap samples
            javaMethodEDT('setText',smoothBootstrapField,num2str(dataStruct.smoothBootstrapSamples));
            %combine profiles as replicates
            javaMethodEDT('setSelected',smoothCombineNamesField,dataStruct.smoothCombineReplicates);
            %kernel type
            javaMethodEDT('setSelectedItem',smoothKernelTypeBox,smoothKernelTypeArray(dataStruct.smoothKernelFull+1));
            %poly deg
            javaMethodEDT('setSelectedItem',smoothKernelDegreeBox,java.lang.String(num2str(dataStruct.smoothPolyDeg)));
            %min bandwidth
            javaMethodEDT('setText',smoothMinBandwidthField,num2str(dataStruct.smoothBandwidthRange(1)));
            %max bandwidth
            javaMethodEDT('setText',smoothMaxBandwidthField,num2str(dataStruct.smoothBandwidthRange(2)));
            %inc bandwidth
            javaMethodEDT('setText',smoothIncBandwidthField,num2str(dataStruct.smoothBandwidthRange(3)));
            %time res
            javaMethodEDT('setText',smoothTimeResField,num2str(dataStruct.smoothTimeResolution));
            %
            %switch options:
            %iterations
            javaMethodEDT('setText',switchIterationsField,num2str(dataStruct.switchIterations));
            %combine profiles as replicates
            javaMethodEDT('setSelected',switchCombineNamesField,dataStruct.switchCombineReplicates);
            %regression method
            javaMethodEDT('setSelectedItem',switchRegressionMethodBox,dataStruct.switchRegressionMethodString);
            %guess initial switch locations
            javaMethodEDT('setSelected',switchGuessSwitchBox,dataStruct.switchIntialSwitchGuess);
            %max switch number
            javaMethodEDT('setText',switchMaxSwitchField,num2str(dataStruct.switchMaxSwitchNum));
            %expected switch number
            javaMethodEDT('setText',switchExpectedSwitchField,num2str(dataStruct.switchExpSwitchNum));
            %min time between switch
            javaMethodEDT('setText',switchMinSwTimeField,num2str(dataStruct.switchDelayTime));
            %burn in period
            javaMethodEDT('setText',switchBurnInField,num2str(dataStruct.switchBurnIn));
            %samples to use
            javaMethodEDT('setText',switchSampleEstNumField,num2str(dataStruct.switchNumSamplesEstPost));
            %baseline function
            javaMethodEDT('setText',switchBaselineFactorField,num2str(dataStruct.switchBaselineFactor));
            %min sw strength
            javaMethodEDT('setText',switchMinSwStrField,num2str(dataStruct.switchMinSwStr));
            %max switch deviation factor
            javaMethodEDT('setText',switchMaxSwDevField,num2str(dataStruct.switchMaxSwDev));
            %samples to plot
            javaMethodEDT('setText',switchSamplePlotNumField,num2str(dataStruct.switchNumSamplePlot));
            %
            %output dir
            %NOT USED AT THE MOMENT
            %outputLocationField = javaObjectEDT('javax.swing.JTextField',outputDir,50);
            %output name
            javaMethodEDT('setText',outputNameField,dataStruct.outputName);
        end
        
        
        function [] = updateDataStructParams
            %general params
            name = char(javaMethodEDT('getText',outputNameField));
            dataStruct.outputName = checkValidFileChar(name);
            javaMethodEDT('setText',outputNameField,dataStruct.outputName);
            smoothScriptFileName = [pwd filesep 'scripts' filesep 'run_' name '_smooth.m'];
            dataStruct.smoothScriptFile = smoothScriptFileName;
            switchScriptFileName = [pwd filesep 'scripts' filesep 'run_' name '_switch.m'];
            dataStruct.switchScriptFile = switchScriptFileName;
            
            degDistType = char(javaMethodEDT('getSelectedItem',degRateDistBox));
            dataStruct.degDistTypeSelected = degDistType;
            dataStruct.degDefaultSelected = char(javaMethodEDT('getSelectedItem',degRateDefaultBox));
            mRNA1str = char(mRNADegField1.getText);
            mRNA2str = char(mRNADegField2.getText);
            if strfind(dataStruct.dataType,'protein')
                protein1str = char(proteinDegField1.getText);
                protein2str = char(proteinDegField2.getText);
            end
            
            if javaMethodEDT('getSelectedIndex',degRateDefaultBox) > 0
                %want to use a default set
                degRates = [str2double(mRNA1str) str2double(mRNA2str) 0 0];
                %                 fprintf(scriptFID,['degRate_mRNA = [' char(mRNADegField1.getText) ' ' char(mRNADegField2.getText) '];\n']);
                if strfind(dataStruct.dataType,'protein')
                    degRates(3) = str2double(protein1str);
                    degRates(4) = str2double(protein2str);
                    %                     fprintf(scriptFID,['degRate_protein = [' char(proteinDegField1.getText) ' ' char(proteinDegField2.getText) '];\n']);
                end
            else
                num_degRate_mRNA_1 = str2double(mRNA1str);
                num_degRate_mRNA_2 = str2double(mRNA2str);
                if isempty(num_degRate_mRNA_1) || isnan(num_degRate_mRNA_1) || ~isscalar(num_degRate_mRNA_1) || num_degRate_mRNA_1 <= 0 || isempty(num_degRate_mRNA_2) || isnan(num_degRate_mRNA_2) || ~isscalar(num_degRate_mRNA_2) || num_degRate_mRNA_2 <= 0
                    %we have a missing/incorrect degradation rate... use of the Arabidopsis defaults
                    if strfind(dataStruct.dataType,'protein')
                        %luc defaults
                        disp('ReTros: Invalid mRNA degradation input.  Using luciferase reporter defaults');
                        if strcmpi(degDistType,'gamma')
                            degRateMRNA_1 = (LUC_MRNA_MEAN * LUC_MRNA_MEAN) / (LUC_MRNA_SDEV * LUC_MRNA_SDEV);
                            degRateMRNA_2 = (LUC_MRNA_SDEV * LUC_MRNA_SDEV) / LUC_MRNA_MEAN;
                        else
                            degRateMRNA_1 = LUC_MRNA_MEAN;
                            degRateMRNA_2 = LUC_MRNA_SDEV;
                        end
                    else
                        disp('ReTros: Invalid mRNA degradation input.  Using Arabidopsis defaults');
                        if strcmpi(degDistType,'normal')
                            degRateMRNA_1 = ARAB_MRNA_ALPHA * ARAB_MRNA_BETA;
                            degRateMRNA_2 = sqrt(ARAB_MRNA_ALPHA * ARAB_MRNA_BETA * ARAB_MRNA_BETA);
                        else
                            degRateMRNA_1 = ARAB_MRNA_ALPHA;
                            degRateMRNA_2 = ARAB_MRNA_BETA;
                        end
                    end
                else
                    degRateMRNA_1 = num_degRate_mRNA_1;
                    degRateMRNA_2 = num_degRate_mRNA_2;
                end
                degRates = [degRateMRNA_1 degRateMRNA_2 0 0];
                if strfind(dataStruct.dataType,'protein')
                    num_degRate_protein_1 = str2double(protein1str);
                    num_degRate_protein_2 = str2double(protein2str);
                    if isempty(num_degRate_protein_1) || isnan(num_degRate_protein_1) || ~isscalar(num_degRate_protein_1) || num_degRate_protein_1 <= 0 || isempty(num_degRate_protein_2) || isnan(num_degRate_protein_2) || ~isscalar(num_degRate_protein_2) || num_degRate_protein_2 <= 0
                        %we have a missing/incorrect degradation rate... use of the Arabidopsis defaults
                        disp('ReTros: Invalid protein degradation input.  Using luciferase defaults');
                        if strcmpi(degDistType,'gamma')
                            degRateProtein_1 = (LUC_PROTEIN_MEAN * LUC_PROTEIN_MEAN) / (LUC_PROTEIN_SDEV * LUC_PROTEIN_SDEV);
                            degRateProtein_2 = (LUC_PROTEIN_SDEV * LUC_PROTEIN_SDEV) / LUC_PROTEIN_SDEV;
                        else
                            degRateProtein_1 = ARAB_PROTEIN_ALPHA;
                            degRateProtein_2 = ARAB_PROTEIN_BETA;
                        end
                    else
                        degRateProtein_1 = num_degRate_protein_1;
                        degRateProtein_2 = num_degRate_protein_2;
                    end
                    degRates(3) = degRateProtein_1;
                    degRates(4) = degRateProtein_2;
                end
            end
            dataStruct.degRates = degRates;
            dataStruct.useDegRateFile = degFileCheckbox.isSelected;
            if dataStruct.useDegRateFile
                dataStruct.degRateFile = char(degFileField.getText);
            else
                dataStruct.degRateFile = '';
            end
            
            %smooth specific
            if strcmpi(degDistType,'gamma')
                newDegRates(1) = degRates(1) * degRates(2);
                newDegRates(2) = sqrt(degRates(1) * degRates(2) * degRates(2));
                newDegRates(3) = 0;
                newDegRates(4) = 0;
                if strfind(dataStruct.dataType,'protein')
                    newDegRates(3) = degRates(3) * degRates(4);
                    newDegRates(4) = sqrt(degRates(3) * degRates(4) * degRates(4));
                end
                dataStruct.smoothDegRates = newDegRates;
            else
                dataStruct.smoothDegRates = degRates;
            end
            
            bootstrapSamples = str2double(char(smoothBootstrapField.getText));
            if isempty(bootstrapSamples) || isnan(bootstrapSamples) || bootstrapSamples < 1 || ~isscalar(bootstrapSamples) || (round(bootstrapSamples) ~= bootstrapSamples)
                disp(['ReTrOS-smooth: invalid value for ''bootstrapSamples''. Using default of ''' int2str(smoothBootstrapSamplesDefault) '''']);
                bootstrapSamples = smoothBootstrapSamplesDefault;
            end
            if bootstrapSamples > 100000
                disp('ReTrOS-smooth: warning value for ''bootstrapSamples'' is very large, consider reducing value');
            end
            dataStruct.smoothBootstrapSamples = bootstrapSamples;
            
            dataStruct.smoothCombineReplicates = smoothCombineNamesField.isSelected;
            dataStruct.smoothKernelFull = smoothKernelTypeBox.getSelectedIndex;
            switch smoothKernelTypeBox.getSelectedIndex
                case 0
                    kernel = 'ga';
                case 1
                    kernel = 'ep';
                case 2
                    kernel = 'tr';
            end
            dataStruct.smoothKernel = kernel;
            dataStruct.smoothPolyDeg = smoothKernelDegreeBox.getSelectedIndex + 1;
            %bandwidth range
            minBandStr = char(smoothMinBandwidthField.getText);
            minBand = str2double(minBandStr);
            if isempty(minBand) || isnan(minBand) || minBand <= 0 || ~isscalar(minBand)
                disp(['ReTrOS-smooth: invalid minimum bandwidth. Using default of ''' num2str(smoothBandwidthMinDefault) '''']);
                minBand = smoothBandwidthMinDefault;
            end
            maxBandStr = char(smoothMaxBandwidthField.getText);
            maxBand = str2double(maxBandStr);
            if isempty(maxBand) || isnan(maxBand) || maxBand <= 0 || ~isscalar(maxBand)
                disp(['ReTrOS-smooth: invalid maximum bandwidth. Using default of ''' num2str(smoothBandwidthMaxDefault) '''']);
                maxBand = smoothBandwidthMaxDefault;
            end
            incBandStr = char(smoothIncBandwidthField.getText);
            incBand = str2double(incBandStr);
            if isempty(incBand) || isnan(incBand) || incBand <= 0 || ~isscalar(incBand)
                disp(['ReTrOS-smooth: invalid bandwidth increment. Using default of ''' num2str(smoothBandwidthIncDefault) '''']);
                incBand = smoothBandwidthIncDefault;
            end
            if maxBand < minBand
                disp('ReTrOS-smooth: minimum bandwidth larger than maximum. Flipping values');
                temp = minBand;
                minBand = maxBand;
                maxBand = temp;
            end
            if incBand > ((maxBand - minBand)/2)
                disp(['ReTrOS-smooth: bandwidth increment too large. Using ''' num2str((maxBand - minBand) / 2) '''']);
                incBand = (maxBand - minBand) / 2;
            end
            dataStruct.smoothBandwidthRange = [minBand maxBand incBand];
            
            %time resolution
            timeResStr = char(smoothTimeResField.getText);
            timeRes = str2double(timeResStr);
            if isempty(timeRes) || isnan(timeRes) || timeRes < 0.01 || ~isscalar(timeRes)
                disp(['ReTrOS-smooth: invalid value for ''timeResolution''. Using default of ' num2str(smoothTimeResolutionDefault)]);
                timeRes = smoothTimeResolutionDefault;
            end
            dataStruct.smoothTimeResolution = timeRes;
            
            
            
            %switch specific
            if strcmpi(degDistType,'normal')
                mRNADegRates(1) = (degRates(1) * degRates(1)) / (degRates(2) * degRates(2));
                mRNADegRates(2) = (degRates(2) * degRates(2)) / degRates(1);
                if strfind(dataStruct.dataType,'protein')
                    proteinDegRates(1) = (degRates(3) * degRates(3)) / (degRates(4) * degRates(4));
                    proteinDegRates(2) = (degRates(4) * degRates(4)) / degRates(3);
                end
                dataStruct.switchmRNADegRates = mRNADegRates;
                dataStruct.switchProteinDegRates = proteinDegRates;
            else
                dataStruct.switchmRNADegRates = degRates(1:2);
                dataStruct.switchProteinDegRates = degRates(3:4);
            end
            dataStruct.switchCombineReplicates = switchCombineNamesField.isSelected;
            iterationsStr = char(switchIterationsField.getText);
            iterations = str2double(iterationsStr);
            if isempty(iterations) || isnan(iterations) || iterations < 1000 || ~isscalar(iterations) || (round(iterations) ~= iterations)
                disp(['ReTrOS-switch: invalid value for ''MCMCiterations''. Using default of ''' int2str(switchIterationsDefault) '''']);
                iterations = switchIterationsDefault;
            end
            if iterations > 10000000
                disp('[ReTrOS-switch: warning value for ''MCMCiterations'' is very large, consider reducing value');
            end
            dataStruct.switchIterations = iterations;
            regressionMethod = switchRegressionMethodBox.getSelectedItem;
            dataStruct.switchRegressionMethodString = regressionMethod;
            dataStruct.switchRegressionMethod = strrep(char(regressionMethod),' ','');
            %guess initial switches
            dataStruct.switchIntialSwitchGuess = switchGuessSwitchBox.isSelected;
            %switch delay time
            switchDelayStr = char(switchMinSwTimeField.getText);
            switchDelay = str2double(switchDelayStr);
            if isempty(switchDelay) || isnan(switchDelay) || switchDelay < 0 || ~isscalar(switchDelay) || (switchDelay > (dataStruct.timescale(end)/2))
                switchDelay = min(dataStruct.timescale(2:end) - dataStruct.timescale(1:(end-1)));
                disp(['ReTrOS-switch: invalid value for ''switchDelayTime''. Using default of ''' num2str(switchDelay) '''']);
            end
            dataStruct.switchDelayTime = switchDelay;
            %maximum switch number
            maxSwitchNumStr = char(switchMaxSwitchField.getText);
            maxSwitchNum = str2double(maxSwitchNumStr);
            if isempty(maxSwitchNum) || isnan(maxSwitchNum) || maxSwitchNum < 0 || ~isscalar(maxSwitchNum) || (round(maxSwitchNum) ~= maxSwitchNum)
                switchDelay = min(dataStruct.timescale(2:end) - dataStruct.timescale(1:(end-1)));
                maxSwitchesDelay = floor((dataStruct.timescale(end)-dataStruct.timescale(1)- 2*dataStruct.switchDelayTime) / dataStruct.switchDelayTime);
                maxSwitchesDelay = min(maxSwitchesDelay,length(dataStruct.timescale)-2);
                disp(['ReTrOS-switch: invalid value for ''maximumSwitchNum''. Using default of ''' int2str(maxSwitchesDelay) '''']);
                maxSwitchNum = maxSwitchesDelay;
            end
            dataStruct.switchMaxSwitchNum = maxSwitchNum;
            %expected switch number
            expSwitchNumStr = char(switchExpectedSwitchField.getText);
            expSwitchNum = str2double(expSwitchNumStr);
            if isempty(expSwitchNum) || isnan(expSwitchNum) || expSwitchNum < 0 || ~isscalar(expSwitchNum) || (round(expSwitchNum) ~= expSwitchNum)
                expSwitchNum = ceil(maxSwitchNum/2);
                disp(['ReTrOS-switch: invalid value for ''expectedSwitchNum''. Using default of ''' int2str(expSwitchNum) '''']);
            end
            dataStruct.switchExpSwitchNum = expSwitchNum;
            %burn-in period
            burnInStr = char(switchBurnInField.getText);
            burnIn = str2double(burnInStr);
            if isempty(burnIn) || isnan(burnIn) || burnIn < 0 || ~isscalar(burnIn) || burnIn > 0.99
                burnIn = switchBurnInDefault;
                disp(['ReTrOS-switch: invalid value for ''switchDelayTime''. Using default of ''' num2str(burnIn) '''']);
            end
            dataStruct.switchBurnIn = burnIn;
            %num samples for posteriors
            numSamplesEstStr = char(switchSampleEstNumField.getText);
            numSamplesEst = str2double(numSamplesEstStr);
            if isempty(numSamplesEst) || isnan(numSamplesEst) || numSamplesEst < 1000 || ~isscalar(numSamplesEst) || (round(numSamplesEst) ~= numSamplesEst)
                numSamplesEst = switchPosteriorSampleNumDefault;
                disp(['ReTrOS-switch: invalid value for ''numSamplesEstPost''. Using default of ''' num2str(numSamplesEst) '''']);
            end
            dataStruct.switchNumSamplesEstPost = numSamplesEst;
            %switch baseline factor
            baselineStr = char(switchBaselineFactorField.getText);
            baseline = str2double(baselineStr);
            if isempty(baseline) || isnan(baseline) || baseline < 0 || ~isscalar(baseline)
                baseline = switchBaselineFactorDefault;
                disp(['ReTrOS-switch: invalid value for ''switchBaselineFactor''. Using default of ''' num2str(baseline) '''']);
            end
            dataStruct.switchBaselineFactor = baseline;
            %minimum switch strength
            minSwStrStr = char(switchMinSwStrField.getText);
            minSwStr = str2double(minSwStrStr);
            if isempty(minSwStr) || isnan(minSwStr) || minSwStr < 0 || ~isscalar(minSwStr) || minSwStr > 0.95
                minSwStr = switchMinimumSwStrDefault;
                disp(['ReTrOS-switch: invalid value for ''switchMinimumSwStr''. Using default of ''' num2str(minSwStr) '''']);
            end
            dataStruct.switchMinSwStr = minSwStr;
            %maximum switch deviation factor
            maxSwDevStr = char(switchMaxSwDevField.getText);
            maxSwDevFactor = str2double(maxSwDevStr);
            if isempty(maxSwDevFactor) || isnan(maxSwDevFactor) || maxSwDevFactor < 0 || ~isscalar(maxSwDevFactor)
                maxSwDevFactor = switchMaxSwDeviationFactorDefault;
                disp(['ReTrOS-switch: invalid value for ''maxSwitchDeviationFactor''. Using default of ''' num2str(maxSwDevFactor) '''']);
            end
            dataStruct.switchMaxSwDev = maxSwDevFactor;
            %num samples to plot
            numSamplesPlotStr = char(switchSamplePlotNumField.getText);
            numSamplesPlot = str2double(numSamplesPlotStr);
            if isempty(numSamplesPlot) || isnan(numSamplesPlot) || numSamplesPlot < 1000 || ~isscalar(numSamplesPlot) || (round(numSamplesPlot) ~= numSamplesPlot)
                numSamplesPlot = switchPlotSampleNumDefault;
                disp(['ReTrOS-switch: invalid value for ''numSamplesPlot''. Using default of ''' num2str(numSamplesPlot) '''']);
            end
            dataStruct.switchNumSamplePlot = numSamplesPlot;
            
            
            
            
            dataStruct.hasAlgorithmParameters = true;
        end
        
        %EVENTS
        
        function [] = degRateFileChooseEvent(varargin)
            degRateFileChooser = javaObjectEDT('javax.swing.JFileChooser',pwd);
            returnVal = javaMethodEDT('showOpenDialog',degRateFileChooser,frame);
            if returnVal == javax.swing.JFileChooser.APPROVE_OPTION
                currentSelectedFile = javaMethodEDT('getAbsolutePath',javaMethodEDT('getSelectedFile',degRateFileChooser));
                degFileField.setText(currentSelectedFile);
            end
        end
        
        function [] = backButtonEvent(varargin)
            dataSelectUI(frame);
        end
        
        function [] = nextButtonEvent(varargin)
            %             javaMethodEDT('dispose',frame);
            %             acceptedInputs = true;
            %             disp('goodbye');
            %             close(hidFig);
            updateDataStructParams;
            
            
            runAlgorithmUI(frame);
            
            %             % TEMP FOR TESTING
            %             %smooth
            %             name = char(javaMethodEDT('getText',outputNameField));
            %             scriptFileName = [pwd filesep 'scripts' filesep 'run_' checkValidFileChar(name) '_smooth.m'];
            %             clear(scriptFileName);
            %             run(scriptFileName);
            %             %switch
            %             name = char(javaMethodEDT('getText',outputNameField));
            %             scriptFileName = [pwd filesep 'scripts' filesep 'run_' checkValidFileChar(name) '_switch.m'];
            %             clear(scriptFileName);
            %             run(scriptFileName);
            %
        end
        
        function [] = cancelButtonEvent(varargin)
            %         java.system.GC;
            dataStruct = [];
            acceptedInputs = false;
            javaMethodEDT('dispose',frame);
            close(hidFig);
        end
        
    end

    function [] = runAlgorithmUI(frame)
        
        javaMethodEDT('setTitle',frame,'ReTrOS: Run Algorithm');
        contentPane = javaMethodEDT('getContentPane',frame);
        javaMethodEDT('removeAll',contentPane);
        
        %added 29/7/15 - DJ
        optionsPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridLayout',4,1));
        optionsPanelBorder = javaMethodEDT('createTitledBorder','javax.swing.BorderFactory','Run options:');
        javaMethodEDT('setBorder',optionsPanel,optionsPanelBorder);
        runParallelCheckBox = javaObjectEDT('javax.swing.JCheckBox','Use ''local'' parallel profile?',true);
        javaMethodEDT('setToolTipText',runParallelCheckBox,'<html>Use Parallel Computing Toolbox (PCT) if available<br>Allows multiple expression profiles to run in parallel using ''local'' profile</html>');
        runParallelCheckBoxH = handle(runParallelCheckBox,'CallbackProperties');
        set(runParallelCheckBoxH,'ActionPerformedCallback',@updateParallel);
        workersPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridLayout',1,2));
        javaMethodEDT('add',workersPanel,javaObjectEDT('javax.swing.JLabel','Number of parallel workers'));
        workersComboBox = javaObjectEDT('javax.swing.JComboBox');
        javaMethodEDT('setToolTipText',workersComboBox,'<html>Number of parallel worker threads to use<br>Default is maximum number of cores on local system</html>');
        javaMethodEDT('add',workersPanel,workersComboBox);
        v = ver('distcomp');
        if isempty(v)
            javaMethodEDT('setEnabled',runParallelCheckBox,false);
            javaMethodEDT('setSelected',runParallelCheckBox,false);
            javaMethodEDT('addItem',workersComboBox,'1');
            javaMethodEDT('setEnabled',workersComboBox,false);
        else
            c = parcluster('local');
            for aaa = 1:c.NumWorkers
                javaMethodEDT('addItem',workersComboBox,num2str(aaa));
            end
            javaMethodEDT('setSelectedIndex',workersComboBox,c.NumWorkers-1);
        end
        
        randomSeedPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridLayout',1,2));
        javaMethodEDT('add',randomSeedPanel,javaObjectEDT('javax.swing.JLabel','Random number seed'));
        randomSeedField = javaObjectEDT('javax.swing.JTextField',num2str(randomSeedDefault));
        javaMethodEDT('setToolTipText',randomSeedField,'<html>Set random number generator seed to use<br>Setting the seed allows results to be reproduced</html>');
        javaMethodEDT('add',randomSeedPanel,randomSeedField);
        
        randomSeedCheckBox = javaObjectEDT('javax.swing.JCheckBox','Generate random seed?',false);
        javaMethodEDT('setToolTipText',randomSeedCheckBox,'<html>Generate a seed for random number generator based on current time<br>If checked a random seed will be used (essential if running multiple times)</html>');
        randomSeedCheckBoxH = handle(randomSeedCheckBox,'CallbackProperties');
        set(randomSeedCheckBoxH,'ActionPerformedCallback',@updateRandomSeed);
        
        javaMethodEDT('add',optionsPanel,runParallelCheckBox);
        javaMethodEDT('add',optionsPanel,workersPanel);
        javaMethodEDT('add',optionsPanel,randomSeedPanel);
        javaMethodEDT('add',optionsPanel,randomSeedCheckBox);
        
        runPanel = javaObjectEDT('javax.swing.JPanel',javaObjectEDT('java.awt.GridLayout',3,1));
        runPanelBorder = javaMethodEDT('createTitledBorder','javax.swing.BorderFactory','Run algorithm:');
        javaMethodEDT('setBorder',runPanel,runPanelBorder);
        runSmoothCheckBox = javaObjectEDT('javax.swing.JCheckBox','Run ReTrOS-smooth?');
        javaMethodEDT('setToolTipText',runSmoothCheckBox,'Run ReTrOS-smooth using currently selected parameters');
        runSwitchCheckBox = javaObjectEDT('javax.swing.JCheckBox','Run ReTrOS-switch?');
        javaMethodEDT('setToolTipText',runSwitchCheckBox,'Run ReTrOS-switch using currently selected parameters');
        createRunScriptsCheckBox = javaObjectEDT('javax.swing.JCheckBox','Create run scripts?',false);
        javaMethodEDT('setToolTipText',createRunScriptsCheckBox,'Generate scripts to run ReTrOS with currently selected parameters');
        disp('TEMP - DISABLED SCRIPT GENERATION');
        javaMethodEDT('setEnabled',createRunScriptsCheckBox,false);
        javaMethodEDT('add',runPanel,runSmoothCheckBox);
        javaMethodEDT('add',runPanel,runSwitchCheckBox);
        javaMethodEDT('add',runPanel,createRunScriptsCheckBox);
        
        %buttons (next/back/cancel)
        buttonsPanel = javaObjectEDT('javax.swing.JPanel');
        cancelButton = javaObjectEDT('javax.swing.JButton','Close');
        runButton = javaObjectEDT('javax.swing.JButton','Run');
        % javaMethodEDT('setEnabled',nextButton,false);
        backButton = javaObjectEDT('javax.swing.JButton','Back');
        
        cancelButtonH = handle(cancelButton,'CallbackProperties');
        set(cancelButtonH,'ActionPerformedCallback',@cancelButtonEvent);
        runButtonH = handle(runButton,'CallbackProperties');
        set(runButtonH,'ActionPerformedCallback',@runButtonEvent);
        backButtonH = handle(backButton,'CallbackProperties');
        set(backButtonH,'ActionPerformedCallback',@backButtonEvent);
        
        javaMethodEDT('setLayout',buttonsPanel,javaObjectEDT('java.awt.FlowLayout'));
        javaMethodEDT('add',buttonsPanel,backButton);
        javaMethodEDT('add',buttonsPanel,cancelButton);
        javaMethodEDT('add',buttonsPanel,runButton);
        
        %add to content pane
        javaMethodEDT('add',contentPane,optionsPanel,java.awt.BorderLayout.NORTH);
        javaMethodEDT('add',contentPane,runPanel,java.awt.BorderLayout.CENTER);
        javaMethodEDT('add',contentPane,buttonsPanel,java.awt.BorderLayout.SOUTH);
        
        %display frame
        javaMethodEDT('pack',frame);
        
        function [] = runReTrOSsmooth
            
            addpath('smooth');
            
            %Parameters:
            inputFile = dataStruct.file;
            outputFile = dataStruct.outputName;
            combineReplicates = dataStruct.smoothCombineReplicates;
            nameColumnIdx = dataStruct.nameColumnSelected;
            
            params = struct;
            params.expressionType = dataStruct.dataType;
            params.bootstrapSamples = dataStruct.smoothBootstrapSamples;
            params.detrend = dataStruct.dataDetrend;
            params.kernel = dataStruct.smoothKernel;
            params.polynomialDegree = dataStruct.smoothPolyDeg;
            
            params.reporterDegRates = dataStruct.smoothDegRates;
            
            %bandwidth range
            params.bandwidthRange = dataStruct.smoothBandwidthRange;
            
            %time resolution
            params.timeResolution = dataStruct.smoothTimeResolution;
            
            useDegRateFile = dataStruct.useDegRateFile;
            degRateFile = dataStruct.degRateFile;
            
            if useDegRateFile
                disp(['Using degradation rate file: ' degRateFile]);
                d = importdata(degRateFile,'\t');
                importedDegRates = d.data;
                importedDegRateNames = d.textdata;
            else
                importedDegRates = [];
                importedDegRateNames = [];
            end
            
            %use parallel computation toolbox functions?
            useParallel = runParallelCheckBox.isSelected; % should we check the parallel pool creation at this point?
            if useParallel
                try
                    %get number of workers from user
                    numWorkers = workersComboBox.getSelectedIndex + 1;
                    pp = parpool('local',numWorkers);
                    disp('Using parpool ''local''');
                catch
                    disp('Cannot create parpool ''local''');
                    useParallel = false;
                end
            end
            
            if ~randomSeedCheckBox.isSelected
                try
                    seed = str2double(char(randomSeedField.getText));
                    
                    if isnan(seed)
                        throw(MException);
                    end
                    seed = uint32(floor(seed));
                    if seed < 0 || ~isscalar(seed) || seed>=2^32
                        throw(MException);
                    end
                catch
                    disp('Invalid seed - must be a non-negative 32-bit integer');
                    disp('Using default seed: 0');
                    seed = 0;
                end
            else
                seed = uint32(floor(sum(clock)*100));
            end
            disp(['Random number generator seed: ', num2str(seed)]);
            
            %Create output directory
            mkdir('output');
            mkdir('output',outputFile);
            mkdir(['output' filesep outputFile],'plots_smooth');
            
            %Import expression data
            contents = importdata(inputFile,'\t');
            geneNames = contents.textdata(2:end,nameColumnIdx);
            time = contents.data(1,:);
            expressionData = contents.data(2:end,:);
            tps = unique(time);
            maxReps = 0;
            for x = 1:tps
                currReps = sum(tps(x) == time);
                if currReps > maxReps
                    maxReps = currReps;
                end
            end
            
            %Main loop
            RETROS_RUNNING = true;
            if combineReplicates
                uniqueNames = uniqueRetainOrder(geneNames);
                ReTrOSsmooth_fits = cell(length(uniqueNames),1);
                plotNames = cell(length(uniqueNames),1);
                NoG = length(uniqueNames);
                
                if useParallel
                    randStreams = RandStream.create('mrg32k3a','NumStreams',NoG,'CellOutput',true,'Seed',seed);
                    
                    parfor x = 1:length(uniqueNames)
                        if RETROS_RUNNING
                            p = params;
                            gene = uniqueNames{x};
                            geneIdx = find(strcmpi(gene,geneNames));
                            data = nan(1,length(time) * length(geneIdx));
                            for y = 1:length(geneIdx)
                                data(((y-1)*length(time)+1):(y*length(time))) = expressionData(geneIdx(y),:);
                            end
                            
                            p.maximumReplicates = length(geneIdx) * maxReps;
                            
                            natDegMean = [];
                            natDegVar = [];
                            if useDegRateFile
                                degIdx = find(strcmpi(gene,importedDegRateNames(:,1)),1);
                                if ~isempty(degIdx)
                                    mRNADegMeanIdx = find(strcmpi('mRNA deg rate mean',importedDegRateNames(1,:)),1);
                                    mRNADegSDevIdx = find(strcmpi('mRNA deg rate s.dev',importedDegRateNames(1,:)),1);
                                    mRNADegMean = importedDegRates(degIdx-1,mRNADegMeanIdx-1);
                                    mRNADegSDev = importedDegRates(degIdx-1,mRNADegSDevIdx-1);
                                    if ~isempty(mRNADegMean) && ~isempty(mRNADegSDev)
                                        if ~isempty(strfind(params.expressionType,'mRNA'))
                                            p.reporterDegRates(1) = mRNADegMean;
                                            p.reporterDegRates(2) = mRNADegSDev;
                                        else
                                            natDegMean = mRNADegMean;
                                            natDegVar = mRNADegSDev * mRNADegSDev;
                                        end
                                    end
                                    
                                    if strfind(params.expressionType,'protein')
                                        proteinDegMeanIdx = find(strcmpi('protein deg rate mean',importedDegRateNames(1,:)),1);
                                        proteinDegSDevIdx = find(strcmpi('protein deg rate s.dev',importedDegRateNames(1,:)),1);
                                        proteinDegMean = importedDegRates(degIdx-1,proteinDegMeanIdx-1);
                                        proteinDegSDev = importedDegRates(degIdx-1,proteinDegSDevIdx-1);
                                        if ~isempty(proteinDegMean) && ~isempty(proteinDegSDev)
                                            p.reporterDegRates(3) = proteinDegMean;
                                            p.reporterDegRates(4) = proteinDegSDev;
                                        end
                                    end
                                end
                            end
                            
                            fprintf('%d/%d: "%s" currently calculating (%d replicates): %s\n',x,NoG,gene,p.maximumReplicates,datestr(now));
                            currFig = [];
                            try
                                warning('off','all');
                                RandStream.setGlobalStream(randStreams{x});
                                reset(randStreams{x});
                                currFig = figure('name',['ReTrOS-smooth: ' outputFile],'numbertitle','off');
                                ReTrOSsmooth_fits{x} = ReTrOSsmooth(p,data',repmat(time,1,length(geneIdx))',natDegMean,natDegVar,gene,currFig);
                                plotNames{x} = ['output' filesep outputFile filesep 'plots_smooth' filesep checkValidFileChar(gene)];
                                set(currFig,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',plotNames{x});
                                print('-dpng','-r150',plotNames{x});
                                if ~isempty(currFig)
                                    close(currFig)
                                end
                            catch exception
                                disp('Error during calculation:');
                                disp(getReport(exception));
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            end
                        end
                    end
                    delete(pp);
                else
                    %serial exectution
                    RandStream.setGlobalStream(RandStream.create('mt19937ar','Seed',seed));
                    for x = 1:length(uniqueNames)
                        if RETROS_RUNNING
                            p = params;
                            gene = uniqueNames{x};
                            geneIdx = find(strcmpi(gene,geneNames));
                            data = nan(1,length(time) * length(geneIdx));
                            for y = 1:length(geneIdx)
                                data(((y-1)*length(time)+1):(y*length(time))) = expressionData(geneIdx(y),:);
                            end
                            
                            p.maximumReplicates = length(geneIdx) * maxReps;
                            
                            natDegMean = [];
                            natDegVar = [];
                            if useDegRateFile
                                degIdx = find(strcmpi(gene,importedDegRateNames(:,1)),1);
                                if ~isempty(degIdx)
                                    mRNADegMeanIdx = find(strcmpi('mRNA deg rate mean',importedDegRateNames(1,:)),1);
                                    mRNADegSDevIdx = find(strcmpi('mRNA deg rate s.dev',importedDegRateNames(1,:)),1);
                                    mRNADegMean = importedDegRates(degIdx-1,mRNADegMeanIdx-1);
                                    mRNADegSDev = importedDegRates(degIdx-1,mRNADegSDevIdx-1);
                                    if ~isempty(mRNADegMean) && ~isempty(mRNADegSDev)
                                        if ~isempty(strfind(params.expressionType,'mRNA'))
                                            p.reporterDegRates(1) = mRNADegMean;
                                            p.reporterDegRates(2) = mRNADegSDev;
                                        else
                                            natDegMean = mRNADegMean;
                                            natDegVar = mRNADegSDev * mRNADegSDev;
                                        end
                                    end
                                    
                                    if strfind(params.expressionType,'protein')
                                        proteinDegMeanIdx = find(strcmpi('protein deg rate mean',importedDegRateNames(1,:)),1);
                                        proteinDegSDevIdx = find(strcmpi('protein deg rate s.dev',importedDegRateNames(1,:)),1);
                                        proteinDegMean = importedDegRates(degIdx-1,proteinDegMeanIdx-1);
                                        proteinDegSDev = importedDegRates(degIdx-1,proteinDegSDevIdx-1);
                                        if ~isempty(proteinDegMean) && ~isempty(proteinDegSDev)
                                            p.reporterDegRates(3) = proteinDegMean;
                                            p.reporterDegRates(4) = proteinDegSDev;
                                        end
                                    end
                                end
                            end
                            
                            fprintf('%d/%d: "%s" currently calculating (%d replicates): %s\n',x,NoG,gene,p.maximumReplicates,datestr(now));
                            try
                                rng(seed);
                                currFig = figure('name',['ReTrOS-smooth: ' outputFile],'numbertitle','off');
                                ReTrOSsmooth_fits{x} = ReTrOSsmooth(p,data',repmat(time,1,length(geneIdx))',natDegMean,natDegVar,gene,currFig);
                                plotNames{x} = ['output' filesep outputFile filesep 'plots_smooth' filesep checkValidFileChar(gene)];
                                set(currFig,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',plotNames{x});
                                print('-dpng','-r150',plotNames{x});
                                if ~isempty(currFig)
                                    close(currFig)
                                end
                            catch exception
                                disp('Error during calculation:');
                                disp(getReport(exception));
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            end
                        end
                    end
                end
                geneNames = uniqueNames;
            else
                ReTrOSsmooth_fits = cell(length(geneNames),1);
                plotNames = cell(length(geneNames),1);
                geneNamesOrig = geneNames;
                NoG = length(geneNames);
                if useParallel
                    randStreams = RandStream.create('mrg32k3a','NumStreams',NoG,'CellOutput',true,'Seed',seed);
                    
                    parfor x = 1:NoG
                        if RETROS_RUNNING %doesnt terminate execution in parallel
                            p = params;
                            p.maximumReplicates = maxReps;
                            
                            gene = geneNames{x};
                            
                            natDegMean = [];
                            natDegVar = [];
                            if useDegRateFile
                                degIdx = find(strcmpi(gene,importedDegRateNames(:,1)),1);
                                if ~isempty(degIdx)
                                    mRNADegMeanIdx = find(strcmpi('mRNA deg rate mean',importedDegRateNames(1,:)),1);
                                    mRNADegSDevIdx = find(strcmpi('mRNA deg rate s.dev',importedDegRateNames(1,:)),1);
                                    mRNADegMean = importedDegRates(degIdx-1,mRNADegMeanIdx-1);
                                    mRNADegSDev = importedDegRates(degIdx-1,mRNADegSDevIdx-1);
                                    if ~isempty(mRNADegMean) && ~isempty(mRNADegSDev)
                                        if ~isempty(strfind(params.expressionType,'mRNA'))
                                            p.reporterDegRates(1) = mRNADegMean;
                                            p.reporterDegRates(2) = mRNADegSDev;
                                        else
                                            natDegMean = mRNADegMean;
                                            natDegVar = mRNADegSDev * mRNADegSDev;
                                        end
                                    end
                                    
                                    if strfind(params.expressionType,'protein')
                                        proteinDegMeanIdx = find(strcmpi('protein deg rate mean',importedDegRateNames(1,:)),1);
                                        proteinDegSDevIdx = find(strcmpi('protein deg rate s.dev',importedDegRateNames(1,:)),1);
                                        proteinDegMean = importedDegRates(degIdx-1,proteinDegMeanIdx-1);
                                        proteinDegSDev = importedDegRates(degIdx-1,proteinDegSDevIdx-1);
                                        if ~isempty(proteinDegMean) && ~isempty(proteinDegSDev)
                                            p.reporterDegRates(3) = proteinDegMean;
                                            p.reporterDegRates(4) = proteinDegSDev;
                                        end
                                    end
                                end
                            end
                            
                            if sum(strcmpi(gene,geneNamesOrig)) > 1
                                geneNum = sum(strcmpi(gene,geneNamesOrig(1:x)));
                                gene = [gene '_' num2str(geneNum)];
                            end
                            geneNames{x} = gene;
                            
                            fprintf('%d/%d: "%s" currently calculating (%d replicates): %s\n',x,NoG,gene,p.maximumReplicates,datestr(now));
                            currFig = [];
                            try
                                warning('off','all');
                                RandStream.setGlobalStream(randStreams{x});
                                reset(randStreams{x});
                                currFig = figure('name',['ReTrOS-smooth: ' outputFile],'numbertitle','off');
                                ReTrOSsmooth_fits{x} = ReTrOSsmooth(p,expressionData(x,:)',time',natDegMean,natDegVar,gene,currFig);
                                plotNames{x} = ['output' filesep outputFile filesep 'plots_smooth' filesep checkValidFileChar(gene)];
                                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',plotNames{x});
                                print('-dpng','-r150',plotNames{x});
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            catch exception
                                disp('Error during calculation:');
                                disp(getReport(exception));
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            end
                        end
                    end
                    delete(pp);
                else
                    %serial exectution
                    RandStream.setGlobalStream(RandStream.create('mt19937ar','Seed',seed));
                    for x = 1:NoG
                        if RETROS_RUNNING %doesnt terminate execution in parallel
                            p = params;
                            p.maximumReplicates = maxReps;
                            
                            gene = geneNames{x};
                            
                            natDegMean = [];
                            natDegVar = [];
                            if useDegRateFile
                                degIdx = find(strcmpi(gene,importedDegRateNames(:,1)),1);
                                if ~isempty(degIdx)
                                    mRNADegMeanIdx = find(strcmpi('mRNA deg rate mean',importedDegRateNames(1,:)),1);
                                    mRNADegSDevIdx = find(strcmpi('mRNA deg rate s.dev',importedDegRateNames(1,:)),1);
                                    mRNADegMean = importedDegRates(degIdx-1,mRNADegMeanIdx-1);
                                    mRNADegSDev = importedDegRates(degIdx-1,mRNADegSDevIdx-1);
                                    if ~isempty(mRNADegMean) && ~isempty(mRNADegSDev)
                                        if ~isempty(strfind(params.expressionType,'mRNA'))
                                            p.reporterDegRates(1) = mRNADegMean;
                                            p.reporterDegRates(2) = mRNADegSDev;
                                        else
                                            natDegMean = mRNADegMean;
                                            natDegVar = mRNADegSDev * mRNADegSDev;
                                        end
                                    end
                                    
                                    if strfind(params.expressionType,'protein')
                                        proteinDegMeanIdx = find(strcmpi('protein deg rate mean',importedDegRateNames(1,:)),1);
                                        proteinDegSDevIdx = find(strcmpi('protein deg rate s.dev',importedDegRateNames(1,:)),1);
                                        proteinDegMean = importedDegRates(degIdx-1,proteinDegMeanIdx-1);
                                        proteinDegSDev = importedDegRates(degIdx-1,proteinDegSDevIdx-1);
                                        if ~isempty(proteinDegMean) && ~isempty(proteinDegSDev)
                                            p.reporterDegRates(3) = proteinDegMean;
                                            p.reporterDegRates(4) = proteinDegSDev;
                                        end
                                    end
                                end
                            end
                            
                            if sum(strcmpi(gene,geneNamesOrig)) > 1
                                geneNum = sum(strcmpi(gene,geneNamesOrig(1:x)));
                                gene = [gene '_' num2str(geneNum)];
                            end
                            geneNames{x} = gene;
                            
                            fprintf('%d/%d: "%s" currently calculating (%d replicates): %s\n',x,NoG,gene,p.maximumReplicates,datestr(now));
                            try
                                rng(seed);
                                currFig = figure('name',['ReTrOS-smooth: ' outputFile],'numbertitle','off');
                                ReTrOSsmooth_fits{x} = ReTrOSsmooth(p,expressionData(x,:)',time',natDegMean,natDegVar,gene,currFig);
                                plotNames{x} = ['output' filesep outputFile filesep 'plots_smooth' filesep checkValidFileChar(gene)];
                                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',plotNames{x});
                                print('-dpng','-r150',plotNames{x});
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            catch exception
                                disp('Error during calculation:');
                                disp(getReport(exception));
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            end
                        end
                    end
                end
            end
            
            disp('Creating output');
            save(['output' filesep outputFile filesep outputFile '_smooth.mat'],'ReTrOSsmooth_fits','geneNames','inputFile');
            
            if ~ispc
                pdfFile = ['output' filesep outputFile filesep outputFile '_smooth.pdf'];
                paramsFile = ['output' filesep outputFile filesep 'params.txt'];
                paramsFid = fopen(paramsFile, 'w');
                fprintf(paramsFid,'-dNOPAUSE\n');
                fprintf(paramsFid,'-dBATCH\n');
                fprintf(paramsFid,'-sOutputFile="%s"\n', pdfFile);
                fprintf(paramsFid,'-sDEVICE=%s\n', 'pdfwrite');
                for m = 1:length(plotNames)
                    if ~isempty(plotNames{m})
                        fprintf(paramsFid,['"' plotNames{m} '.pdf"\n']);
                    end
                end
                fclose(paramsFid);
                [status, result] = system(['/usr/local/bin/gs @"' paramsFile '"']);
                delete(paramsFile);
            end
            
            disp(['Finished (' datestr(clock) ')']);
            
            RETROS_RUNNING = false;
            rmpath('smooth');
        end
        
        function [] = runReTrOSswitch
            addpath('switch');
            %parameters:
            %input file
            inputFile = dataStruct.file ;
            %output name
            outputFile = dataStruct.outputName;
            combineReplicates = dataStruct.switchCombineReplicates;
            nameColumnIdx = dataStruct.nameColumnSelected;
            
            %expression type
            expressionType = dataStruct.dataType;
            dataDetrend = dataStruct.dataDetrend;
            regressionMethod = dataStruct.switchRegressionMethod;
            mRNADegRates = dataStruct.switchmRNADegRates;
            proteinDegRates = dataStruct.switchProteinDegRates;
            MCMCiterations = dataStruct.switchIterations;
            initialSwitchGuess = dataStruct.switchIntialSwitchGuess;
            maximumSwitchNum = dataStruct.switchMaxSwitchNum;
            expectedSwitchNum = dataStruct.switchExpSwitchNum;
            switchDelayTime = dataStruct.switchDelayTime;
            burnIn = dataStruct.switchBurnIn;
            numSamplesEstPost = dataStruct.switchNumSamplesEstPost;
            baselineFactor = dataStruct.switchBaselineFactor;
            maxSwDevFactor = dataStruct.switchMaxSwDev;
            minSwStr = dataStruct.switchMinSwStr;
            numSamplePlot = dataStruct.switchNumSamplePlot;
            
            %use parallel computation toolbox functions?
            useParallel =  runParallelCheckBox.isSelected;
            if useParallel
                try
                    %get number of workers from user
                    numWorkers = workersComboBox.getSelectedIndex + 1;
                    pp = parpool('local',numWorkers);
                    disp('Using parpool ''local''');
                catch
                    disp('Cannot create parpool ''local''');
                    useParallel = false;
                end
            end
            
            if ~randomSeedCheckBox.isSelected
                try
                    seed = str2double(char(randomSeedField.getText));
                    
                    if isnan(seed)
                        throw(MException);
                    end
                    seed = uint32(floor(seed));
                    if seed < 0 || ~isscalar(seed) || seed>=2^32
                        throw(MException);
                    end
                catch
                    disp('Invalid seed - must be a non-negative 32-bit integer');
                    disp('Using default seed: 0');
                    seed = 0;
                end
            else
                seed = uint32(floor(sum(clock)*100));
            end
            disp(['Random number generator seed: ', num2str(seed)]);
            
            mkdir('output');
            mkdir('output',outputFile);
            mkdir(['output' filesep outputFile],'plots_switch');
            
            contents = importdata(inputFile,'\t');
            geneNames = contents.textdata(2:end,nameColumnIdx);
            time = contents.data(1,:);
            expressionData = contents.data(2:end,:);
            timescale = unique(time);
            
            maxReps = 0;
            for x = 1:length(timescale)
                if sum(timescale(x) == time) > maxReps
                    maxReps = sum(timescale(x) == time);
                end
            end
            %Main loop
            RETROS_RUNNING = true;
            if combineReplicates
                uniqueNames = uniqueRetainOrder(geneNames);
                NoG = length(uniqueNames);
                ReTrOSswitch_fits = cell(NoG,5);
                plotNames = cell(length(uniqueNames),2);
                
                if useParallel
                    randStreams = RandStream.create('mrg32k3a','NumStreams',NoG,'CellOutput',true,'Seed',seed);
                    
                    parfor x = 1:NoG
                        if RETROS_RUNNING  %doesnt terminate execution in parallel computation
                            gene = uniqueNames{x};
                            geneIdx = find(strcmpi(gene,geneNames));
                            data = nan(1,length(timescale) * maxReps * length(geneIdx));
                            for y = 1:length(geneIdx)
                                offset = (y-1) * maxReps * length(timescale);
                                for z = 1:length(timescale)
                                    d = expressionData(geneIdx(y),timescale(z)==time);
                                    for w = 1:length(d)
                                        data(offset + (w-1) * length(timescale) + z) = d(w);
                                    end
                                end
                            end
                            fprintf('%d/%d: "%s" currently calculating (%d replicates): %s\n',x,NoG,gene,maxReps*length(geneIdx),datestr(now));
                            currFig = [];
                            %run switch model
                            %mRNA
                            try
                                warning('off','all');
                                RandStream.setGlobalStream(randStreams{x});
                                reset(randStreams{x});
                                currFig = figure('name',['ReTrOS-switch: ',outputFile,' - ',gene],'numbertitle','off');
                                if strfind(expressionType,'mRNA')
                                    [fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_mRNA(gene,data,timescale,maxReps*length(geneIdx),MCMCiterations,'detrendData',dataDetrend,'regressionMethod',regressionMethod,'mRNADegradationRatePrior',mRNADegRates,'useInitialSwitchSuggest',initialSwitchGuess,'maxSwitches',maximumSwitchNum,'switchNumberPrior',expectedSwitchNum,'switchDelay',switchDelayTime,'burnIn',burnIn,'maxSamplesToUse',numSamplesEstPost,'switchBaselineFactor',baselineFactor,'minSwitchStrength',minSwStr,'maxSwitchDeviationFactor',maxSwDevFactor,'maxSamplesToPlot',numSamplePlot,'figHandle',currFig);
                                    %protein
                                else
                                    [fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_protein(gene,data,timescale,maxReps*length(geneIdx),MCMCiterations,'detrendData',dataDetrend,'regressionMethod',regressionMethod,'mRNADegradationRatePrior',mRNADegRates,'proteinDegradationRatePrior',proteinDegRates,'useInitialSwitchSuggest',initialSwitchGuess,'maxSwitches',maximumSwitchNum,'switchNumberPrior',expectedSwitchNum,'switchDelay',switchDelayTime,'burnIn',burnIn,'maxSamplesToUse',numSamplesEstPost,'switchBaselineFactor',baselineFactor,'minSwitchStrength',minSwStr,'maxSwitchDeviationFactor',maxSwDevFactor,'maxSamplesToPlot',numSamplePlot,'figHandle',currFig);
                                end
                                ReTrOSswitch_fits(x,:) = {gene, fitData, fitSamples, fitPosteriors, fitParams};
                                
                                summaryPlotName = ['output' filesep outputFile filesep 'plots_switch' filesep checkValidFileChar(gene)];
                                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',summaryPlotName);
                                print('-dpng','-r150',summaryplotName);
                                
                                plotSwitchModelsBasic(fitPosteriors,fitSamples,fitParams,fitData,currFig);
                                modelPlotName = ['output' filesep outputFile filesep 'plots_switch' filesep checkValidFileChar(gene) '_models'];
                                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',modelPlotName);
                                print('-dpng','-r150',modelPlotName);
                                plotNames(x,:) = {summaryPlotName, modelPlotName};
                                
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            catch exception
                                disp('Error during calculation:');
                                disp(getReport(exception));
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            end
                        end
                    end
                    delete(pp);
                else
                    %serial execution of models
                    RandStream.setGlobalStream(RandStream.create('mt19937ar','Seed',seed));
                    for x = 1:NoG
                        if RETROS_RUNNING
                            gene = uniqueNames{x};
                            geneIdx = find(strcmpi(gene,geneNames));
                            data = nan(1,length(timescale) * maxReps * length(geneIdx));
                            for y = 1:length(geneIdx)
                                offset = (y-1) * maxReps * length(timescale);
                                for z = 1:length(timescale)
                                    d = expressionData(geneIdx(y),timescale(z)==time);
                                    for w = 1:length(d)
                                        data(offset + (w-1) * length(timescale) + z) = d(w);
                                    end
                                end
                            end
                            fprintf('%d/%d: "%s" currently calculating (%d replicates): %s\n',x,NoG,gene,maxReps*length(geneIdx),datestr(now));
                            %run switch model
                            %mRNA
                            try
                                rng(seed);
                                currFig = figure('name',['ReTrOS-switch: ',outputFile,' - ',gene],'numbertitle','off');
                                if strfind(expressionType,'mRNA')
                                    [fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_mRNA(gene,data,timescale,maxReps*length(geneIdx),MCMCiterations,'detrendData',dataDetrend,'regressionMethod',regressionMethod,'mRNADegradationRatePrior',mRNADegRates,'useInitialSwitchSuggest',initialSwitchGuess,'maxSwitches',maximumSwitchNum,'switchNumberPrior',expectedSwitchNum,'switchDelay',switchDelayTime,'burnIn',burnIn,'maxSamplesToUse',numSamplesEstPost,'switchBaselineFactor',baselineFactor,'minSwitchStrength',minSwStr,'maxSwitchDeviationFactor',maxSwDevFactor,'maxSamplesToPlot',numSamplePlot,'figHandle',currFig);
                                    %protein
                                else
                                    [fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_protein(gene,data,timescale,maxReps*length(geneIdx),MCMCiterations,'detrendData',dataDetrend,'regressionMethod',regressionMethod,'mRNADegradationRatePrior',mRNADegRates,'proteinDegradationRatePrior',proteinDegRates,'useInitialSwitchSuggest',initialSwitchGuess,'maxSwitches',maximumSwitchNum,'switchNumberPrior',expectedSwitchNum,'switchDelay',switchDelayTime,'burnIn',burnIn,'maxSamplesToUse',numSamplesEstPost,'switchBaselineFactor',baselineFactor,'minSwitchStrength',minSwStr,'maxSwitchDeviationFactor',maxSwDevFactor,'maxSamplesToPlot',numSamplePlot,'figHandle',currFig);
                                end
                                ReTrOSswitch_fits(x,:) = {gene, fitData, fitSamples, fitPosteriors, fitParams};
                                
                                summaryPlotName = ['output' filesep outputFile filesep 'plots_switch' filesep checkValidFileChar(gene)];
                                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',summaryPlotName);
                                print('-dpng','-r150',summaryPlotName);
                                
                                plotSwitchModelsBasic(fitPosteriors,fitSamples,fitParams,fitData,currFig);
                                modelPlotName = ['output' filesep outputFile filesep 'plots_switch' filesep checkValidFileChar(gene) '_models'];
                                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',modelPlotName);
                                print('-dpng','-r150',modelPlotName);
                                plotNames(x,:) = {summaryPlotName, modelPlotName};
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            catch exception
                                disp('Error during calculation:');
                                disp(getReport(exception));
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            end
                        end
                    end
                end
                geneNames = uniqueNames;
            else
                NoG = size(expressionData,1);
                ReTrOSswitch_fits = cell(NoG,5);
                plotNames = cell(NoG,2);
                data = nan(NoG,length(timescale * maxReps));
                for x = 1:length(timescale)
                    t = find(time == timescale(x));
                    for y = 1:length(t)
                        data(:,(y-1) * length(timescale) + x) = expressionData(:,t(y));
                    end
                end
                geneNamesOrig = geneNames;
                
                if useParallel
                    randStreams = RandStream.create('mrg32k3a','NumStreams',NoG,'CellOutput',true,'Seed',seed);
                    
                    parfor x = 1:NoG
                        if RETROS_RUNNING % doesnt work with parallel for early termination
                            gene = geneNames{x};
                            if sum(strcmpi(geneNamesOrig,gene)) > 1
                                geneNum = sum(strcmpi(gene,geneNamesOrig(1:x)));
                                gene = [gene '_' num2str(geneNum)];
                            end
                            geneNames{x} = gene;
                            fprintf('%d/%d: "%s" currently calculating (%d replicates): %s\n',x,NoG,gene,maxReps,datestr(now));
                            currFig = [];
                            %run switch model
                            %mRNA
                            try
                                warning('off','all');
                                RandStream.setGlobalStream(randStreams{x});
                                reset(randStreams{x});
                                currFig = figure('name',['ReTrOS-switch: ',outputFile,' - ',gene],'numbertitle','off');
                                if strfind(expressionType,'mRNA')
                                    [fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_mRNA(gene,data(x,:),timescale,maxReps,MCMCiterations,'detrendData',dataDetrend,'regressionMethod',regressionMethod,'mRNADegradationRatePrior',mRNADegRates,'useInitialSwitchSuggest',initialSwitchGuess,'maxSwitches',maximumSwitchNum,'switchNumberPrior',expectedSwitchNum,'switchDelay',switchDelayTime,'burnIn',burnIn,'maxSamplesToUse',numSamplesEstPost,'switchBaselineFactor',baselineFactor,'minSwitchStrength',minSwStr,'maxSwitchDeviationFactor',maxSwDevFactor,'maxSamplesToPlot',numSamplePlot,'figHandle',currFig);
                                    %protein
                                else
                                    [fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_protein(gene,data(x,:),timescale,maxReps,MCMCiterations,'detrendData',dataDetrend,'regressionMethod',regressionMethod,'mRNADegradationRatePrior',mRNADegRates,'proteinDegradationRatePrior',proteinDegRates,'useInitialSwitchSuggest',initialSwitchGuess,'maxSwitches',maximumSwitchNum,'switchNumberPrior',expectedSwitchNum,'switchDelay',switchDelayTime,'burnIn',burnIn,'maxSamplesToUse',numSamplesEstPost,'switchBaselineFactor',baselineFactor,'minSwitchStrength',minSwStr,'maxSwitchDeviationFactor',maxSwDevFactor,'maxSamplesToPlot',numSamplePlot,'figHandle',currFig);
                                end
                                
                                ReTrOSswitch_fits(x,:) = {gene, fitData, fitSamples, fitPosteriors, fitParams};
                                
                                summaryPlotName = ['output' filesep outputFile filesep 'plots_switch' filesep checkValidFileChar(gene)];
                                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',summaryPlotName);
                                print('-dpng','-r150',summaryPlotName);
                                
                                plotSwitchModelsBasic(fitPosteriors,fitSamples,fitParams,fitData,currFig);
                                modelPlotName = ['output' filesep outputFile filesep 'plots_switch' filesep checkValidFileChar(gene) '_models'];
                                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',modelPlotName);
                                print('-dpng','-r150',modelPlotName);
                                
                                plotNames(x,:) = {summaryPlotName, modelPlotName};
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            catch exception
                                disp('Error during calculation:');
                                disp(getReport(exception));
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            end
                        end
                    end
                    delete(pp);
                else
                    %serial execution of loop
                    RandStream.setGlobalStream(RandStream.create('mt19937ar','Seed',seed));
                    for x = 1:NoG
                        if RETROS_RUNNING
                            gene = geneNames{x};
                            if sum(strcmpi(geneNamesOrig,gene)) > 1
                                geneNum = sum(strcmpi(gene,geneNamesOrig(1:x)));
                                gene = [gene '_' num2str(geneNum)];
                            end
                            geneNames{x} = gene;
                            fprintf('%d/%d: "%s" currently calculating (%d replicates): %s\n',x,NoG,gene,maxReps,datestr(now));
                            %run switch model
                            %mRNA
                            try
                                rng(seed);
                                currFig = figure('name',['ReTrOS-switch: ',outputFile,' - ',gene],'numbertitle','off');
                                if strfind(expressionType,'mRNA')
                                    [fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_mRNA(gene,data(x,:),timescale,maxReps,MCMCiterations,'detrendData',dataDetrend,'regressionMethod',regressionMethod,'mRNADegradationRatePrior',mRNADegRates,'useInitialSwitchSuggest',initialSwitchGuess,'maxSwitches',maximumSwitchNum,'switchNumberPrior',expectedSwitchNum,'switchDelay',switchDelayTime,'burnIn',burnIn,'maxSamplesToUse',numSamplesEstPost,'switchBaselineFactor',baselineFactor,'minSwitchStrength',minSwStr,'maxSwitchDeviationFactor',maxSwDevFactor,'maxSamplesToPlot',numSamplePlot,'figHandle',currFig);
                                    %protein
                                else
                                    [fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_protein(gene,data(x,:),timescale,maxReps,MCMCiterations,'detrendData',dataDetrend,'regressionMethod',regressionMethod,'mRNADegradationRatePrior',mRNADegRates,'proteinDegradationRatePrior',proteinDegRates,'useInitialSwitchSuggest',initialSwitchGuess,'maxSwitches',maximumSwitchNum,'switchNumberPrior',expectedSwitchNum,'switchDelay',switchDelayTime,'burnIn',burnIn,'maxSamplesToUse',numSamplesEstPost,'switchBaselineFactor',baselineFactor,'minSwitchStrength',minSwStr,'maxSwitchDeviationFactor',maxSwDevFactor,'maxSamplesToPlot',numSamplePlot,'figHandle',currFig);
                                end
                                
                                ReTrOSswitch_fits(x,:) = {gene, fitData, fitSamples, fitPosteriors, fitParams};
                                
                                summaryPlotName = ['output' filesep outputFile filesep 'plots_switch' filesep checkValidFileChar(gene)];
                                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',summaryPlotName);
                                print('-dpng','-r150',summaryPlotName);
                                
                                plotSwitchModelsBasic(fitPosteriors,fitSamples,fitParams,fitData,currFig);
                                modelPlotName = ['output' filesep outputFile filesep 'plots_switch' filesep checkValidFileChar(gene) '_models'];
                                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                                print('-dpdf','-r150',modelPlotName);
                                print('-dpng','-r150',modelPlotName);
                                
                                plotNames(x,:) = {summaryPlotName, modelPlotName};
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            catch exception
                                disp('Error during calculation:');
                                disp(getReport(exception));
                                if ~isempty(currFig)
                                    close(currFig);
                                end
                            end
                        end
                    end
                end
            end
            disp('Creating output');
            %output switch times and switch times plot
            swTimesFile = ['output' filesep outputFile filesep outputFile '_switchTimes.txt'];
            swTimesFID = fopen(swTimesFile, 'w');
            fprintf(swTimesFID,'gene\tswitch type\tswitch mean time\tswitch standard deviation\tswitch strength\n');
            for x = 1:size(ReTrOSswitch_fits,1)
                posteriors = ReTrOSswitch_fits{x,4};
                gene = ReTrOSswitch_fits{x,1};
                if ~isempty(posteriors)
                    swTimes = posteriors.switchDistribution.means;
                    swSigmas = posteriors.switchDistribution.standardDeviations;
                    swStrengths = posteriors.switchDistribution.strengths;
                    rates = posteriors.modelFit.taus_medianModel;
                    for z = 1:length(swTimes)
                        if (rates(z+1)-rates(z)) > 0
                            type = 'increase';
                        else
                            type = 'decrease';
                        end
                        fprintf(swTimesFID,[gene '\t' type '\t' num2str(swTimes(z)) '\t' num2str(swSigmas(z)) '\t' num2str(swStrengths(z)) '\n']);
                    end
                end
            end
            fclose(swTimesFID);
            heatmapH = figure;
            outputSwPlot = plotSwitchHeatmap(swTimesFile,timescale,heatmapH);
            outputSwPlotName = [];
            if ~isempty(outputSwPlot)
                outputSwPlotName = ['output' filesep outputFile filesep 'plots_switch' filesep 'switchTimes'];
                set(gcf,'PaperUnits','centimeters','PaperSize',[29.7,21],'PaperPosition',[0,0,29.7,21]);
                print('-dpdf','-r150',outputSwPlotName);
                print('-dpng','-r150',outputSwPlotName);
            end
            if ~isempty(heatmapH)
                close(heatmapH);
            end
            %output degradation rates
            degRatesFile = ['output' filesep outputFile filesep outputFile '_degRates.txt'];
            degRatesFID = fopen(degRatesFile,'w');
            if strfind(expressionType,'mRNA')
                fprintf(degRatesFID,'gene\tmRNA deg rate mean\tmRNA deg rate s.dev\n');
                for x = 1:size(ReTrOSswitch_fits,1)
                    posteriors = ReTrOSswitch_fits{x,4};
                    gene = ReTrOSswitch_fits{x,1};
                    if ~isempty(posteriors)
                        fprintf(degRatesFID,[gene '\t' num2str(posteriors.degradationRate_mRNA.mean) '\t' num2str(posteriors.degradationRate_mRNA.standardDeviation) '\n']);
                    end
                end
            else
                fprintf(degRatesFID,'gene\tmRNA deg rate mean\tmRNA deg rate s.dev\tprotein deg rate mean\tprotein deg rate s.dev\n');
                for x = 1:size(ReTrOSswitch_fits,1)
                    posteriors = ReTrOSswitch_fits{x,4};
                    gene = ReTrOSswitch_fits{x,1};
                    if ~isempty(posteriors)
                        fprintf(degRatesFID,[gene '\t' num2str(posteriors.degradationRate_mRNA.mean) '\t' num2str(posteriors.degradationRate_mRNA.standardDeviation) '\t' num2str(posteriors.degradationRate_protein.mean) '\t' num2str(posteriors.degradationRate_protein.standardDeviation) '\n']);
                    end
                end
            end
            fclose(degRatesFID);
            %save outputs as MAT file
            save(['output' filesep outputFile filesep outputFile '_switch.mat'],'inputFile','geneNames','ReTrOSswitch_fits');
            %if on a unix machine, try to combine plots using GhostScript
            if ~ispc
                pdfFile = ['output' filesep outputFile filesep outputFile '_switch.pdf'];
                paramsFile = ['output' filesep outputFile filesep outputFile '_params.txt'];
                paramsFid = fopen(paramsFile, 'w');
                fprintf(paramsFid,'-dNOPAUSE\n');
                fprintf(paramsFid,'-dBATCH\n');
                fprintf(paramsFid,'-sOutputFile="%s"\n', pdfFile);
                fprintf(paramsFid,'-sDEVICE=%s\n', 'pdfwrite');
                for m = 1:size(ReTrOSswitch_fits,1)
                    if ~isempty(plotNames{m,1})
                        fprintf(paramsFid,['"' plotNames{m,1} '.pdf"\n']);
                    end
                    if ~isempty(plotNames{m,2})
                        fprintf(paramsFid,['"' plotNames{m,2} '.pdf"\n']);
                    end
                end
                if ~isempty(outputSwPlotName)
                    fprintf(paramsFid,['"' outputSwPlotName '.pdf"\n']);
                end
                fclose(paramsFid);
                [status, result] = system(['/usr/local/bin/gs @"' paramsFile '"']);
                delete(paramsFile);
            end
            disp(['Finished (' datestr(clock) ')']);
            RETROS_RUNNING = false;
            rmpath('switch');
        end
        
        %Create scripts functions
        function [] = createSmoothScript
            
            if exist(dataStruct.smoothScriptFile,'file')
                disp(['ReTrOS-smooth: file "' dataStruct.smoothScriptFile '" already exists - overwriting']);
            end
            scriptFID = fopen(dataStruct.smoothScriptFile,'w');
            fprintf(scriptFID,['function [ReTrOSsmooth_fits] = run_' dataStruct.outputName '_smooth\n']);
            fprintf(scriptFID,['\n%%Generated ReTrOS-smooth script (' datestr(clock) ')\n']);
            fprintf(scriptFID,'\ncurrentWarnings = warning(''off'',''all'');\n');
            fprintf(scriptFID,'addpath([''..'' filesep ''smooth'']);\n');
            fprintf(scriptFID,'\n%%Parameters:\n');
            fprintf(scriptFID,['inputFile = ''' dataStruct.file ''';\n']);
            fprintf(scriptFID,['outputFile = ''' dataStruct.outputName ''';\n']);
            fprintf(scriptFID,['combineReplicates = ' logical2str(dataStruct.smoothCombineReplicates) ';\n']);
            fprintf(scriptFID,['nameColumnIdx = ' num2str(dataStruct.nameColumnSelected) ';\n']);
            
            fprintf(scriptFID,'\nparams = struct;\n');
            fprintf(scriptFID,['params.expressionType = ''' dataStruct.dataType ''';\n']);
            fprintf(scriptFID,['params.bootstrapSamples = ' num2str(dataStruct.smoothBootstrapSamples) ';\n']);
            fprintf(scriptFID,['params.detrend = ''' dataStruct.dataDetrend ''';\n']);
            fprintf(scriptFID,['params.kernel = ''' dataStruct.smoothKernel ''';\n']);
            fprintf(scriptFID,['params.polynomialDegree = ' num2str(dataStruct.smoothPolyDeg) ';\n']);
            fprintf(scriptFID,['params.reporterDegRates = [' num2str(dataStruct.smoothDegRates(1)) ' ' num2str(dataStruct.smoothDegRates(2)) ' ' num2str(dataStruct.smoothDegRates(3)) ' ' num2str(dataStruct.smoothDegRates(4)) '];\n']);
            fprintf(scriptFID,['params.bandwidthRange = [' num2str(dataStruct.smoothBandwidthRange(1)) ' ' num2str(dataStruct.smoothBandwidthRange(2)) ' ' num2str(dataStruct.smoothBandwidthRange(3)) '];\n']);
            fprintf(scriptFID,['params.timeResolution = ' num2str(dataStruct.smoothTimeResolution) ';\n']);
            fprintf(scriptFID,['params.useDegRateFile = ' logical2str(dataStruct.useDegRateFile) ';\n']);
            fprintf(scriptFID,['params.degRateFile = ''' dataStruct.degRateFile ''';\n']);
            
            fprintf(scriptFID,'\nif params.useDegRateFile\n');
            fprintf(scriptFID,'\tdisp([''Using degradation rate: '' params.degRateFile]);\n');
            fprintf(scriptFID,'\td = importdata(params.degRateFile,''\\t'');\n');
            fprintf(scriptFID,'\timportedDegRates = d.data;\n');
            fprintf(scriptFID,'\timportedDegRateNames = d.textdata;\n');
            fprintf(scriptFID,'end\n');
            
            fprintf(scriptFID,'\n%%Create output directory\n');
            fprintf(scriptFID,'mkdir([''..'' filesep ''output'']);\n');
            fprintf(scriptFID,'mkdir([''..'' filesep ''output''],outputFile);\n');
            fprintf(scriptFID,'mkdir([''..'' filesep ''output'' filesep outputFile],''plots_smooth'');\n');
            
            fprintf(scriptFID,'\n%%Import expression data\n');
            fprintf(scriptFID,'contents = importdata(inputFile,''\\t'');\n');
            fprintf(scriptFID,'geneNames = contents.textdata(2:end,nameColumnIdx);\n');
            fprintf(scriptFID,'time = contents.data(1,:);\n');
            fprintf(scriptFID,'expressionData = contents.data(2:end,:);\n');
            fprintf(scriptFID,'tps = unique(time);\n');
            fprintf(scriptFID,'maxReps = 0;\n');
            fprintf(scriptFID,'for x = 1:tps\n');
            fprintf(scriptFID,'\tcurrReps = sum(tps(x) == time);\n');
            fprintf(scriptFID,'\tif currReps > maxReps\n');
            fprintf(scriptFID,'\tmaxReps = currReps;\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'end\n');
            
            fprintf(scriptFID,'\nfigH = figure(''name'',[''ReTrOS-smooth: '' outputFile],''numbertitle'',''off'');\n');
            fprintf(scriptFID,'%%Main loop\n');
            fprintf(scriptFID,'if combineReplicates\n');
            fprintf(scriptFID,'\tuniqueNames = uniqueRetainOrder(geneNames);\n');
            fprintf(scriptFID,'\tReTrOSsmooth_fits = cell(length(uniqueNames),1);\n');
            fprintf(scriptFID,'\tplotNames = cell(length(uniqueNames),1);\n');
            fprintf(scriptFID,'\tfor x = 1:length(uniqueNames)\n');
            fprintf(scriptFID,'\t\tp = params;\n');
            fprintf(scriptFID,'\t\tgene = uniqueNames{x};\n');
            fprintf(scriptFID,'\t\tgeneIdx = find(strcmpi(gene,geneNames));\n');
            fprintf(scriptFID,'\t\tdata = nan(1,length(time) * length(geneIdx));\n');
            fprintf(scriptFID,'\t\tfor y = 1:length(geneIdx)\n');
            fprintf(scriptFID,'\t\t\tdata(((y-1)*length(time)+1):(y*length(time))) = expressionData(geneIdx(y),:);\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\n\t\tp.maximumReplicates = length(geneIdx) * maxReps;\n');
            fprintf(scriptFID,'\n\t\tnatDegMean = [];\n');
            fprintf(scriptFID,'\n\t\tnatDegVar = [];\n');
            fprintf(scriptFID,'\t\tif params.useDegRateFile\n');
            fprintf(scriptFID,'\t\t\tdegIdx = find(strcmpi(gene,importedDegRateNames(:,1)),1);\n');
            fprintf(scriptFID,'\t\t\tif ~isempty(degIdx)\n');
            fprintf(scriptFID,'\t\t\t\tmRNADegMeanIdx = find(strcmpi(''mRNA deg rate mean'',importedDegRateNames(1,:)),1);\n');
            fprintf(scriptFID,'\t\t\t\tmRNADegSDevIdx = find(strcmpi(''mRNA deg rate s.dev'',importedDegRateNames(1,:)),1);\n');
            fprintf(scriptFID,'\t\t\t\tmRNADegMean = importedDegRates(degIdx-1,mRNADegMeanIdx-1);\n');
            fprintf(scriptFID,'\t\t\t\tmRNADegSDev = importedDegRates(degIdx-1,mRNADegSDevIdx-1);\n');
            fprintf(scriptFID,'\t\t\t\tif ~isempty(mRNADegMean) && ~isempty(mRNADegSDev)\n');
            fprintf(scriptFID,'\t\t\t\t\tif ~isempty(strfind(params.expressionType,''mRNA''))\n');
            fprintf(scriptFID,'\t\t\t\t\t\tp.reporterDegRates(1) = mRNADegMean;\n');
            fprintf(scriptFID,'\t\t\t\t\t\tp.reporterDegRates(2) = mRNADegSDev;\n');
            fprintf(scriptFID,'\t\t\t\t\telse\n');
            fprintf(scriptFID,'\t\t\t\t\t\tnatDegMean = mRNADegMean;\n');
            fprintf(scriptFID,'\t\t\t\t\t\tnatDegVar = mRNADegSDev * mRNADegSDev;\n');
            fprintf(scriptFID,'\t\t\t\t\tend\n');
            fprintf(scriptFID,'\t\t\t\tend\n');
            fprintf(scriptFID,'\n\t\t\t\tif strfind(params.expressionType,''protein'')\n');
            fprintf(scriptFID,'\t\t\t\t\tproteinDegMeanIdx = find(strcmpi(''protein deg rate mean'',importedDegRateNames(1,:)),1);\n');
            fprintf(scriptFID,'\t\t\t\t\tproteinDegSDevIdx = find(strcmpi(''protein deg rate s.dev'',importedDegRateNames(1,:)),1);\n');
            fprintf(scriptFID,'\t\t\t\t\tproteinDegMean = importedDegRates(degIdx-1,proteinDegMeanIdx-1);\n');
            fprintf(scriptFID,'\t\t\t\t\tproteinDegSDev = importedDegRates(degIdx-1,proteinDegSDevIdx-1);\n');
            fprintf(scriptFID,'\t\t\t\t\tif ~isempty(proteinDegMean) && ~isempty(proteinDegSDev)\n');
            fprintf(scriptFID,'\t\t\t\t\t\tp.reporterDegRates(3) = proteinDegMean;\n');
            fprintf(scriptFID,'\t\t\t\t\t\tp.reporterDegRates(4) = proteinDegSDev;\n');
            fprintf(scriptFID,'\t\t\t\t\tend\n');
            fprintf(scriptFID,'\t\t\t\tend\n');
            fprintf(scriptFID,'\t\t\tend\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\n\t\tfprintf(''%%d/%%d: "%%s" currently calculating (%%d replicates): %%s\\n'',x,length(uniqueNames),gene,p.maximumReplicates,datestr(now));\n');
            fprintf(scriptFID,'\t\tReTrOSsmooth_fits{x} = ReTrOSsmooth(p,data'',repmat(time,1,length(geneIdx))'',natDegMean,natDegVar,gene);\n');
            fprintf(scriptFID,'\t\tplotNames{x} = [''..'' filesep ''output'' filesep outputFile filesep ''plots_smooth'' filesep checkValidFileChar(gene)];\n');
            fprintf(scriptFID,'\t\tset(gcf,''PaperUnits'',''centimeters'',''PaperSize'',[29.7,21],''PaperPosition'',[0,0,29.7,21]);\n');
            fprintf(scriptFID,'\t\tprint(''-dpdf'',''-r150'',plotNames{x});\n');
            fprintf(scriptFID,'\t\tprint(''-dpng'',''-r150'',plotNames{x});\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'\tgeneNames = uniqueNames;\n');
            fprintf(scriptFID,'else\n');
            fprintf(scriptFID,'\tReTrOSsmooth_fits = cell(length(geneNames),1);\n');
            fprintf(scriptFID,'\tplotNames = cell(length(geneNames),1);\n');
            fprintf(scriptFID,'\tgeneNamesOrig = geneNames;\n');
            fprintf(scriptFID,'\tfor x = 1:length(geneNames)\n');
            fprintf(scriptFID,'\t\tp = params;\n');
            fprintf(scriptFID,'\t\tp.maximumReplicates = maxReps;\n');
            fprintf(scriptFID,'\n\t\tgene = geneNames{x};\n');
            fprintf(scriptFID,'\n\t\tnatDegMean = [];\n');
            fprintf(scriptFID,'\t\tnatDegVar = [];\n');
            fprintf(scriptFID,'\t\tif params.useDegRateFile\n');
            fprintf(scriptFID,'\t\t\tdegIdx = find(strcmpi(gene,importedDegRateNames(:,1)),1);\n');
            fprintf(scriptFID,'\t\t\tif ~isempty(degIdx)\n');
            fprintf(scriptFID,'\t\t\t\tmRNADegMeanIdx = find(strcmpi(''mRNA deg rate mean'',importedDegRateNames(1,:)),1);\n');
            fprintf(scriptFID,'\t\t\t\tmRNADegSDevIdx = find(strcmpi(''mRNA deg rate s.dev'',importedDegRateNames(1,:)),1);\n');
            fprintf(scriptFID,'\t\t\t\tmRNADegMean = importedDegRates(degIdx-1,mRNADegMeanIdx-1);\n');
            fprintf(scriptFID,'\t\t\t\tmRNADegSDev = importedDegRates(degIdx-1,mRNADegSDevIdx-1);\n');
            fprintf(scriptFID,'\t\t\t\tif ~isempty(mRNADegMean) && ~isempty(mRNADegSDev)\n');
            fprintf(scriptFID,'\t\t\t\t\tif ~isempty(strfind(params.expressionType,''mRNA''))\n');
            fprintf(scriptFID,'\t\t\t\t\t\tp.reporterDegRates(1) = mRNADegMean;\n');
            fprintf(scriptFID,'\t\t\t\t\t\tp.reporterDegRates(2) = mRNADegSDev;\n');
            fprintf(scriptFID,'\t\t\t\t\telse\n');
            fprintf(scriptFID,'\t\t\t\t\t\tnatDegMean = mRNADegMean;\n');
            fprintf(scriptFID,'\t\t\t\t\t\tnatDegVar = mRNADegSDev * mRNADegSDev;\n');
            fprintf(scriptFID,'\t\t\t\t\tend\n');
            fprintf(scriptFID,'\t\t\t\tend\n');
            fprintf(scriptFID,'\n\t\t\t\tif strfind(params.expressionType,''protein'')\n');
            fprintf(scriptFID,'\t\t\t\t\tproteinDegMeanIdx = find(strcmpi(''protein deg rate mean'',importedDegRateNames(1,:)),1);\n');
            fprintf(scriptFID,'\t\t\t\t\tproteinDegSDevIdx = find(strcmpi(''protein deg rate s.dev'',importedDegRateNames(1,:)),1);\n');
            fprintf(scriptFID,'\t\t\t\t\tproteinDegMean = importedDegRates(degIdx-1,proteinDegMeanIdx-1);\n');
            fprintf(scriptFID,'\t\t\t\t\tproteinDegSDev = importedDegRates(degIdx-1,proteinDegSDevIdx-1);\n');
            fprintf(scriptFID,'\t\t\t\t\tif ~isempty(proteinDegMean) && ~isempty(proteinDegSDev)\n');
            fprintf(scriptFID,'\t\t\t\t\t\tp.reporterDegRates(3) = proteinDegMean;\n');
            fprintf(scriptFID,'\t\t\t\t\t\tp.reporterDegRates(4) = proteinDegSDev;\n');
            fprintf(scriptFID,'\t\t\t\t\tend\n');
            fprintf(scriptFID,'\t\t\t\tend\n');
            fprintf(scriptFID,'\t\t\tend\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\n\t\tif sum(strcmpi(gene,geneNamesOrig)) > 1\n');
            fprintf(scriptFID,'\t\t\tgeneNum = sum(strcmpi(gene,geneNamesOrig(1:x)));\n');
            fprintf(scriptFID,'\t\t\tgene = [gene ''_'' num2str(geneNum)];\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\t\tgeneNames{x} = gene;\n');
            fprintf(scriptFID,'\n\t\tfprintf(''%%d/%%d: "%%s" currently calculating (%%d replicates): %%s\\n'',x,length(geneNames),gene,p.maximumReplicates,datestr(now));\n');
            fprintf(scriptFID,'\t\tReTrOSsmooth_fits{x} = ReTrOSsmooth(p,expressionData(x,:)'',time'',natDegMean,natDegVar,gene);\n');
            fprintf(scriptFID,'\t\tplotNames{x} = [''..'' filesep ''output'' filesep outputFile filesep ''plots_smooth'' filesep checkValidFileChar(gene)];\n');
            fprintf(scriptFID,'\t\tset(gcf,''PaperUnits'',''centimeters'',''PaperSize'',[29.7,21],''PaperPosition'',[0,0,29.7,21]);\n');
            fprintf(scriptFID,'\t\tprint(''-dpdf'',''-r150'',plotNames{x});\n');
            fprintf(scriptFID,'\t\tprint(''-dpng'',''-r150'',plotNames{x});\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'close(figH);\n');
            fprintf(scriptFID,'\ndisp(''Creating output'');\n');
            fprintf(scriptFID,'save([''..'' filesep ''output'' filesep outputFile filesep outputFile ''_smooth.mat''],''ReTrOSsmooth_fits'',''geneNames'',''inputFile'');\n');
            fprintf(scriptFID,'\nif ~ispc\n');
            fprintf(scriptFID,'\tpdfFile = [''..'' filesep ''output'' filesep outputFile filesep outputFile ''_smooth.pdf''];\n');
            fprintf(scriptFID,'\tparamsFile = [''..'' filesep ''output'' filesep outputFile filesep ''params.txt''];\n');
            fprintf(scriptFID,'\tparamsFid = fopen(paramsFile, ''w'');\n');
            fprintf(scriptFID,'\tfprintf(paramsFid,''-dNOPAUSE\\n'');\n');
            fprintf(scriptFID,'\tfprintf(paramsFid,''-dBATCH\\n'');\n');
            fprintf(scriptFID,'\tfprintf(paramsFid,''-sOutputFile="%%s"\\n'', pdfFile);\n');
            fprintf(scriptFID,'\tfprintf(paramsFid,''-sDEVICE=%%s\\n'', ''pdfwrite'');\n');
            fprintf(scriptFID,'\tfor m = 1:length(geneNames)\n');
            fprintf(scriptFID,'\t\tfprintf(paramsFid,[''"'' plotNames{m} ''.pdf"\\n'']);\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'\tfclose(paramsFid);\n');
            fprintf(scriptFID,'\t[status result] = system([''/usr/local/bin/gs @"'' paramsFile ''"'']);\n');
            fprintf(scriptFID,'\tdelete(paramsFile);\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'\ndisp([''Finished ('' datestr(clock) '')'']);\n');
            fprintf(scriptFID,'\nrmpath([''..'' filesep ''smooth'']);\n');
            fprintf(scriptFID,'warning(currentWarnings);\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'\nfunction o = checkValidFileChar(i)\n');
            fprintf(scriptFID,'invalidChar = ''[\\/:*?"<>|;'''']*'';\n');
            fprintf(scriptFID,'o = regexprep(i,invalidChar,''_'');\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'\nfunction u = uniqueRetainOrder(i)\n');
            fprintf(scriptFID,'u = {i{1}};\n');
            fprintf(scriptFID,'for j = 2:length(i)\n');
            fprintf(scriptFID,'\tif sum(strcmpi(u,i{j})) == 0\n');
            fprintf(scriptFID,'\t\tu = [u ; i{j}];\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'end\n');
            fclose(scriptFID);
        end
        
        function [] = createSwitchScript
            if exist(dataStruct.switchScriptFile,'file')
                disp(['ReTrOS-switch: file "' dataStruct.switchScriptFile '" already exists - overwriting']);
            end
            scriptFID = fopen(dataStruct.switchScriptFile,'w');
            fprintf(scriptFID,['function [ReTrOSswitch_fits] = run_' dataStruct.outputName '_switch\n']);
            fprintf(scriptFID,['\n%%Generated ReTrOS-switch script (' datestr(clock) ')\n']);
            fprintf(scriptFID,'\ncurrentWarnings = warning(''off'',''all'');\n');
            fprintf(scriptFID,'addpath([''..'' filesep ''switch'']);\n');
            fprintf(scriptFID,'\n%%Parameters:\n');
            fprintf(scriptFID,['inputFile = ''' dataStruct.file ''';\n']);
            fprintf(scriptFID,['outputFile = ''' dataStruct.outputName ''';\n']);
            fprintf(scriptFID,['combineReplicates = ' logical2str(dataStruct.switchCombineReplicates) ';\n']);
            fprintf(scriptFID,['nameColumnIdx = ' num2str(dataStruct.nameColumnSelected) ';\n']);
            fprintf(scriptFID,['expressionType = ''' dataStruct.dataType ''';\n']);
            fprintf(scriptFID,['dataDetrend = ''' dataStruct.dataDetrend ''';\n']);
            fprintf(scriptFID,['regressionMethod = ''' dataStruct.switchRegressionMethod ''';\n']);
            fprintf(scriptFID,['mRNADegRates = [' num2str(dataStruct.switchmRNADegRates(1)) ' ' num2str(dataStruct.switchmRNADegRates(2)) '];\n']);
            fprintf(scriptFID,['proteinDegRates = [' num2str(dataStruct.switchProteinDegRates(1)) ' ' num2str(dataStruct.switchProteinDegRates(2)) '];\n']);
            fprintf(scriptFID,['MCMCiterations = ' num2str(dataStruct.switchIterations) ';\n']);
            fprintf(scriptFID,['initialSwitchGuess = ' logical2str(dataStruct.switchIntialSwitchGuess) ';\n']);
            fprintf(scriptFID,['maximumSwitchNum = ' num2str(dataStruct.switchMaxSwitchNum) ';\n']);
            fprintf(scriptFID,['expectedSwitchNum = ' num2str(dataStruct.switchExpSwitchNum) ';\n']);
            fprintf(scriptFID,['switchDelayTime = ' num2str(dataStruct.switchDelayTime) ';\n']);
            fprintf(scriptFID,['burnIn = ' num2str(dataStruct.switchBurnIn) ';\n']);
            fprintf(scriptFID,['numSamplesEstPost = ' num2str(dataStruct.switchNumSamplesEstPost) ';\n']);
            fprintf(scriptFID,['baselineFactor = ' num2str(dataStruct.switchBaselineFactor) ';\n']);
            fprintf(scriptFID,['maxSwitchDeviationFactor = ' num2str(dataStruct.switchMaxSwDev) ';\n']);
            fprintf(scriptFID,['minSwStr = ' num2str(dataStruct.switchMinSwStr) ';\n']);
            fprintf(scriptFID,['numSamplePlot = ' num2str(dataStruct.switchNumSamplePlot) ';\n']);
            fprintf(scriptFID,'\nmkdir([''..'' filesep ''output'']);\n');
            fprintf(scriptFID,'mkdir([''..'' filesep ''output''],outputFile);\n');
            fprintf(scriptFID,'mkdir([''..'' filesep ''output'' filesep outputFile],''plots_switch'');\n');
            fprintf(scriptFID,'\ncontents = importdata(inputFile,''\\t'');\n');
            fprintf(scriptFID,'geneNames = contents.textdata(2:end,nameColumnIdx);\n');
            fprintf(scriptFID,'time = contents.data(1,:);\n');
            fprintf(scriptFID,'expressionData = contents.data(2:end,:);\n');
            fprintf(scriptFID,'timescale = unique(time);\n');
            fprintf(scriptFID,'\nmaxReps = 0;\n');
            fprintf(scriptFID,'for x = 1:length(timescale)\n');
            fprintf(scriptFID,'\tif sum(timescale(x) == time) > maxReps\n');
            fprintf(scriptFID,'\t\tmaxReps = sum(timescale(x) == time);\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'figH = figure(''name'',[''ReTrOS-switch: '' outputFile],''numbertitle'',''off'');\n');
            fprintf(scriptFID,'%%Main loop\n');
            fprintf(scriptFID,'if combineReplicates\n');
            fprintf(scriptFID,'\tuniqueNames = uniqueRetainOrder(geneNames);\n');
            fprintf(scriptFID,'\tNoG = length(uniqueNames);\n');
            fprintf(scriptFID,'\tReTrOSswitch_fits = cell(NoG,1);\n');
            fprintf(scriptFID,'\tplotNames = cell(length(uniqueNames)*2,1);\n');
            fprintf(scriptFID,'\tfor x = 1:NoG\n');
            fprintf(scriptFID,'\t\tgene = uniqueNames{x};\n');
            fprintf(scriptFID,'\t\tgeneIdx = find(strcmpi(gene,geneNames));\n');
            fprintf(scriptFID,'\t\tdata = nan(1,length(timescale) * maxReps * length(geneIdx));\n');
            fprintf(scriptFID,'\t\tfor y = 1:length(geneIdx)\n');
            fprintf(scriptFID,'\t\t\toffset = (y-1) * maxReps * length(timescale);\n');
            fprintf(scriptFID,'\t\t\tfor z = 1:length(timescale)\n');
            fprintf(scriptFID,'\t\t\t\td = expressionData(geneIdx(y),timescale(z)==time);\n');
            fprintf(scriptFID,'\t\t\t\tfor w = 1:length(d)\n');
            fprintf(scriptFID,'\t\t\t\t\tdata(offset + (w-1) * length(timescale) + z) = d(w);\n');
            fprintf(scriptFID,'\t\t\t\tend\n');
            fprintf(scriptFID,'\t\t\tend\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\t\tfprintf(''%%d/%%d: "%%s" currently calculating (%%d replicates): %%s\\n'',x,NoG,gene,maxReps*length(geneIdx),datestr(now));\n');
            fprintf(scriptFID,'\t\ttry\n');
            fprintf(scriptFID,'\t\t\tif strfind(expressionType,''mRNA'')\n');
            fprintf(scriptFID,'\t\t\t\t[fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_mRNA(gene,data,timescale,maxReps*length(geneIdx),MCMCiterations,''detrendData'',dataDetrend,''regressionMethod'',regressionMethod,''mRNADegradationRatePrior'',mRNADegRates,''useInitialSwitchSuggest'',initialSwitchGuess,''maxSwitches'',maximumSwitchNum,''switchNumberPrior'',expectedSwitchNum,''switchDelay'',switchDelayTime,''burnIn'',burnIn,''maxSamplesToUse'',numSamplesEstPost,''switchBaselineFactor'',baselineFactor,''minSwitchStrength'',minSwStr,''maxSwitchDeviationFactor'',maxSwitchDeviationFactor,''maxSamplesToPlot'',numSamplePlot);\n');
            fprintf(scriptFID,'\t\t\telse\n');
            fprintf(scriptFID,'\t\t\t\t[fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_protein(gene,data,timescale,maxReps*length(geneIdx),MCMCiterations,''detrendData'',dataDetrend,''regressionMethod'',regressionMethod,''mRNADegradationRatePrior'',mRNADegRates,''proteinDegradationRatePrior'',proteinDegRates,''useInitialSwitchSuggest'',initialSwitchGuess,''maxSwitches'',maximumSwitchNum,''switchNumberPrior'',expectedSwitchNum,''switchDelay'',switchDelayTime,''burnIn'',burnIn,''maxSamplesToUse'',numSamplesEstPost,''switchBaselineFactor'',baselineFactor,''minSwitchStrength'',minSwStr,''maxSwitchDeviationFactor'',maxSwitchDeviationFactor,''maxSamplesToPlot'',numSamplePlot);\n');
            fprintf(scriptFID,'\t\t\tend\n');
            fprintf(scriptFID,'\t\t\tReTrOSswitch_fits{x,1} = gene;\n');
            fprintf(scriptFID,'\t\t\tReTrOSswitch_fits{x,2} = fitData;\n');
            fprintf(scriptFID,'\t\t\tReTrOSswitch_fits{x,3} = fitSamples;\n');
            fprintf(scriptFID,'\t\t\tReTrOSswitch_fits{x,4} = fitPosteriors;\n');
            fprintf(scriptFID,'\t\t\tReTrOSswitch_fits{x,5} = fitParams;\n');
            fprintf(scriptFID,'\t\t\tplotNames{x*2-1} = [''..'' filesep ''output'' filesep outputFile filesep ''plots_switch'' filesep checkValidFileChar(gene)];\n');
            fprintf(scriptFID,'\t\t\tset(gcf,''PaperUnits'',''centimeters'',''PaperSize'',[29.7,21],''PaperPosition'',[0,0,29.7,21]);\n');
            fprintf(scriptFID,'\t\t\tprint(''-dpdf'',''-r150'',plotNames{x*2-1});\n');
            fprintf(scriptFID,'\t\t\tprint(''-dpng'',''-r150'',plotNames{x*2-1});\n');
            fprintf(scriptFID,'\t\t\tplotSwitchModelsBasic(fitPosteriors,fitSamples,fitParams,fitData);\n');
            fprintf(scriptFID,'\t\t\tplotNames{x*2} = [''..'' filesep ''output'' filesep outputFile filesep ''plots_switch'' filesep checkValidFileChar(gene) ''_models''];\n');
            fprintf(scriptFID,'\t\t\tset(gcf,''PaperUnits'',''centimeters'',''PaperSize'',[29.7,21],''PaperPosition'',[0,0,29.7,21]);\n');
            fprintf(scriptFID,'\t\t\tprint(''-dpdf'',''-r150'',plotNames{x*2});\n');
            fprintf(scriptFID,'\t\t\tprint(''-dpng'',''-r150'',plotNames{x*2});\n');
            fprintf(scriptFID,'\t\tcatch exception\n');
            fprintf(scriptFID,'\t\t\tdisp(''Error during calculation:'');\n');
            fprintf(scriptFID,'\t\t\tdisp(getReport(exception));\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'\tgeneNames = uniqueNames;\n');
            fprintf(scriptFID,'else\n');
            fprintf(scriptFID,'\tNoG = size(expressionData,1);\n');
            fprintf(scriptFID,'\tReTrOSswitch_fits = cell(NoG,5);\n');
            fprintf(scriptFID,'\tplotNames = cell(NoG*2,1);\n');
            fprintf(scriptFID,'\tdata = nan(NoG,length(timescale * maxReps));\n');
            fprintf(scriptFID,'\tfor x = 1:length(timescale)\n');
            fprintf(scriptFID,'\t\tt = find(time == timescale(x));\n');
            fprintf(scriptFID,'\t\tfor y = 1:length(t)\n');
            fprintf(scriptFID,'\t\t\tdata(:,(y-1) * length(timescale) + x) = expressionData(:,t(y));\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'\tgeneNamesOrig = geneNames;\n');
            fprintf(scriptFID,'\tfor x = 1:NoG\n');
            fprintf(scriptFID,'\t\tgene = geneNames{x};\n');
            fprintf(scriptFID,'\t\tif sum(strcmpi(geneNamesOrig,gene)) > 1\n');
            fprintf(scriptFID,'\t\t\tgeneNum = sum(strcmpi(gene,geneNamesOrig(1:x)));\n');
            fprintf(scriptFID,'\t\t\tgene = [gene ''_'' num2str(geneNum)];\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\t\tgeneNames{x} = gene;\n');
            fprintf(scriptFID,'\t\tfprintf(''%%d/%%d: "%%s" currently calculating (%%d replicates): %%s\\n'',x,NoG,gene,maxReps,datestr(now));\n');
            fprintf(scriptFID,'\t\ttry\n');
            fprintf(scriptFID,'\t\t\tif strfind(expressionType,''mRNA'')\n');
            fprintf(scriptFID,'\t\t\t\t[fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_mRNA(gene,data(x,:),timescale,maxReps,MCMCiterations,''detrendData'',dataDetrend,''regressionMethod'',regressionMethod,''mRNADegradationRatePrior'',mRNADegRates,''useInitialSwitchSuggest'',initialSwitchGuess,''maxSwitches'',maximumSwitchNum,''switchNumberPrior'',expectedSwitchNum,''switchDelay'',switchDelayTime,''burnIn'',burnIn,''maxSamplesToUse'',numSamplesEstPost,''switchBaselineFactor'',baselineFactor,''minSwitchStrength'',minSwStr,''maxSwitchDeviationFactor'',maxSwitchDeviationFactor,''maxSamplesToPlot'',numSamplePlot);\n');
            fprintf(scriptFID,'\t\t\telse\n');
            fprintf(scriptFID,'\t\t\t\t[fitData, fitSamples, fitPosteriors, fitParams] = ReTrOSswitch_protein(gene,data(x,:),timescale,maxReps,MCMCiterations,''detrendData'',dataDetrend,''regressionMethod'',regressionMethod,''mRNADegradationRatePrior'',mRNADegRates,''proteinDegradationRatePrior'',proteinDegRates,''useInitialSwitchSuggest'',initialSwitchGuess,''maxSwitches'',maximumSwitchNum,''switchNumberPrior'',expectedSwitchNum,''switchDelay'',switchDelayTime,''burnIn'',burnIn,''maxSamplesToUse'',numSamplesEstPost,''switchBaselineFactor'',baselineFactor,''minSwitchStrength'',minSwStr,''maxSwitchDeviationFactor'',maxSwitchDeviationFactor,''maxSamplesToPlot'',numSamplePlot);\n');
            fprintf(scriptFID,'\t\t\tend\n');
            fprintf(scriptFID,'\t\t\tReTrOSswitch_fits{x,1} = gene;\n');
            fprintf(scriptFID,'\t\t\tReTrOSswitch_fits{x,2} = fitData;\n');
            fprintf(scriptFID,'\t\t\tReTrOSswitch_fits{x,3} = fitSamples;\n');
            fprintf(scriptFID,'\t\t\tReTrOSswitch_fits{x,4} = fitPosteriors;\n');
            fprintf(scriptFID,'\t\t\tReTrOSswitch_fits{x,5} = fitParams;\n');
            fprintf(scriptFID,'\t\t\tplotNames{x*2-1} = [''..'' filesep ''output'' filesep outputFile filesep ''plots_switch'' filesep checkValidFileChar(gene)];\n');
            fprintf(scriptFID,'\t\t\tset(gcf,''PaperUnits'',''centimeters'',''PaperSize'',[29.7,21],''PaperPosition'',[0,0,29.7,21]);\n');
            fprintf(scriptFID,'\t\t\tprint(''-dpdf'',''-r150'',plotNames{x*2-1});\n');
            fprintf(scriptFID,'\t\t\tprint(''-dpng'',''-r150'',plotNames{x*2-1});\n');
            fprintf(scriptFID,'\t\t\tplotSwitchModelsBasic(fitPosteriors,fitSamples,fitParams,fitData);\n');
            fprintf(scriptFID,'\t\t\tplotNames{x*2} = [''..'' filesep ''output'' filesep outputFile filesep ''plots_switch'' filesep checkValidFileChar(gene) ''_models''];\n');
            fprintf(scriptFID,'\t\t\tset(gcf,''PaperUnits'',''centimeters'',''PaperSize'',[29.7,21],''PaperPosition'',[0,0,29.7,21]);\n');
            fprintf(scriptFID,'\t\t\tprint(''-dpdf'',''-r150'',plotNames{x*2});\n');
            fprintf(scriptFID,'\t\t\tprint(''-dpng'',''-r150'',plotNames{x*2});\n');
            fprintf(scriptFID,'\t\tcatch exception\n');
            fprintf(scriptFID,'\t\t\tdisp(''Error during calculation:'');\n');
            fprintf(scriptFID,'\t\t\tdisp(getReport(exception));\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'disp(''Creating output'');\n');
            fprintf(scriptFID,'%%output switch times and switch times plot\n');
            fprintf(scriptFID,'swTimesFile = [''..'' filesep ''output'' filesep outputFile filesep outputFile ''_switchTimes.txt''];\n');
            fprintf(scriptFID,'swTimesFID = fopen(swTimesFile, ''w'');\n');
            fprintf(scriptFID,'fprintf(swTimesFID,''gene\\tswitch type\\tswitch mean time\\tswitch standard deviation\\tswitch strength\\n'');\n');
            fprintf(scriptFID,'for x = 1:length(geneNames)\n');
            fprintf(scriptFID,'\tposteriors = ReTrOSswitch_fits{x,4};\n');
            fprintf(scriptFID,'\tgene = ReTrOSswitch_fits{x,1};\n');
            fprintf(scriptFID,'\tif ~isempty(posteriors)\n');
            fprintf(scriptFID,'\t\tswTimes = posteriors.switchDistribution.means;\n');
            fprintf(scriptFID,'\t\tswSigmas = posteriors.switchDistribution.standardDeviations;\n');
            fprintf(scriptFID,'\t\tswStrengths = posteriors.switchDistribution.strengths;\n');
            fprintf(scriptFID,'\t\trates = posteriors.modelFit.taus_medianModel;\n');
            fprintf(scriptFID,'\t\tfor z = 1:length(swTimes)\n');
            fprintf(scriptFID,'\t\t\tif (rates(z+1)-rates(z)) > 0\n');
            fprintf(scriptFID,'\t\t\t\ttype = ''increase'';\n');
            fprintf(scriptFID,'\t\t\telse\n');
            fprintf(scriptFID,'\t\t\t\ttype = ''decrease'';\n');
            fprintf(scriptFID,'\t\t\tend\n');
            fprintf(scriptFID,'\t\t\tfprintf(swTimesFID,[gene ''\\t'' type ''\\t'' num2str(swTimes(z)) ''\\t'' num2str(swSigmas(z)) ''\\t'' num2str(swStrengths(z)) ''\\n'']);\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'fclose(swTimesFID);\n');
            fprintf(scriptFID,'plotNum = length(geneNames)*2;\n');
            fprintf(scriptFID,'outputSwPlot = plotSwitchHeatmap(swTimesFile,timescale);\n');
            fprintf(scriptFID,'if ~isempty(outputSwPlot)\n');
            fprintf(scriptFID,'\tplotNum = plotNum+1;\n');
            fprintf(scriptFID,'\tplotNames{plotNum} = [''..'' filesep ''output'' filesep outputFile filesep ''plots_switch'' filesep ''switchTimes''];\n');
            fprintf(scriptFID,'\tset(gcf,''PaperUnits'',''centimeters'',''PaperSize'',[29.7,21],''PaperPosition'',[0,0,29.7,21]);\n');
            fprintf(scriptFID,'\tprint(''-dpdf'',''-r150'',plotNames{plotNum});\n');
            fprintf(scriptFID,'\tprint(''-dpng'',''-r150'',plotNames{plotNum});\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'close(figH);\n');
            fprintf(scriptFID,'%%output degradation rates\n');
            fprintf(scriptFID,'degRatesFile = [''..'' filesep ''output'' filesep outputFile filesep outputFile ''_degRates.txt''];\n');
            fprintf(scriptFID,'degRatesFID = fopen(degRatesFile,''w'');\n');
            fprintf(scriptFID,'if strfind(expressionType,''mRNA'')\n');
            fprintf(scriptFID,'\tfprintf(degRatesFID,''gene\\tmRNA deg rate mean\\tmRNA deg rate s.dev\\n'');\n');
            fprintf(scriptFID,'\tfor x = 1:length(geneNames)\n');
            fprintf(scriptFID,'\t\tposteriors = ReTrOSswitch_fits{x,4};\n');
            fprintf(scriptFID,'\t\tgene = ReTrOSswitch_fits{x,1};\n');
            fprintf(scriptFID,'\t\tif ~isempty(posteriors)\n');
            fprintf(scriptFID,'\t\t\tfprintf(degRatesFID,[gene ''\\t'' num2str(posteriors.degradationRate_mRNA.mean) ''\\t'' num2str(posteriors.degradationRate_mRNA.standardDeviation) ''\\n'']);\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'else\n');
            fprintf(scriptFID,'\tfprintf(degRatesFID,''gene\\tmRNA deg rate mean\\tmRNA deg rate s.dev\\tprotein deg rate mean\\tprotein deg rate s.dev\\n'');\n');
            fprintf(scriptFID,'\tfor x = 1:length(geneNames)\n');
            fprintf(scriptFID,'\t\tposteriors = ReTrOSswitch_fits{x,4};\n');
            fprintf(scriptFID,'\t\tgene = ReTrOSswitch_fits{x,1};\n');
            fprintf(scriptFID,'\t\tif ~isempty(posteriors)\n');
            fprintf(scriptFID,'\t\t\tfprintf(degRatesFID,[gene ''\\t'' num2str(posteriors.degradationRate_mRNA.mean) ''\\t'' num2str(posteriors.degradationRate_mRNA.standardDeviation) ''\\t'' num2str(posteriors.degradationRate_protein.mean) ''\\t'' num2str(posteriors.degradationRate_protein.standardDeviation) ''\\n'']);\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'fclose(degRatesFID);\n');
            fprintf(scriptFID,'save([''..'' filesep ''output'' filesep outputFile filesep outputFile ''_switch.mat''],''inputFile'',''geneNames'',''ReTrOSswitch_fits'');\n');
            fprintf(scriptFID,'if ~ispc\n');
            fprintf(scriptFID,'\tpdfFile = [''..'' filesep ''output'' filesep outputFile filesep outputFile ''_switch.pdf''];\n');
            fprintf(scriptFID,'\tparamsFile = [''..'' filesep ''output'' filesep outputFile filesep ''params.txt''];\n');
            fprintf(scriptFID,'\tparamsFid = fopen(paramsFile, ''w'');\n');
            fprintf(scriptFID,'\tfprintf(paramsFid,''-dNOPAUSE\\n'');\n');
            fprintf(scriptFID,'\tfprintf(paramsFid,''-dBATCH\\n'');\n');
            fprintf(scriptFID,'\tfprintf(paramsFid,''-sOutputFile="%%s"\\n'', pdfFile);\n');
            fprintf(scriptFID,'\tfprintf(paramsFid,''-sDEVICE=%%s\\n'', ''pdfwrite'');\n');
            fprintf(scriptFID,'\tfor m = 1:plotNum\n');
            fprintf(scriptFID,'\t\tif ~isempty(plotNames{m})\n');
            fprintf(scriptFID,'\t\t\tfprintf(paramsFid,[''"'' plotNames{m} ''.pdf"\\n'']);\n');
            fprintf(scriptFID,'\t\tend\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'\tfclose(paramsFid);\n');
            fprintf(scriptFID,'\t[status, result] = system([''/usr/local/bin/gs @"'' paramsFile ''"'']);\n');
            fprintf(scriptFID,'\tdelete(paramsFile);\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'disp([''Finished ('' datestr(clock) '')'']);\n');
            fprintf(scriptFID,'rmpath([''..'' filesep ''switch'']);\n');
            fprintf(scriptFID,'warning(currentWarnings);\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'\nfunction o = checkValidFileChar(i)\n');
            fprintf(scriptFID,'invalidChar = ''[\\/:*?"<>|;'''']*'';\n');
            fprintf(scriptFID,'o = regexprep(i,invalidChar,''_'');\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'\nfunction u = uniqueRetainOrder(i)\n');
            fprintf(scriptFID,'u = {i{1}};\n');
            fprintf(scriptFID,'for j = 2:length(i)\n');
            fprintf(scriptFID,'\tif sum(strcmpi(u,i{j})) == 0\n');
            fprintf(scriptFID,'\t\tu = [u ; i{j}];\n');
            fprintf(scriptFID,'\tend\n');
            fprintf(scriptFID,'end\n');
            fprintf(scriptFID,'end\n');
            fclose(scriptFID);
        end
        
        function [] = runButtonEvent(varargin)
            if ~RETROS_RUNNING
                if createRunScriptsCheckBox.isSelected
                    disp('SCRIPTS NOT CORRECT');
                    mkdir('scripts');
                    createSmoothScript;
                    createSwitchScript;
                end
                javaMethodEDT('setEnabled',runButton,false);
                javaMethodEDT('setEnabled',backButton,false);
                javaMethodEDT('setText',cancelButton,'Cancel');
                if runSmoothCheckBox.isSelected
                    if createRunScriptsCheckBox.isSelected
                        scriptFileName = dataStruct.smoothScriptFile;
                        disp(['Running ReTrOS-smooth @ ',datestr(clock),' (' scriptFileName ')']);
                    else
                        disp(['Running ReTrOS-smooth @ ',datestr(clock)]);
                    end
                    %                     clear(scriptFileName);
                    %                     run(scriptFileName);
                    runReTrOSsmooth;
                end
                RETROS_CANCELLED = false;
                if runSwitchCheckBox.isSelected
                    if createRunScriptsCheckBox.isSelected
                        scriptFileName = dataStruct.switchScriptFile;
                        disp(['Running ReTrOS-switch @ ',datestr(clock),' (' scriptFileName ')']);
                    else
                        disp(['Running ReTrOS-switch @ ',datestr(clock)]);
                    end
                    %                     clear(scriptFileName);
                    %                     run(scriptFileName);
                    runReTrOSswitch;
                end
                RETROS_CANCELLED = false;
                javaMethodEDT('setEnabled',runButton,true);
                javaMethodEDT('setEnabled',backButton,true);
                javaMethodEDT('setText',cancelButton,'Close');
            end
        end
        
        function [] = backButtonEvent(varargin)
            if ~RETROS_RUNNING
                algorithmSelectUI(frame);
            end
        end
        
        function [] = cancelButtonEvent(varargin)
            %         java.system.GC;
            if RETROS_RUNNING
                RETROS_CANCELLED = true;
                disp('Stopping current ReTrOS execution');
                RETROS_RUNNING = false;
            elseif ~RETROS_CANCELLED
                dataStruct = [];
                acceptedInputs = false;
                javaMethodEDT('dispose',frame);
                close(hidFig);
            end
        end
        
        function [] = updateRandomSeed(varargin)
            if randomSeedCheckBox.isSelected
                javaMethodEDT('setEnabled',randomSeedField,false);
            else
                javaMethodEDT('setEnabled',randomSeedField,true);
            end
        end
        
        function [] = updateParallel(varargin)
            if runParallelCheckBox.isSelected
                javaMethodEDT('setEnabled',workersComboBox,true);
            else
                javaMethodEDT('setEnabled',workersComboBox,false);
                javaMethodEDT('setSelectedIndex',workersComboBox,0);
            end
        end
    end
end


function o = checkValidFileChar(i)
invalidChar = '[\/:*?"<>|;'']*';
o = regexprep(i,invalidChar,'_');
end

function [s] = logical2str(a)
if a
    s = 'true';
else
    s = 'false';
end
end

function u = uniqueRetainOrder(i)
u = {i{1}};
for j = 2:length(i)
    if sum(strcmpi(u,i{j})) == 0
        u = [u ; i{j}];
    end
end
end