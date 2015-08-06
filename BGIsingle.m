% BGI: measures Biofilm growth intensity from a single image
% Written by Curtis Larimer for Pacific Northwest National Laboratory
% Direct questions to curtis.larimer@pnnl.gov

function BGIsingle
close all
clear all

newfile = questdlg('Create a new data file or save to an existing file?','Coverage Map.','New file','Choose existing','New file');
    if isempty(newfile) || strcmpi(newfile,'New file')
        DataFile = 'Biofilm Image Data.xls'; 
    else
        DataFile = uigetfile('*.xls','Select data file');
    end

[fName, Pathname] = uigetfile('*.jpg','Select a biofilm image');
Fullpath = strcat(Pathname,fName);
cd(Pathname);

BGIcallable(fName,DataFile,Fullpath)

end



