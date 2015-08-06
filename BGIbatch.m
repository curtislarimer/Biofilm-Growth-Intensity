% BGI: measures Biofilm growth intensity from a group of images in a
% directory
% Written by Curtis Larimer
function BGIbatch
close all
% clear all

fName0 = uigetdir();
cd(fName0);

fileList = {};
contents = dir('*.JPG');
for i = 1:numel(contents)
    fileList{end+1,1} = contents(i).name;
end

if numel(fileList)==0
    h = errordlg('No files found in directory','','modal');
    uiwait(h)
end

newfile = questdlg('Create a new data file or save to an existing file?','Coverage Map.','New file','Choose existing','New file');
    if isempty(newfile) || strcmpi(newfile,'New file')
        DataFile = 'Biofilm Image Data.xls'; 
    else
        DataFile = uigetfile('*.xls','Select data file');
    end
    
Samplenumcolm = xlsread(DataFile,'A:A');
    if isempty(Samplenumcolm) == 1
        SampleNum = 0;
    else
        SampleNum = max(Samplenumcolm);
    end

coltitles = {'Sample','File','Control File','BGI gray','BGI red','BGI green','BGI blue','Select X1','Select Y1','Select X2','Select Y2','Control X1','Control Y1','Control X2','Control Y2'};
    sheet = 1;
    xlRange1 = 'A1';
    xlswrite(DataFile,coltitles,sheet,xlRange1);

% Start of image analysis
for ii = 1:length(fileList)

Pathname = fName0;
fName = char(fileList(ii));
Fullpath = strcat(Pathname,fName);
cd(Pathname);
BGIcallable(fName,DataFile,Fullpath)

end

end



