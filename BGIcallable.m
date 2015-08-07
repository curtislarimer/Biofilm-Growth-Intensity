% BGI: to be called from BGIbatch.m or BGIsingle.m for each image in a directory

% Written by Curtis Larimer for Pacific Northwest National Laboratory
% Direct questions to curtis.larimer@pnnl.gov


function BGIcallable(fName,DataFile,Fullpath)

% Read data file to identify sample number
Samplenumcolm = xlsread(DataFile,'A:A');
if isempty(Samplenumcolm) == 1
    SampleNum = 0;
else
SampleNum = max(Samplenumcolm);
end

% Select a control image
fName2 = uigetfile('*.jpg','Select a control image');
if ~ischar(fName)
    return;
end

Figure2 = figure(1); %creates the figure that the raw image is added to
biofilm = imread(fName);
set(Figure2,'units','normalized');
    figure(1)
    imshow(biofilm)
    hold on
    title('Select corners of  biofilm region')
    q = [];
%     Begin selecting corners of biofilm region (two clicks for upper left and lower right corners. Area outside this rectangle is not analyzed
    while isempty(q) || strcmpi(q,'no')
        corners1 = round(ginput(2));
        corners(1,1)=corners1(1,1);
        corners(4,1)=corners1(1,1);
        corners(1,2)= corners1(1,2);
        corners(2,2)= corners1(1,2);
        corners(2,1)= corners1(2,1);
        corners(3,1)= corners1(2,1);
        corners(3,2)= corners1(2,2);
        corners(4,2)= corners1(2,2);

        qq = line([corners(:,1);corners(1,1)],[corners(:,2);corners(1,2)],'marker','o','color','g'); %adds a colored box around the user selected region

            q = questdlg('Acceptable?','Line finished.','Yes','No','Yes'); %checks if selection is OK
        if isempty(q) || strcmpi(q,'no') %deletes selection if not ok and returns to reselect
            delete(qq)
        end
    end
%     End of selecting corners of biofilm region

    
%     Begin selecting control region - blank
Figure2 = figure(2); %creates the figure that the control image is added to
set(Figure2,'units','normalized');
    control = imread(fName2);
    figure(2)
    imshow(control)
    title('Select corners of  control region')
    hold on
    q = [];
    
    while isempty(q) || strcmpi(q,'no')
        blank = round(ginput(2));
        qq = line([blank(:,1); blank(2,1); blank(1,1); blank(1,1)],[blank(1,2); blank(:,2); blank(2,2); blank(1,2)],'marker','o','color','r');
        q = questdlg('Acceptable?','Line finished.','Yes','No','Yes');
        if isempty(q) || strcmpi(q,'no')
            delete(qq)
        end
    end

%     End selecting control region
    
    % separate data from biofilm image in selected area
    subSet = [min(corners) max(corners)];
    subSet = round(subSet([2 4 1 3]));
    plate = biofilm(subSet(1):subSet(2),subSet(3):subSet(4),:); %this is all the data inside the selection
    hold off
    
        %converts image to grayscale
    X = rgb2gray(plate);
    
    subSet = [min(blank) max(blank)]; %creates a set with the range of grayscales in the blank control
    bPlate = control(subSet(2):subSet(4),subSet(1):subSet(3),:); %creates an image of just selected area

    [pathstr, imagename, ext]=fileparts(Fullpath); %gets info about the full file path so transformed image can be saved to the same directory

    % ---- Begin BGI analysis of Separated RGB channels
    redplate = plate(:,:,1);
    greenplate = plate(:,:,2);
    blueplate = plate(:,:,3);

% Data from controls for each color channel	
    XBred = bPlate(:,:,1);
    XBgreen = bPlate(:,:,2);
    XBblue = bPlate(:,:,3);

% redplate BGI
N = 10; %Sets the number of levels in the multilevel threshold
[thresh, metric] = multithresh(redplate, N); %returns thresh a 1-by-N vector containing N threshold values calculated using Otsu's thresholding method
    biomax = double(max(max(redplate))); %determines the maximum level in the biofilm image
    biomin = double(min(min(redplate))); %determines the min level in the biofilm image
    biorange = double(biomax - biomin); % range of levels in biofilm image
    bioave = mean(mean(redplate)); % average of levels in biofilm image
    controlmax = double(max(max(XBred))); %determines the max level in the control image
    controlmin = double(min(min(XBred))); %determines the min level in the control
    controlrange = double(controlmax - controlmin); %range of control
    controlave = mean(mean(XBred)); %average of control
    
    Foulcontrast = 'some'; %if this is 'some' then the analysis can proceed normally 
    if abs(bioave-controlave) < controlrange/4 %tests if the biofilm image is very close in average levels to the control
        disp('fouled sample is very close to control sample - red')
        Foulcontrast = 'none'; 
	%Alternate analysis when sample is very similar to control
        redplate = redplate - biomin;
        threshint = 255/10;
        thresh = [threshint 2*threshint 3*threshint 4*threshint 5*threshint 6*threshint 7*threshint 8*threshint 9*threshint];
    end
    
    %normal BGI analysis when image contrast is sufficient

    Multithreshimagered = imquantize(redplate,thresh); %creates a new multilevel thresholded image

    GrayscaleImagered = abs(Multithreshimagered/(N+1)-1); %converts 0-9 scale to 0 to 1 scale

    uniqLevels = unique(GrayscaleImagered(:));
[a,b]=hist(GrayscaleImagered(:),uniqLevels); %creates a histogram with N levels

    a = fliplr(a); %flips histogram data left to right

    Totalbins = sum(a); %total of all counts in the histogram. varies based on size of image

    clrmapred = zeros(10,3); %Create var for redscale colormap
    for j = 1:length(a)
        BGIcalc(1,j) = a(1,j)*(j-1)*.1; %***** Key histogram enhancement step. Increases relative importance of high intensity pixels
        clrmapred(j+1,1) = 1/length(a)*j; %creates redscale colormap
    end

    BGIbins = sum(BGIcalc); %sums counts in all bins after enhancement

	%Calculates BGI
    if Foulcontrast == 'none'
        OBGIred = BGIbins/Totalbins*100;
    else
        OBGIred = 100-BGIbins/Totalbins*100;
    end

     clear a b uniqLevels j BGIbins BGIcalc Totalbins Foulcontrast

    %creates title with BGI for final image
    formatSpec = '%10.2f';
    redBGIstring = strcat('BGIred = ',num2str(OBGIred,formatSpec),'%');

    % repeats above process for greenplate BGI
N = 10;
[thresh, metric] = multithresh(greenplate, N);
 biomax = double(max(max(greenplate)));
    biomin = double(min(min(greenplate)));
    biorange = double(biomax - biomin);
    bioave = mean(mean(greenplate));
    controlmax = double(max(max(XBgreen)));
    controlmin = double(min(min(XBgreen)));
    controlrange = double(controlmax - controlmin);
    controlave = mean(mean(XBgreen));

    Foulcontrast = 'some';    
    if abs(bioave-controlave) < controlrange/4
        disp('fouled sample is very close to control sample - green')
        Foulcontrast = 'none';
        greenplate = greenplate - biomin;
        threshint = 255/10;
        thresh = [threshint 2*threshint 3*threshint 4*threshint 5*threshint 6*threshint 7*threshint 8*threshint 9*threshint];
    end
    
    Multithreshimagegreen = imquantize(greenplate,thresh);

    GrayscaleImagegreen = abs(Multithreshimagegreen/(N+1)-1);

    uniqLevels = unique(GrayscaleImagegreen(:));
[a,b]=hist(GrayscaleImagegreen(:),uniqLevels);
    a = fliplr(a);
    Totalbins = sum(a);

    clrmapgreen = zeros(10,3);
    for j = 1:length(a)
        BGIcalc(1,j) = a(1,j)*(j-1)*.1;
        clrmapgreen(j+1,2) = 1/length(a)*j;
    end
    BGIbins = sum(BGIcalc);

    if Foulcontrast == 'none'
        OBGIgreen = BGIbins/Totalbins*100;
    else
        OBGIgreen = 100 - BGIbins/Totalbins*100;
    end
    clear a b uniqLevels j BGIbins BGIcalc Totalbins Foulcontrast
    formatSpec = '%10.2f';
    greenBGIstring = strcat('BGIgreen = ',num2str(OBGIgreen,formatSpec),'%');

    
    % repeats above process for blueplate BGI
N = 10;
[thresh, metric] = multithresh(blueplate, N);
 biomax = double(max(max(blueplate)));
    biomin = double(min(min(blueplate)));
    biorange = double(biomax - biomin);
    bioave = mean(mean(blueplate));
    controlmax = double(max(max(XBblue)));
    controlmin = double(min(min(XBblue)));
    controlrange = double(controlmax - controlmin);
    controlave = mean(mean(XBblue));
    
    Foulcontrast = 'some';    
    if abs(bioave-controlave) < controlrange/4
        disp('fouled sample is very close to control sample - blue')
        Foulcontrast = 'none';
        blueplate = blueplate - biomin;
        threshint = 255/10;
        thresh = [threshint 2*threshint 3*threshint 4*threshint 5*threshint 6*threshint 7*threshint 8*threshint 9*threshint];
    end
    
    Multithreshimageblue = imquantize(blueplate,thresh);

    GrayscaleImageblue = abs(Multithreshimageblue/(N+1)-1);

    uniqLevels = unique(GrayscaleImageblue(:));
[a,b]=hist(GrayscaleImageblue(:),uniqLevels);

    a = fliplr(a);

    Totalbins = sum(a);

clrmapblue = zeros(10,3);
    for j = 1:length(a)
        BGIcalc(1,j) = a(1,j)*(j-1)*.1;
        clrmapblue(j+1,3) = 1/length(a)*j;
    end
    BGIbins = sum(BGIcalc);

    if Foulcontrast == 'none'
        OBGIblue = BGIbins/Totalbins*100;
    else
        OBGIblue = 100 - BGIbins/Totalbins*100;
    end
    clear a b uniqLevels j BGIbins BGIcalc Totalbins
    formatSpec = '%10.2f';
    blueBGIstring = strcat('BGIblue = ',num2str(OBGIblue,formatSpec),'%');
   
    
    % ---- End of RGB analysis

    
    % repeats above process for grayscale BGI
N = 10;
[thresh, metric] = multithresh(X, N);

    biomax = double(max(max(X)));
    biomin = double(min(min(X)));
    biorange = double(biomax - biomin);
    bioave = mean(mean(X));
    controlmax = double(max(max(XBred)));
    controlmin = double(min(min(XBred)));
    controlrange = double(controlmax - controlmin);
    controlave = mean(mean(XBred));
    
        Foulcontrast = 'some';
    if abs(bioave-controlave) < controlrange/4
        disp('fouled sample is very close to control sample - gray')
        Foulcontrast = 'none';
        X = X - biomin;
        threshint = 255/10;
        thresh = [threshint 2*threshint 3*threshint 4*threshint 5*threshint 6*threshint 7*threshint 8*threshint 9*threshint];
    end
  
    Multithreshimage = imquantize(X,thresh);

    GrayscaleImage = abs(Multithreshimage/(N+1)-1);

    
    uniqLevels = unique(GrayscaleImage(:));
    [a,b]=hist(GrayscaleImage(:),uniqLevels);
 
    a = fliplr(a);
    
    Totalbins = sum(a);
    clrmapgray = zeros(10,3);
    for j = 1:length(a)
        BGIcalc(1,j) = a(1,j)*(j-1)*.1;
        clrmapgray(j+1,1) = 1/length(a)*j;
        clrmapgray(j+1,2) = 1/length(a)*j;
        clrmapgray(j+1,3) = 1/length(a)*j;
    end
    BGIbins = sum(BGIcalc);

    if Foulcontrast == 'none'
        OBGI = BGIbins/Totalbins*100;
    else
        OBGI = 100 - BGIbins/Totalbins*100;
    end
    clear a b uniqLevels j BGIbins BGIcalc Totalbins
    BGIstring = strcat('BGI gray = ',num2str(OBGI,formatSpec),'%');

    % ----- End of BGI analysis


    
    clearvars a
     % Creates a BGI figure with enhanced images from each color channel and gray
    figure103 = figure(103);
    set(figure103,'units','normalized','position',[.05 .5 .45 .4])
    imagesc(GrayscaleImage), title(BGIstring)
    axis off
    colormap(clrmapgray)
        saveas(gcf,strcat(imagename,'BGIgray'),'tif'); %saves figure to a file
    figure104 = figure(104);
    set(figure104,'units','normalized','position',[.5 .5 .45 .4])
    imagesc(GrayscaleImagered), title(redBGIstring)
    axis off
    colormap(clrmapred)
        saveas(gcf,strcat(imagename,'BGIred'),'tif'); %saves figure to a file
    figure105 = figure(105);
    set(figure105,'units','normalized','position',[.05 .05 .45 .4])
    imagesc(GrayscaleImagegreen), title(greenBGIstring)
    axis off
    colormap(clrmapgreen)
        saveas(gcf,strcat(imagename,'BGIgreen'),'tif'); %saves figure to a file
    figure106 = figure(106);
    set(figure106,'units','normalized','position',[.5 .05 .45 .4])
    imagesc(GrayscaleImageblue), title(blueBGIstring)
    axis off
    colormap(clrmapblue)
        saveas(gcf,strcat(imagename,'BGIblue'),'tif'); %saves figure to a file
        
    
    
      qq = questdlg('Would you like to accept the data?','Coverage Map.','Yes','No','Yes');
    if isempty(qq) || strcmpi(qq,'no')
        disp('Error: Data not accepted')
        msgbox('Error: Data not accepted')
        close all
    else
      SampleNum(1,1) = SampleNum(1,1) + 1;

	% Records data to an Excel file
      sheet = 1;
      A = {SampleNum,fName,fName2,OBGI,OBGIred,OBGIgreen,OBGIblue,corners1(1,1),corners1(1,2),corners1(2,1),corners1(2,2),blank(1,1),blank(1,2),blank(2,1),blank(2,2)};
      xlRange = strcat('A',num2str(SampleNum+1));
      xlswrite(DataFile,A,sheet,xlRange);  
      close all

    end
end



