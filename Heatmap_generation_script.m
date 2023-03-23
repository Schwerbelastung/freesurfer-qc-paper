% Heatmap generation script for the heatmaps
% Authors: Baran Aydogan (coordinate conversion), Vesa Vahermaa (rest of
% the script)

% Make sure that Freesurfer executables are in the path
% Make sure that Nifti toolbox is in the path (It is here in case it is not 
% already installed: https://se.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
addpath('/.../MatLab_Libraries/NIfTI_20140122');
addpath('/.../MatLab_Libraries/create_mosaic_1.0.m'); % Potentially unnecessary
addpath('/Applications/freesurfer/7.1.1/matlab');
addpath('/Applications/freesurfer/7.1.1/fsfast/toolbox');
setenv( 'PATH',[getenv('PATH') ':/Applications/freesurfer/7.1.1/bin:/Applications/freesurfer/7.1.1/fsfast/bin:/Applications/freesurfer/7.1.1/mni/bin:/Users/vesavahermaa/opt/anaconda3/bin:/Users/vesavahermaa/opt/anaconda3/condabin:/Applications/freesurfer/7.1.1/bin:/Applications/freesurfer/7.1.1/fsfast/bin:/Applications/freesurfer/7.1.1/mni/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/opt/X11/bin'])
setenv FREESURFER_HOME /Applications/freesurfer/7.1.1

qc_folder       = '/Users/...'; % QC working folder
main_fs_folder  = '/Users/...';
out_folder      = '/Users/...'; % THIS SHOULD BE CREATED

%%


% Import CSV files for fixes

fixsubjects = readtable('/.../freesurfer-corr-export-names.csv','PreserveVariableNames',false);
fixvalues = readmatrix([qc_folder,'/freesurfer-corr-export-values.csv']);
number_of_errors =  int32(height(fixvalues));



%% Convert coordinates


for i=1:size(fixvalues,1)
    

    subject         = char(fixsubjects{i,1});
    RAS             = [fixvalues(i,1), fixvalues(i,2), fixvalues(i,3)]; % Example


    subjects_fs_folder  = [main_fs_folder,'/',subject];
    origMgz             = [subjects_fs_folder, '/mri/orig.mgz'];
    talaM3z             = [subjects_fs_folder, '/mri/transforms/talairach.m3z'];
    out                 = [out_folder, '/subject'];

    for n=1:size(RAS,1)

        % Convert RAS -> Freesurfer space
        norig = [out_folder, '/subject_fNorig_', num2str(n)];
        %mkdir(out_folder, ['subject_fNorig_', num2str(n)]);
        eval(['!mri_info --vox2ras ' origMgz '>' norig]);
        Norig = dlmread(norig);
        p = round(inv(Norig)*[RAS(n,1) RAS(n,2) RAS(n,3) 1]')+1; %#ok<MINV>
        p = p(1:3); % This is in real space. Need to multiply with srow (x,y,z) or its inverse, from test.hdr.hist (written on 24.11. during meeting)

        % Mark the point as 1 in Freesurfer space
        tmp                     = [out_folder, '/subject_tmp_', num2str(n) '.nii.gz'];
        eval(['!mri_convert ' origMgz ' ' tmp]);
        nii                     = load_untouch_nii(tmp);
        nii.img                 = nii.img*0;
        nii.img(p(1),p(2),p(3)) = 1;
        nii.hdr.dime.datatype   = 16;
        save_untouch_nii(nii,tmp);

        % Mark the point as 1 in MNI305 space
        eval(['!mri_convert ' tmp ' --apply_transform ' talaM3z ' -rt cubic --out_data_type float -oc 0 0 0 ' out '_' num2str(n) '.nii.gz']);
        nii               = load_untouch_nii([out '_' num2str(n) '.nii.gz']);
        [max_val,max_ind] = max(nii.img(:));
        nii.img     = nii.img*0;
        nii.img(max_ind)  = 1; 
        save_untouch_nii(nii,[out '_' num2str(i) '.nii.gz']);

        delete(tmp);
        delete(norig);
    end
end

%%

%########################
%### Sum error images  ##
%########################


imagesToSum = dir(out_folder + "/*.nii.gz");
imageFileNames = fullfile(out_folder, {imagesToSum.name});

% Get initial image

combinedErrorImage = load_untouch_nii(char(imageFileNames(1))); 

for i=2:size(imagesToSum)
    combinedErrorImage.img = combinedErrorImage.img + load_untouch_nii(char(imageFileNames(i))).img;
end

outputFilename_summedErrors = '/.../MNI305_ErrorHeatMap/GeneratedImages/CombinedErrors/mni.nii';

save_untouch_nii(combinedErrorImage,outputFilename_summedErrors);

%%

%###########################
%### Gaussian smoothing  ###
%###########################

imageToSmooth = load_untouch_nii(outputFilename_summedErrors);

imageToSmooth.img = imgaussfilt3(imageToSmooth.img,3); % Please modify the SIGMA parameter accordingly. SIGMA here is in the unit of voxel. If it is 2, then it is "2 voxels".

outputFilename_smoothedErrors = '/.../MNI305_ErrorHeatMap/GeneratedImages/SmoothedErrors/mni.nii';

save_untouch_nii(imageToSmooth,outputFilename_smoothedErrors);


%%

%######################################
%### Create combined overlay image  ###
%######################################

% Import MNI305 brain

mni305 = load_untouch_nii('/.../MNI305_ErrorHeatMap/GeneratedImages/MNI305/mni.nii');

% Import smoothed heatmap

heatmap = load_untouch_nii(outputFilename_smoothedErrors);

% Combine images

%% This takes a long time: skip if it's not needed

% %#########################################################
% %### Convert APARC-images to original anatomical space ###
% %#########################################################
% 
% % Set up Freesurfer dir
% 
% setenv SUBJECTS_DIR /...
% 
% for i=1:height(fixsubjects)
%     user = char(fixsubjects{i,1});
%     aseg_input_path = qc_folder + "/" + user + "/mri/aparc+aseg.mgz";
%     aseg_output_path = qc_folder + "/" + user + "/mri/aparc-in-rawavg.mgz";
%     rawavg_path = qc_folder + "/" + user + "/mri/rawavg.mgz";
%     disp("Processing APARC space conversion of image " + i + " of " + height(fixsubjects) + ":" + user);
%     eval(['!export SUBJECTS_DIR=/...']);
%     eval(["!mri_label2vol --seg " + aseg_input_path + " --temp " + rawavg_path + " --o " + aseg_output_path + " --regheader " + aseg_input_path]);
%     disp("Processing of subject " + user + " done.");
% end
% 
% %##############################################
% %### Convert APARC-images MGZ to NII images ###
% %##############################################
% 
% for i=1:height(fixsubjects)
%     user = char(fixsubjects{i,1});
%     aseg_input_path = qc_folder + "/" + user + "/mri/aparc-in-rawavg.mgz";
%     aseg_output_path = qc_folder + "/" + user + "/mri/aparc-in-rawavg.nii.gz";
%     disp("Processing APARC conversion from MGZ to NII for image " + i + " of " + height(fixsubjects) + ":" + user);
%     eval(['!mri_convert --out_orientation RAS ' + aseg_input_path + ' ' + aseg_output_path]);
%     disp("Processing of subject " + user + " done.");
% end


%%

% ##############################################
% ### Convert APARC-images MGZ to NII images ###
% ##############################################

for i=1:height(fixsubjects)
    user = char(fixsubjects{i,1});
    aseg_input_path = qc_folder + "/" + user + "/mri/aparc+aseg.mgz";
    aseg_output_path = qc_folder + "/" + user + "/mri/aparc+aseg.nii.gz";
    disp("Processing image " + i + " of " + height(fixsubjects) + ":" + user);
    eval(['!mri_convert ' + aseg_input_path + ' ' + aseg_output_path]);
    disp("Processing of subject " + user + " done.");
end




%%

%#########################################################################
%### Find out what brain areas are affected  (not used in final paper) ###
%#########################################################################

% Initialize error area array

errorAreas = zeros(1,number_of_errors);


for i=1:number_of_errors
    
    disp("Processing recorded error " + i + " of " + number_of_errors + "...");
    
    % Select user and related error (round error to get full integer)
    
    user = char(fixsubjects{i,1}); % Get current user
    error = fixvalues(i,:); % Note that error is not yet an integer
    
    % Fix error coordinates (must be positive) (TODO: This 
    %error = error + [128,128,128];
    %error = [error(1),error(2),256-error(3)];

    % Determine aparc file path and target unzip path
    
    aparc_aseg_path = qc_folder + "/" + user + "/mri/aparc+aseg.nii.gz";
    unzip_path = '/.../Unzip_Folder_For_Matlab';
    
    % Unzip and load the NII file

    gunzip(aparc_aseg_path,unzip_path)
    aparc_path = append(unzip_path,'/aparc+aseg.nii');
    aparc_image = load_untouch_nii(aparc_path).img;
    
    % Read origin file to figure out offset
    
    orig_path_mgz = qc_folder + "/" + user + "/mri/orig.mgz";
    aparc_image_mgz = MRIread(char(orig_path_mgz));
    
    % Read offset data from orig image
    
    orig_offset_r = aparc_image_mgz.c_r;
    orig_offset_a = aparc_image_mgz.c_a;
    orig_offset_s = aparc_image_mgz.c_s;
    
    % Translate error with offset info
    
    error2 = round([128 - error(1) + orig_offset_r, 128 - error(3) + orig_offset_s, 128 + error(2) - orig_offset_a]);
            
    % Get label: Apparently X and Y coordinates are reversed in aparc
    
    errorLabel = aparc_image(error2(1),error2(2),error2(3));
    
    % Check if error label is 0
    
    flag = 0; % Initialize flag for breaking out of nested loop
    
    while errorLabel == 0
        disp("Image " + i + " came up with label 0. Finding non-zero value..");
        
        side_length = 7; % Determine the size of the lookup cube, e.g., 3x3x3 voxels. Must be an odd number!
        start_index = 0-((side_length-1)/2)-1; % Determine which voxel is checked first
        
        for x=1:side_length
           for y=1:side_length
               for z=1:side_length
                   errorLabel = aparc_image(error2(1)+start_index+x,error2(2)+start_index+y,error2(3)+start_index+z);
                   %disp([error2(1)+start_index+x,error2(2)+start_index+y,error2(3)+start_index+z]);
                   
                   if errorLabel ~= 0
                       
                       disp("Non-zero error label found for error " + i + ", user " + user +": " + errorLabel + "!");
                       disp("Original error coordinates were: " + error(1) + "," + error(2) + "," + error(3));
                       flag=1;
                       
                       break
                       
                   end
                   
                   %
                   if flag == 1
                       break
                   end
                   %
                   
               end
               
               %
               if flag == 1
                   break
               end
               %
               
           end
           
           %
           if flag == 1
               break
           end
           %
           
        end
          
        if errorLabel == 0
            
            disp("Non-zero error label could not be found for error " + i + ", user " + user +": leaving as " + errorLabel + ".");
            disp("Original error coordinates were: " + error(1) + "," + error(2) + "," + error(3));
            
        end
        
        
        
        break
    end
        
    
    % Add label to array
    
    errorAreas(i) = errorLabel;
    
    % Delete file
    
    delete(aparc_path);
    
%     figure
%     imagesc(aparc_image)
%     imagesc(flip(squeeze(aparc_image(:,140,:))'))
    
    disp("Processing of recorded error " + i + " done.");
    disp("----------");
    
end


edges = unique(errorAreas).';
counts = histc(errorAreas(:), edges);

errorAreasTable = [edges counts];

errorAreasTable;



%aseg_input_path = qc_folder + "/" + user + "/mri/aseg.mgz";
%aseg_output_path = qc_folder + "/" + user + "/mri/aseg.nii.gz";

% Convert aseg

% eval(['!mri_convert --out_orientation RAS ' + aseg_input_path + ' ' + aseg_output_path])
%%

% %###############################################
% %### Create histogram of error distribution  ###
% %###############################################

% Note: This section is dependent on the files loaded at the start of the
% script

% Note 2: THIS PART REQUIRES THE MANUAL ADDITION OF ZERO ERROR SUBJECTS.
% This can't automatically be determined from the data loaded, since the
% subjects with 0 errors aren't included in the error location count
% script. Check the value from the error recording Excel sheet!

% Last changed: 26.9.2022
zeroErrorSubjects = 19;

% Load required data and convert it to proper format
subjectsArray = table2array(fixsubjects);

% Calculate the number of errors for each subject
[vectorErrorCount, vectorSubjectID] = groupcounts(subjectsArray);

% Create table out of the subjects with errors and 
errorCountTable = table(vectorSubjectID, vectorErrorCount);

% Create vector counting error occurrences
errorOccurrences = table2array(errorCountTable(:,2));

% Add zero error count to the array
% NOTE: MAKE SURE THE COUNT IS RIGHT!

for i=1:zeroErrorSubjects
    errorOccurrences(end+1) = 0;
end

% Determine max number of errors in an image

maxErrors = max(unique(errorOccurrences));

% Initialize and populate barplot table

errorBarplotTable = zeros(maxErrors+1,2); % +1 to account for zero
errorBarplotTable(1,2) = zeroErrorSubjects;

for i=1:maxErrors
    errorBarplotTable(i+1,1) = i;
    errorBarplotTable(i+1,2) = sum(errorOccurrences(:) == i);
end

errorBarDataX = 0:maxErrors;
errorBarDataY = zeros(1,maxErrors+1);
    
for i=1:maxErrors+1
    errorBarDataY(i) = errorBarplotTable(i,2);
end

% Create bar plot

b = bar(errorBarDataX,errorBarDataY, 'FaceColor', [0.27,0.48,0.76]);

% Leave room on top

ylim([0 35]);

% Add numbers on top of bars

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')

% Add titles

title('Distribution of errors across all inspected MRI images (n=167)');
xlabel('Number of identified errors in an inspected image');
ylabel('Frequency of images with X identified errors');



