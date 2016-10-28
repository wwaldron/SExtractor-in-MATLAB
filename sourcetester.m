%% Test for Uncommon source function.
% This is designed to be used in conjunction with the finduncommonsources
% function.


%% Clear
close all
clear
clc


%% Declare Variables
% Notes on FITS Files here
% ib9ga1010_drc_slice.fits = WFC3,  275 nm
% iblwa1010_drz_slice.fits = WFC3, 1600 nm
% jb9g01010_drc_slice.fits =  ACS,  814 nm
% jb9g01020_drc_slice.fits =  ACS,  475 nm

% % The FITS Files of interest (full path recommended)
% fitsFiles = {'/home/wwaldron/Pictures/tmp/sliced_images/ib9ga1010_drc_slice.fits',...
% '/home/wwaldron/Pictures/tmp/sliced_images/iblwa1010_drz_slice.fits',...
% '/home/wwaldron/Pictures/tmp/sliced_images/jb9g01010_drc_slice.fits',...
% '/home/wwaldron/Pictures/tmp/sliced_images/jb9g01020_drc_slice.fits'};
% 
% % Helper files
% confFile  = '/home/wwaldron/Pictures/tmp/SExtractorFiles/eso_137-001.conf';
% paramFile = '/home/wwaldron/Pictures/tmp/SExtractorFiles/eso_137-001.param';
% convFile  = '/usr/share/sextractor/default.conv';
% nnFile    = '/usr/share/sextractor/default.nnw';

% The FITS Files of interest (full path recommended)
fitsFiles = {'/home/wwaldron/Documents/DoctoralResearch-Dissertation/Images/09-01-16_MAST_DATA/sliced_images/ib9ga1010_drc_slice.fits',...
'/home/wwaldron/Documents/DoctoralResearch-Dissertation/Images/09-01-16_MAST_DATA/sliced_images/iblwa1010_drz_slice.fits',...
'/home/wwaldron/Documents/DoctoralResearch-Dissertation/Images/09-01-16_MAST_DATA/sliced_images/jb9g01010_drc_slice.fits',...
'/home/wwaldron/Documents/DoctoralResearch-Dissertation/Images/09-01-16_MAST_DATA/sliced_images/jb9g01020_drc_slice.fits'};

% Helper files
confFile  = '/home/wwaldron/Documents/DoctoralResearch-Dissertation/Images/09-01-16_MAST_DATA/ESO_137-001_SExtractor_Files/mingSExtractor2.conf';
paramFile = '/home/wwaldron/Documents/DoctoralResearch-Dissertation/Images/09-01-16_MAST_DATA/ESO_137-001_SExtractor_Files/mingSExtractor.param';
convFile  = '/home/wwaldron/Software/SExtractor/ConvolutionKernals/default.conv';
nnFile    = '/home/wwaldron/Software/SExtractor/ConvolutionKernals/default.nnw';

% Name-Value Pair Arguments
nvPairs.SExtractorCommand = 'sextractor';
nvPairs.AnalysisThreshold = 10;
nvPairs.DetectThreshold   = [6,6,6,6];
nvPairs.MaxAngSepTrueSrc  = 0.2;
nvPairs.Gain              = [1.54,2.5,2.0,2.0];
nvPairs.SeeingFWHM        = [0.4 ,0.1,0.5,0.5];
nvPairs.DetectMinArea     = 4;
nvPairs.PixelScale        = 0;
nvPairs.RegOutputFormat   = 'SAOImageCirc';

data = finduncommonsources(fitsFiles,confFile,paramFile,convFile,nnFile,...
    nvPairs);

clearvars -except data