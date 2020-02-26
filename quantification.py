# This script:
# 1.) Defines cortical areas on each section based on segmentation masks and plots it.
# 2.) Defines cortical depth on each section based on segmentation masks and plots it.
# 3.) Defines radial distance on each section based on cortical depth and plots it.
# 4.) Quantifies the expression of adrenergic receptors in astrocytes as a function of cortical depth
#     at three radial positions and plots it.

import os
os.chdir('/home/jovyan/MH_DDD/')
from fnmatch import fnmatch
import pickle
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt; plt.style.use('ggplot')
import yaml
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.geometry import LineString
from matplotlib.pyplot import figure
from nested_lookup import nested_lookup
import seaborn as sns

# Load the data and segmentation masks:

# Get all evaluation files:
root = '../data/OB_ADR/'
pattern = "Objects_Population - All Nuclei.txt"
allFiles = []
measurementNames = []
slideNames = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            allFiles.append(os.path.join(path, name))
            measurementNames.append(str.split(allFiles[-1], '/')[3])
            slideNames.append(str.split(measurementNames[-1], '__')[0])
            
pattern = "Objects_Population - Astrocytes.txt"
allFiles_Astrocytes = []
measurementNames_Astrocytes = []
slideNames_Astrocytes = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            allFiles_Astrocytes.append(os.path.join(path, name))
            measurementNames_Astrocytes.append(str.split(allFiles_Astrocytes[-1], '/')[3])
            slideNames_Astrocytes.append(str.split(measurementNames[-1], '__')[0])

sectionNumbers = np.array((2,1,3,1,2,2,1,3))
numberOfSections = len(sectionNumbers)
# Plot the quanitification:
f, axis = plt.subplots(numberOfSections,5, figsize=(35,numberOfSections*7))
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['axes.facecolor'] = 'white'
dotSize = 0.5
count = 0
for slide in range(len(allFiles)):
    
    print(slide)
    sectionNumber = sectionNumbers[slide]
    
    # Import nuclei data:
    nuclei_all = pd.read_csv(allFiles[slide], sep = '\t' , skiprows = 8, header = 1, error_bad_lines=False)
    nuclei_data = np.asarray(nuclei_all[['Position X [µm]', 'Position Y [µm]']])

   # Astrocyte data:
    astrocytes_all = pd.read_csv(allFiles_Astrocytes[slide], sep = '\t' , skiprows = 8, header = 1)
    astrocytes_data = np.asarray(astrocytes_all[['Position X [µm]', 'Position Y [µm]']])
            
    with open("/home/jovyan/data/OB_ADR/SegmentationMasks/" + slideNames[slide] + '_Segmentation.yml', 'r') as stream:
            segmentationMasks = yaml.safe_load(stream)          
    offset = pd.read_csv("/home/jovyan/data/OB_ADR/SegmentationMasks/" + slideNames[slide] + '_offset.txt', header = None)
    offset = np.array(offset)[0]
    # Correct offset for one slide:
    if slide == 0:
        offset[1] = offset[1] - 0.5
    
    cortex = np.zeros(np.shape(nuclei_all)[0])
    positions =  nested_lookup('position', segmentationMasks)
    wm_masks = nested_lookup('wm', segmentationMasks)
    pia_masks = nested_lookup('pia', segmentationMasks)
    # Left Hemisphere:
    firstSplit = str.split(positions[(2*sectionNumber)-2][0], ';')
    polygon = Polygon(np.asarray([str.split(firstSplit[i], ' ') for i in range(len(firstSplit))]).astype(np.float))
    for j in range(len(cortex)):
        point = Point(nuclei_data[j,0]/1000+offset[0], nuclei_data[j,1]/1000+offset[1])
        if polygon.contains(point):
            cortex[j] = 1  
    # Right Hemisphere:
    firstSplit = str.split(positions[(2*sectionNumber)-1][0], ';')
    polygon = Polygon(np.asarray([str.split(firstSplit[i], ' ') for i in range(len(firstSplit))]).astype(np.float))
    for j in range(len(cortex)):
        point = Point(nuclei_data[j,0]/1000+offset[0], nuclei_data[j,1]/1000+offset[1])
        if polygon.contains(point):
            cortex[j] = 2  
            
    # Do the same for astrocytes:
    cortex_A = np.zeros(np.shape(astrocytes_all)[0])
    # Left Hemisphere:
    firstSplit = str.split(positions[(2*sectionNumber)-2][0], ';')
    polygon = Polygon(np.asarray([str.split(firstSplit[i], ' ') for i in range(len(firstSplit))]).astype(np.float))
    for j in range(len(cortex_A)):
        point = Point(astrocytes_data[j,0]/1000+offset[0], astrocytes_data[j,1]/1000+offset[1])
        if polygon.contains(point):
            cortex_A[j] = 1  
    # Right Hemisphere:
    firstSplit = str.split(positions[(2*sectionNumber)-1][0], ';')
    polygon = Polygon(np.asarray([str.split(firstSplit[i], ' ') for i in range(len(firstSplit))]).astype(np.float))
    for j in range(len(cortex_A)):
        point = Point(astrocytes_data[j,0]/1000+offset[0], astrocytes_data[j,1]/1000+offset[1])
        if polygon.contains(point):
            cortex_A[j] = 2  
    
    colours = np.repeat('black', np.shape(nuclei_data)[0])
    colours[cortex == 1] = 'blue'
    colours[cortex == 2] = 'red' 
    
    # Also mark white matter surface and bubbles in figures:
    wml = str.split(wm_masks[(2*sectionNumber)-2][0], ';')
    wml = LineString(np.asarray([str.split(wml[i], ' ') for i in range(len(wml))]).astype(np.float))
    x_wml, y_wml = wml.xy
    x_wml = [(x_wml[i]-offset[0])*1000 for i in range(len(x_wml))]
    y_wml = [(y_wml[i]-offset[1])*1000 for i in range(len(y_wml))]
    pial = str.split(pia_masks[(2*sectionNumber)-2][0], ';')
    pial = LineString(np.asarray([str.split(pial[i], ' ') for i in range(len(pial))]).astype(np.float)) 
    x_pial, y_pial = pial.xy
    x_pial = [(x_pial[i]-offset[0])*1000 for i in range(len(x_pial))]
    y_pial = [(y_pial[i]-offset[1])*1000 for i in range(len(y_pial))]
    wmr = str.split(wm_masks[(2*sectionNumber)-1][0], ';')
    wmr = LineString(np.asarray([str.split(wmr[i], ' ') for i in range(len(wmr))]).astype(np.float))
    x_wmr, y_wmr = wmr.xy
    x_wmr = [(x_wmr[i]-offset[0])*1000 for i in range(len(x_wmr))]
    y_wmr = [(y_wmr[i]-offset[1])*1000 for i in range(len(y_wmr))]
    piar = str.split(pia_masks[(2*sectionNumber)-1][0], ';')
    piar = LineString(np.asarray([str.split(piar[i], ' ') for i in range(len(piar))]).astype(np.float)) 
    x_piar, y_piar = piar.xy
    x_piar = [(x_piar[i]-offset[0])*1000 for i in range(len(x_piar))]
    y_piar = [(y_piar[i]-offset[1])*1000 for i in range(len(y_piar))]
    
    # Assign one nuclei to right cortex if none are there to fix bug:
    if sum(cortex == 1) == 0:
        cortex[1] = 1
        cortex[2] = 1
        
    if sum(cortex_A == 1) == 0:
        cortex_A[1] = 1
        cortex_A[2] = 1
    
    # Compute new coordinates for left hemisphere:

    radialDistanceL = np.repeat(None, len(nuclei_data[:,0]))
    for j in np.nditer(np.where(cortex == 1)):
        point = Point(nuclei_data[:,0][j]/1000+offset[0], nuclei_data[:,1][j]/1000 + offset[1])
        radialDistanceL[j] = 1-wml.project(point, normalized = True)
    wmDistanceL = np.repeat(None, len(nuclei_data[:,0]))
    for j in np.nditer(np.where(cortex == 1)):
        point = Point(nuclei_data[:,0][j]/1000+offset[0], nuclei_data[:,1][j]/1000 + offset[1])
        wmDistanceL[j] = wml.distance(point)
    pialDistanceL = np.repeat(None, len(nuclei_data[:,0]))
    for j in np.nditer(np.where(cortex == 1)):
        point = Point(nuclei_data[:,0][j]/1000+offset[0], nuclei_data[:,1][j]/1000 + offset[1])
        pialDistanceL[j] = pial.distance(point)
    norm_corticalDepthL = np.repeat(None, len(nuclei_data[:,0]))
    for j in np.nditer(np.where(cortex == 1)):
        norm_corticalDepthL[j] = wmDistanceL[j]/(wmDistanceL[j] + pialDistanceL[j])

    # Compute new coordinates for right hemisphere:

    radialDistanceR = np.repeat(None, len(nuclei_data[:,0]))
    for j in np.nditer(np.where(cortex == 2)):
        point = Point(nuclei_data[:,0][j]/1000+offset[0], nuclei_data[:,1][j]/1000 + offset[1])
        radialDistanceR[j] = wmr.project(point, normalized = True)
    wmDistanceR = np.repeat(None, len(nuclei_data[:,0]))
    for j in np.nditer(np.where(cortex == 2)):
        point = Point(nuclei_data[:,0][j]/1000+offset[0], nuclei_data[:,1][j]/1000 + offset[1])
        wmDistanceR[j] = wmr.distance(point)
    pialDistanceR = np.repeat(None, len(nuclei_data[:,0]))
    for j in np.nditer(np.where(cortex == 2)):
        point = Point(nuclei_data[:,0][j]/1000+offset[0], nuclei_data[:,1][j]/1000 + offset[1])
        pialDistanceR[j] = piar.distance(point)
    norm_corticalDepthR = np.repeat(None, len(nuclei_data[:,0]))
    for j in np.nditer(np.where(cortex == 2)):
        norm_corticalDepthR[j] = wmDistanceR[j]/(wmDistanceR[j] + pialDistanceR[j])

    # We want the radialDistanceR to increase with x coordinate:
    if np.corrcoef(nuclei_data[:,0][cortex == 2], radialDistanceR[cortex == 2].astype(float))[1,0] < 0:
        radialDistanceR[cortex == 2] = 1-radialDistanceR[cortex == 2]
    # We want the radialDistanceL to decrease with x coordinate:
    if np.corrcoef(nuclei_data[:,0][cortex == 1], radialDistanceL[cortex == 1].astype(float))[1,0] > 0:
        radialDistanceL[cortex == 1] = 1-radialDistanceL[cortex == 1]
        
    # Finally get number of spots in each astrocyte as a function on cortical depth in mid-cortex_A:
    
        # Compute new coordinates for left hemisphere, Astrocytes_only:

    radialDistance_A_L = np.repeat(None, len(astrocytes_data[:,0]))
    for j in np.nditer(np.where(cortex_A == 1)):
        point = Point(astrocytes_data[:,0][j]/1000+offset[0], astrocytes_data[:,1][j]/1000 + offset[1])
        radialDistance_A_L[j] = 1-wml.project(point, normalized = True)
    wmDistanceL = np.repeat(None, len(astrocytes_data[:,0]))
    for j in np.nditer(np.where(cortex_A == 1)):
        point = Point(astrocytes_data[:,0][j]/1000+offset[0], astrocytes_data[:,1][j]/1000 + offset[1])
        wmDistanceL[j] = wml.distance(point)
    pialDistanceL = np.repeat(None, len(astrocytes_data[:,0]))
    for j in np.nditer(np.where(cortex_A == 1)):
        point = Point(astrocytes_data[:,0][j]/1000+offset[0], astrocytes_data[:,1][j]/1000 + offset[1])
        pialDistanceL[j] = pial.distance(point)
    norm_corticalDepth_A_L = np.repeat(None, len(astrocytes_data[:,0]))
    for j in np.nditer(np.where(cortex_A == 1)):
        norm_corticalDepth_A_L[j] = wmDistanceL[j]/(wmDistanceL[j] + pialDistanceL[j])

    # Compute new coordinates for right hemisphere, Astrocytes_only:

    radialDistance_A_R = np.repeat(None, len(astrocytes_data[:,0]))
    for j in np.nditer(np.where(cortex_A == 2)):
        point = Point(astrocytes_data[:,0][j]/1000+offset[0], astrocytes_data[:,1][j]/1000 + offset[1])
        radialDistance_A_R[j] = wmr.project(point, normalized = True)
    wmDistanceR = np.repeat(None, len(astrocytes_data[:,0]))
    for j in np.nditer(np.where(cortex_A == 2)):
        point = Point(astrocytes_data[:,0][j]/1000+offset[0], astrocytes_data[:,1][j]/1000 + offset[1])
        wmDistanceR[j] = wmr.distance(point)
    pialDistanceR = np.repeat(None, len(astrocytes_data[:,0]))
    for j in np.nditer(np.where(cortex_A == 2)):
        point = Point(astrocytes_data[:,0][j]/1000+offset[0], astrocytes_data[:,1][j]/1000 + offset[1])
        pialDistanceR[j] = piar.distance(point)
    norm_corticalDepth_A_R = np.repeat(None, len(astrocytes_data[:,0]))
    for j in np.nditer(np.where(cortex_A == 2)):
        norm_corticalDepth_A_R[j] = wmDistanceR[j]/(wmDistanceR[j] + pialDistanceR[j])

    # We want the radialDistance_A_R to increase with x coordinatE:
    if np.corrcoef(astrocytes_data[:,0][cortex_A == 2], radialDistance_A_R[cortex_A == 2].astype(float))[1,0] < 0:
        radialDistance_A_R[cortex_A == 2] = 1-radialDistance_A_R[cortex_A == 2]
    # We want the radialDistance_A_L to decrease with x coordinate:
    if np.corrcoef(astrocytes_data[:,0][cortex_A == 1], radialDistance_A_L[cortex_A == 1].astype(float))[1,0] > 0:
        radialDistance_A_L[cortex_A == 1] = 1-radialDistance_A_L[cortex_A == 1]
    
    # Get mean expression of all astrocytes in ten cortical depth bins and with radial distance between 0.33,0.66:
    
    spots_520 = np.repeat(np.NaN,10)
    spots_570 = np.repeat(np.NaN,10)
    spots_650 = np.repeat(np.NaN,10)
    
    for x in range(10):
        subset_L = [radialDistance_A_L[cortex_A == 1][i] >= 0.33 and
                    radialDistance_A_L[cortex_A == 1][i] <= 0.66 and
                    norm_corticalDepth_A_L[cortex_A == 1][i] <= (x+1)*0.1 and
                    norm_corticalDepth_A_L[cortex_A == 1][i] >= x*0.1 for i in range(sum(cortex_A == 1))]
        
        subset_R = [radialDistance_A_R[cortex_A == 2][i] >= 0.33 and
                    radialDistance_A_R[cortex_A == 2][i] <= 0.66 and
                    norm_corticalDepth_A_R[cortex_A == 2][i] <= (x+1)*0.1 and
                    norm_corticalDepth_A_R[cortex_A == 2][i] >= x*0.1 for i in range(sum(cortex_A == 2))]
    
        spots_520[x] = np.mean((np.concatenate((np.array(astrocytes_all['Astrocytes - Number of Spots 520- per Cell '][cortex_A == 1][subset_L]), np.array(astrocytes_all['Astrocytes - Number of Spots 520- per Cell '][cortex_A == 2][subset_R])))))
        spots_570[x] = np.mean((np.concatenate((np.array(astrocytes_all['Astrocytes - Number of Spots 570- per Cell '][cortex_A == 1][subset_L]), np.array(astrocytes_all['Astrocytes - Number of Spots 570- per Cell '][cortex_A == 2][subset_R])))))
        spots_650[x] = np.mean((np.concatenate((np.array(astrocytes_all['Astrocytes - Number of Spots 650- per cell '][cortex_A == 1][subset_L]), np.array(astrocytes_all['Astrocytes - Number of Spots 650- per cell '][cortex_A == 2][subset_R])))))
        
    cmap = sns.cubehelix_palette(as_cmap=True)
    
    axis[count,0].scatter(nuclei_data[:,0], nuclei_data[:,1], c = colours, s = 0.05)
    axis[count,0].plot(x_pial, y_pial, color = 'blue')
    axis[count,0].plot(x_wml, y_wml, color = 'orange')  
    axis[count,0].plot(x_piar, y_piar, color = 'blue')
    axis[count,0].plot(x_wmr, y_wmr, color = 'orange')    
    axis[count,0].set_title('Left (red) and Right (green) Cortex Segmentation \n' + slideNames[slide] + ' Section' + str(sectionNumber))
    axis[count,1].scatter(astrocytes_data[:,0], astrocytes_data[:,1], c = 'orange', s = 0.5)  
    axis[count,1].set_title('Astrocyte Positions \n' + slideNames[slide] + ' Section' + str(sectionNumber))
    axis[count,2].scatter(nuclei_data[:,0][cortex == 1], nuclei_data[:,1][cortex == 1], c=norm_corticalDepthL[cortex == 1], s=0.05, cmap=cmap)
    axis[count,2].scatter(nuclei_data[:,0][cortex == 2], nuclei_data[:,1][cortex == 2], c=norm_corticalDepthR[cortex == 2], s=0.05, cmap=cmap)
    axis[count,2].set_title('Cortical Depth \n' + slideNames[slide] + ' Section' + str(sectionNumber))
    axis[count,3].scatter(nuclei_data[:,0][cortex == 1], nuclei_data[:,1][cortex == 1], c=radialDistanceL[cortex == 1], s=0.05, cmap=cmap)
    axis[count,3].scatter(nuclei_data[:,0][cortex == 2], nuclei_data[:,1][cortex == 2], c=radialDistanceR[cortex == 2], s=0.05, cmap=cmap)
    axis[count,3].set_title('Radial Distance \n' + slideNames[slide] + ' Section' + str(sectionNumber))
    axis[count,4].scatter(np.array((0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95)), spots_520, label = '520', c = 'red')
    axis[count,4].scatter(np.array((0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95)), spots_570, label = '570', c = 'blue')
    axis[count,4].scatter(np.array((0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95)), spots_650, label = '650', c = 'yellow')
    axis[count,4].set_xlabel('Cortical Depth')
    axis[count,4].set_ylabel('Spot Counts')
    axis[count,4].legend()
    axis[count,4].set_title('Spot counts as function on cortical depth in mid-cortex')
    count = count + 1
f.savefig('test.png')    

# Defines cortical areas on each section based on segmentation masks and plots it.

