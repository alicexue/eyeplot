import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import csv
import os
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

"""
Written by Alice Xue, August 2017
Modified in July 2018
Modified in November 2019
Modified in December 2020

Analyzes *.asc eyetracking output files converted from Eyelink *.edf files
Retrieves total fixation duration and total number of fixations on each interest area, which must be defined in the *.asc files 

Fixation durations are calculated as end time - start time
EyeLink DataViewer adds 1 or 2 milliseconds to that value for each fixation
Total fixations differ from the EyeLink DataViewer's calculations by a multiple of the number of fixations

(0, 0) for fixation position in EDF is top left corner.

TRIALID values in EDF file must be numbers (can start with 0 or 1)
Interest areas must be defined correctly
"""

class EyeData:
    def __init__(self, fixationData, iareaData, displayCoords):
        self.fixationData = fixationData # data frame
        self.iareaData = iareaData # data frame
        self.displayCoords = displayCoords # x1, y1, x2, y2


def get_trial_times(asc_name, trial_marker_start='TRIALID', trial_marker_end='TRIAL_RESULT'):
    """
    Parses asc_name line by line to get start and end times of each trial
    Returns 2 dictionaries
        trialStartTimes {dict}: keys are trial numbers (based on TRIALID in asc file), values are start times
        trialEndTimes {dict}: keys are trial numbers (based on TRIALID in asc file), values are end times
    """
    
    trialStartTimes = {}
    trialEndTimes = {}
    with open(asc_name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in reader:
            if len(row) > 0:
                r = np.asarray(row)
                if r[0] == 'MSG' and ('TRIALID' in r[1]):
                    currentTrial = int(r[1].split()[-1])
                if r[0] == 'MSG' and (trial_marker_start in r[1]) and currentTrial not in trialEndTimes.keys():
                    trialStartTimes[currentTrial] = int(r[1].split()[0])
                if r[0] == 'MSG' and (trial_marker_end in r[1]) and currentTrial not in trialEndTimes.keys():
                    trialEndTimes[currentTrial] = int(r[1].split()[0])
    return trialStartTimes, trialEndTimes



def get_iarea_coordinates(asc_name):
    """
    Parses asc_name line by line to get IAREAs for each trial
    Returns allTrialsIAREA {dataframe}: stores trial number, IAREA id, and the IAREA coordinates
    """
    
    allTrialsIAREAList = []
    with open(asc_name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in reader:
            if len(row) > 0:
                r = np.asarray(row)
                if r[0] == 'MSG' and ('TRIALID' in r[1]):
                    currentTrial = int(r[1].split()[-1])
                elif r[0] == 'MSG' and ("IAREA" in r[1]):
                    iareaMessage = r[1].split()
                    iareaType = iareaMessage[3]
                    iareaCoords = iareaMessage[5:9]
                    iareaName = iareaMessage[-1]
                    currIAREA = {}
                    currIAREA['id'] = iareaName
                    currIAREA['type'] = iareaType
                    currIAREA['LeftX'] = float(iareaCoords[0])
                    currIAREA['UpperY'] = float(iareaCoords[1])
                    currIAREA['RightX'] = float(iareaCoords[2])
                    currIAREA['LowerY'] = float(iareaCoords[3])
                    currIAREA['trialN'] = currentTrial
                    allTrialsIAREAList.append(currIAREA)
    allTrialsIAREA = pd.DataFrame(allTrialsIAREAList)
    allTrialsIAREA['trialN'] = allTrialsIAREA['trialN'].astype(int)
    return allTrialsIAREA


def process_edf(asc_name, trial_marker_start='TRIALID', trial_marker_end='TRIAL_RESULT'):
    """
    Parses asc_name line by line to get fixation data for each trial
    Saves each fixation to a dataframe and stores its avg x and y coords, duration, and whether it's on an interest area
    """
    
    csv_name = asc_name[:-4] + '_FixationDataSummary.csv'
    trialStartTimes, trialEndTimes = get_trial_times(asc_name, trial_marker_start, trial_marker_end)
    allTrialsIAREA = get_iarea_coordinates(asc_name)

    print("Processing asc file now...")
    allTrialsFixationDataList = []

    if os.path.exists(asc_name):
        with open(asc_name, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
            trialFixationData = []

            trialNums = trialStartTimes.keys()
            sorted(trialNums)
            trialNumbers = list(trialNums)
            fixationStartTrial = trialNumbers[0]
            currentTrial = trialNumbers[0]
            for row in reader:
                isTrial = False

                if len(row) > 0:
                    r = np.asarray(row)
                    if r[0] == 'MSG' and ("TRIALID" in r[1]):
                        currentTrial = int(r[1].split(' ')[-1])
                    elif r[0] == 'MSG' and 'DISPLAY_COORDS' in r[1]:
                        displaycoords = r[1].split()[3:]
                        for i in range(0, len(displaycoords)):
                            displaycoords[i] = int(displaycoords[i])
                    elif 'SFIX' in r[0]:
                        fixationStartTrial = currentTrial
                    elif 'EFIX' in r[0]:
                        fixationStarted = False
                        floatR = []
                        for info in r:
                            try:
                                floatR.append(float(info))
                            except ValueError:
                                floatR.append(info)
                        fixationStartTime = float(floatR[0].split()[2])
                        fixationEndTime = float(floatR[1])
                        # key is trial number, value is duration
                        duration = {}
                        fixationType = {}

                        fixationTrial = currentTrial

                        # fixation starts before trial does and ends during the trial
                        if fixationStartTime <= trialStartTimes[currentTrial] <= fixationEndTime <= trialEndTimes[
                            currentTrial]:
                            duration[fixationTrial] = fixationEndTime - trialStartTimes[currentTrial]
                            fixationType[fixationTrial] = 'started before trial'
                            # fixation starts during previous trial and ends during current trial - have not taken into account fixations that span multiple trials
                            if currentTrial != trialNumbers[0] and trialEndTimes[
                                currentTrial - 1] >= fixationStartTime >= trialStartTimes[currentTrial - 1]:
                                fixationTrial = currentTrial - 1
                                duration[fixationTrial] = trialEndTimes[currentTrial - 1] - fixationStartTime
                                fixationType[fixationTrial] = 'started in previous trial'
                        # fixation occurs during trial
                        elif fixationStartTime >= trialStartTimes[currentTrial] and fixationEndTime <= trialEndTimes[
                            currentTrial]:
                            duration[fixationTrial] = fixationEndTime - fixationStartTime
                            fixationType[fixationTrial] = 'started and ended during trial'
                        # fixation begins during the trial and ends after the trial
                        elif fixationEndTime >= trialEndTimes[currentTrial] >= fixationStartTime >= trialStartTimes[
                            currentTrial]:
                            duration[fixationTrial] = trialEndTimes[currentTrial] - fixationStartTime
                            fixationType[fixationTrial] = 'ended after trial'
                        if fixationStartTime <= trialStartTimes[currentTrial]:
                            foundFixationTrial = False
                            # fixation started before a trial and ended after the trial - not the current trial but an earlier trial
                            # need to account for fixation in the earlier trial AND in the current trial
                            for trialN in reversed(trialNumbers[:currentTrial]):
                                if fixationStartTime <= trialStartTimes[trialN] and fixationEndTime >= trialEndTimes[
                                    trialN]:
                                    duration[trialN] = trialEndTimes[trialN] - trialStartTimes[trialN]
                                    fixationType[trialN] = 'started in previous trial and ended after trial'
                                    foundFixationTrial = True

                        for fixTrial in duration:
                            if duration[fixTrial] > 0:
                                fixationsOnComponents = []
                                iarea_name_list = allTrialsIAREA.loc[(allTrialsIAREA['trialN'] == fixTrial)][
                                    'id'].values
                                for cName in iarea_name_list:
                                    currIAREA = allTrialsIAREA.loc[
                                        (allTrialsIAREA['trialN'] == fixTrial) & (allTrialsIAREA['id'] == cName)]
                                    i = currIAREA.index.values[0]
                                    leftX = currIAREA.at[(i, 'LeftX')]
                                    lowerY = currIAREA.at[(i, 'LowerY')]
                                    rightX = currIAREA.at[(i, 'RightX')]
                                    upperY = currIAREA.at[(i, 'UpperY')]
                                    iareaType = currIAREA.at[(i, 'type')]
                                    avgFixationX = floatR[3]
                                    avgFixationY = floatR[4]

                                    fixationOnComponent = check_fixation_in_bounds(iareaType, avgFixationX,
                                                                                   avgFixationY, leftX, rightX, lowerY,
                                                                                   upperY)

                                    """
                                    b1=avgFixationX >= leftX and avgFixationX <= rightX
                                    b2=(avgFixationY >= upperY and avgFixationY <= lowerY) or (avgFixationY <= upperY and avgFixationY >= lowerY)
                                    # upper and lower Y's might be flipped
                                    if b1 and b2:
                                        fixationOnComponent = True
                                    """
                                    fixationsOnComponents.append(fixationOnComponent)

                                fixationData = {}
                                fixationData['duration'] = duration[fixTrial]
                                fixationData['avg x pos'] = floatR[3]
                                fixationData['avg y pos'] = floatR[4]
                                fixationData['type'] = fixationType[fixTrial]
                                for i in range(0, len(iarea_name_list)):
                                    fixationData['fix on ' + iarea_name_list[i]] = fixationsOnComponents[i]
                                fixationData['trialN'] = fixTrial
                                allTrialsFixationDataList.append(fixationData)

        allTrialsFixationData = pd.DataFrame(allTrialsFixationDataList)
        print("Data processing completed.")
        return EyeData(allTrialsFixationData, allTrialsIAREA, displaycoords)
    else:
        print("File not found: %s" % (asc_name))


def get_screen_width_height(displayCoords):
    # displayCoords: [x1, y1, x2, y2]
    SCREENWIDTH = displayCoords[2] - displayCoords[0] + 1
    SCREENHEIGHT = displayCoords[3] - displayCoords[1] + 1
    return SCREENWIDTH, SCREENHEIGHT
        

def plot_trial_fixations(eyedata):
    
    """
    If using interactive python, can plot each fixation and the interest areas
    """
    
    allTrialsFixationData = eyedata.fixationData
    allTrialsIAREA = eyedata.iareaData
    displaycoords = eyedata.displayCoords
    SCREENWIDTH, SCREENHEIGHT = get_screen_width_height(eyedata.displayCoords)
    scale = 3
    
    print("Plotting fixations now...")
    trialNums = np.unique(allTrialsFixationData['trialN'].values).astype(int)
    for trial in trialNums:
        if trialNums[0] == 0:
            i = trial + 1
        else:
            i = trial
        currTrial = allTrialsFixationData.loc[allTrialsFixationData['trialN'] == trial]
        plt.figure(i)
        plt.plot([SCREENWIDTH / 2, SCREENWIDTH / 2], [0, SCREENHEIGHT], linewidth=0.5, color='k')
        plt.plot([0, SCREENWIDTH], [SCREENHEIGHT / 2, SCREENHEIGHT / 2], linewidth=0.5, color='k')

        plt.scatter(x=currTrial['avg x pos'], y=-currTrial['avg y pos'] + SCREENHEIGHT, color='lightslategray', s=currTrial['duration'] * scale,
                    alpha=0.5, edgecolor='none')
        
        axes = plt.gca()
        axes.set_title('Trial ' + str(trial))
        axes.set_xlim(0, SCREENWIDTH)
        axes.set_ylim(0, SCREENHEIGHT)

        iarea_name_list = allTrialsIAREA.loc[(allTrialsIAREA['trialN'] == trial)]['id'].values
        for cName in iarea_name_list:
            currIAREA = allTrialsIAREA.loc[(allTrialsIAREA['trialN'] == trial) & (allTrialsIAREA['id'] == cName)]
            n = currIAREA.index.values[0]
            leftX = currIAREA.at[(n, 'LeftX')]
            lowerY = - currIAREA.at[(n, 'LowerY')] + SCREENHEIGHT
            rightX = currIAREA.at[(n, 'RightX')]
            upperY = - currIAREA.at[(n, 'UpperY')] + SCREENHEIGHT

            iareaRect = ptch.Rectangle((leftX, lowerY), rightX - leftX, upperY - lowerY, linewidth=1, edgecolor='r',
                                       facecolor='none')
            axes.add_patch(iareaRect)
        screenRect = ptch.Rectangle((0, 0), SCREENWIDTH, SCREENHEIGHT, linewidth=1, edgecolor='r', facecolor='none')
        axes.add_patch(screenRect)
        # based on https://jakevdp.github.io/PythonDataScienceHandbook/04.06-customizing-legends.html
        for duration in [50*scale, 100*scale, 200*scale]:
            plt.scatter([], [], c='lightslategray', alpha=0.3, s=duration, edgecolor='none',
                        label=str(duration/scale) + ' ms')
        plt.legend(scatterpoints=1, frameon=False, labelspacing=2, bbox_to_anchor=(1, 1), title='Fixation duration')
        plt.show()

        
def plot_all_fixations(title, eyedata):
    allTrialsFixationData = eyedata.fixationData
    allTrialsIAREA = eyedata.iareaData
    SCREENWIDTH, SCREENHEIGHT = get_screen_width_height(eyedata.displayCoords)
    scale = 3
    
    print("Plotting fixations now...")
    trialNums=np.unique(allTrialsFixationData['trialN'].values).astype(int)
    
    plt.figure(title)
    axes = plt.gca()
    plt.plot([SCREENWIDTH/2,SCREENWIDTH/2], [0,SCREENHEIGHT], linewidth = 0.5, color = 'k')
    plt.plot([0,SCREENWIDTH], [SCREENHEIGHT/2,SCREENHEIGHT/2], linewidth = 0.5, color = 'k')

    for trial in trialNums:
        if trialNums[0] == 0:
            i=trial+1
        else:
            i=trial
        currTrial = allTrialsFixationData.loc[allTrialsFixationData['trialN']==trial]
        
        fig = plt.scatter(x=currTrial['avg x pos'], y=-currTrial['avg y pos']+SCREENHEIGHT, s=currTrial['duration']*scale, alpha = 0.5, edgecolor = 'none');
        
    axes.set_title(title)
    axes.set_xlim(0,SCREENWIDTH)
    axes.set_ylim(0,SCREENHEIGHT)

    iarea_name_list=allTrialsIAREA.loc[(allTrialsIAREA['trialN']==trial)]['id'].values
    for cName in iarea_name_list:
        currIAREA = allTrialsIAREA.loc[(allTrialsIAREA['trialN']==trial) & (allTrialsIAREA['id']==cName)]
        n=currIAREA.index.values[0]
        leftX = currIAREA.at[(n,'LeftX')]
        lowerY = - currIAREA.at[(n,'LowerY')] + SCREENHEIGHT
        rightX = currIAREA.at[(n,'RightX')]
        upperY = - currIAREA.at[(n,'UpperY')] + SCREENHEIGHT

        iareaRect = ptch.Rectangle((leftX,lowerY),rightX-leftX,upperY-lowerY,linewidth=1,edgecolor='r',facecolor='none')
        axes.add_patch(iareaRect)
    screenRect = ptch.Rectangle((0,0),SCREENWIDTH,SCREENHEIGHT,linewidth=1,edgecolor='r',facecolor='none')
    axes.add_patch(screenRect)
    
    # based on https://jakevdp.github.io/PythonDataScienceHandbook/04.06-customizing-legends.html
    for duration in [50*scale, 100*scale, 200*scale]:
        plt.scatter([], [], c='gray', alpha=0.3, s=duration, edgecolor='none',
                    label=str(duration/scale) + ' ms')
    plt.legend(scatterpoints=1, frameon=False, labelspacing=2, bbox_to_anchor=(1, 1), title='Fixation duration')
    
    plt.show()
    

def plot_ordered_trial_fixations(eyedata):
    
    """
    If using interactive python, can plot each fixation and the interest areas
    Also plot arrows in between fixations
    """
    
    allTrialsFixationData = eyedata.fixationData
    allTrialsIAREA = eyedata.iareaData
    displaycoords = eyedata.displayCoords
    SCREENWIDTH, SCREENHEIGHT = get_screen_width_height(eyedata.displayCoords)
    scale = 3
    
    print("Plotting fixations now...")
    trialNums = np.unique(allTrialsFixationData['trialN'].values).astype(int)
    for trial in trialNums:
        if trialNums[0] == 0:
            i = trial + 1
        else:
            i = trial
        currTrial = allTrialsFixationData.loc[allTrialsFixationData['trialN'] == trial]
        fig = plt.figure(i)
        plt.plot([SCREENWIDTH / 2, SCREENWIDTH / 2], [0, SCREENHEIGHT], linewidth=0.5, color='k')
        plt.plot([0, SCREENWIDTH], [SCREENHEIGHT / 2, SCREENHEIGHT / 2], linewidth=0.5, color='k')
        
        plt.scatter(x=currTrial['avg x pos'], y=-currTrial['avg y pos'] + SCREENHEIGHT, c="lightslategray", s=currTrial['duration'] * scale,
                    alpha=0.5, edgecolor='none')
        
        axes = plt.gca()
        axes.set_title('Trial ' + str(trial))
        axes.set_xlim(0, SCREENWIDTH)
        axes.set_ylim(0, SCREENHEIGHT)
        
        colors = mpl.cm.get_cmap('Blues', len(currTrial))
        # plot arrows in between fixations
        xvals = currTrial['avg x pos'].values
        yvals = currTrial['avg y pos'].values
        for j in range(1,len(xvals)):
            x1 = xvals[j-1]
            y1 = yvals[j-1]
            x2 = xvals[j]
            y2 = yvals[j]
            axes.arrow(x1, -y1+SCREENHEIGHT, x2-x1, y1-y2, head_width=20, head_length=20, fc=colors(j), ec=colors(j))

        iarea_name_list = allTrialsIAREA.loc[(allTrialsIAREA['trialN'] == trial)]['id'].values
        for cName in iarea_name_list:
            currIAREA = allTrialsIAREA.loc[(allTrialsIAREA['trialN'] == trial) & (allTrialsIAREA['id'] == cName)]
            n = currIAREA.index.values[0]
            leftX = currIAREA.at[(n, 'LeftX')]
            lowerY = - currIAREA.at[(n, 'LowerY')] + SCREENHEIGHT
            rightX = currIAREA.at[(n, 'RightX')]
            upperY = - currIAREA.at[(n, 'UpperY')] + SCREENHEIGHT

            iareaRect = ptch.Rectangle((leftX, lowerY), rightX - leftX, upperY - lowerY, linewidth=1, edgecolor='r',
                                       facecolor='none')
            axes.add_patch(iareaRect)
        screenRect = ptch.Rectangle((0, 0), SCREENWIDTH, SCREENHEIGHT, linewidth=1, edgecolor='r', facecolor='none')
        axes.add_patch(screenRect)
        
        
        # based on https://jakevdp.github.io/PythonDataScienceHandbook/04.06-customizing-legends.html
        for duration in [50*scale, 100*scale, 200*scale]:
            plt.scatter([], [], c='lightslategray', alpha=0.3, s=duration, edgecolor='none',
                        label=str(round(duration/scale)) + ' ms')
        plt.legend(scatterpoints=1, frameon=False, labelspacing=2, bbox_to_anchor=(1, 1), title='Fixation duration')
        
        # https://matplotlib.org/3.3.3/gallery/axes_grid1/demo_colorbar_with_axes_divider.html#sphx-glr-gallery-axes-grid1-demo-colorbar-with-axes-divider-py
        # https://matplotlib.org/3.1.0/gallery/axes_grid1/demo_colorbar_with_inset_locator.html
        cmap = mpl.cm.get_cmap('Blues', len(currTrial))
        axins = inset_axes(axes,
                   width="25%",  # width = 5% of parent_bbox width
                   height="10%",  # height : 50%
                   bbox_to_anchor=(0.3,-0.7,1,1),
                   bbox_transform=axes.transAxes,
                   borderpad=0,
                   )
        
        bounds = range(1,len(currTrial)+1)
        cb = fig.colorbar(mpl.cm.ScalarMappable(cmap=colors), cax=axins, boundaries=bounds, orientation="horizontal")
        cb.set_label('Fixation number')
        
        plt.show()
    
def write_fixation_data_to_csv(asc_name, csv_name, trial_marker_start='TRIALID', trial_marker_end='TRIAL_RESULT'):
    """
    Gets fixation data from process_csv()
    Summarizes fixation data and stores in csv
    """
    if asc_name.endswith('.asc') and os.path.exists(asc_name):

        allTrialsFixationData, allTrialsIAREA, displaycoords = process_edf(asc_name, trial_marker_start,
                                                                           trial_marker_end)
        trialStartTimes, trialEndTimes = get_trial_times(asc_name, trial_marker_start, trial_marker_end)

        allTrialsFixationData.to_csv(csv_name[:-4] + '_all_fixations.csv')
        trialTimes = pd.DataFrame(trialStartTimes.items(), columns=['trialN', 'startTime'])
        trialTimes['endTime'] = trialTimes['trialN'].map(trialEndTimes)
        trialTimes['duration'] = trialTimes['endTime'] - trialTimes['startTime']
        df = pd.DataFrame()
        df['trialN'] = np.unique(allTrialsFixationData['trialN']).astype(int)

        iarea_name_list = np.unique(allTrialsIAREA['id'])
        for iarea in iarea_name_list:
            if 'fix on ' + iarea in allTrialsFixationData.columns:
                df[iarea + 'TotalFixation'] = trialTimes['trialN'].map(
                    allTrialsFixationData.loc[allTrialsFixationData['fix on ' + iarea] == True].groupby('trialN')[
                        'duration'].sum())
                df[iarea + 'FixationCount'] = trialTimes['trialN'].map(
                    allTrialsFixationData.loc[allTrialsFixationData['fix on ' + iarea] == True].groupby('trialN')[
                        'duration'].count())
                df.fillna(0, inplace=True)
                df[iarea + 'FixationCount'] = df[iarea + 'FixationCount'].values.astype(int)
        df.to_csv(csv_name, index=False)
        print("Data written to %s." % (csv_name))
    else:
        print("File not found: %s" % (asc_name))


def check_fixation_in_bounds(iareaType, fixX, fixY, leftX, rightX, lowerY, upperY):
    """
    Determine whether fixation is within interest area
    """
    
    # in case leftX and rightX were swapped
    if leftX > rightX:
        tmp = leftX
        leftX = rightX
        rightX = tmp
        
    # in case upperY and lowerY were swapped
    if lowerY > upperY:
        tmp = lowerY
        lowerY = upperY
        upperY = tmp
        
    b1 = leftX <= fixX <= rightX # within left and right bounds of rectangle interest area
    b2 = upperY >= fixY >= lowerY # within upper and lower bounds of rectangle interest area
    if b1 and b2:
        if iareaType == 'ELLIPSE':
            return check_fixation_in_ellipse(fixX, fixY, leftX, rightX, lowerY, upperY)
        return True
    return False


def check_fixation_in_ellipse(fixX, fixY, leftX, rightX, lowerY, upperY):
    centerX = (rightX - leftX) / 2 + rightX
    centerY = (upperY - lowerY) / 2 + lowerY
    a = (rightX - leftX) / 2
    b = (upperY - lowerY) / 2
    if ((fixX - centerX) ** 2) / (a ** 2) + ((fixY - centerY) ** 2) / (b ** 2) <= 1:
        return True
    return False
