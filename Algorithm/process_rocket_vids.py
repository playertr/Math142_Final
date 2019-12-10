#!/usr/bin/env python3
'''
    File name: process_rocket_vids.py

    Author: Tim Player
    Date created: 04/28/2019
    Python Version: 3

    Note: debug functionality through pdb.
'''

from imutils.video import VideoStream
import imutils
import cv2
import numpy as np
import pdb
import os
from functools import reduce
import csv

def main():
    """
    Selects two camera videos and runs through them, identifying the horizon and printing it to an output file.
    Prints the horizons from both cameras to the output-horizon file. Note: it will probably overwrite the file if you run it for both flights. The format for this file is, as defined at the end of the main loop, [fnum, x0, y0, x1, y1, x0, y0, x1, y1], where the first two points are from Cam1 and the second two points are from Cam2. 
    """
    # designate filestrings depending on which flight we are looking at
    flight = 'flight2'
    if flight == 'flight3':
        cam1filestring = r"Cam1_third_launch/MINI0001.AVI"
        cam2filestring = r"Cam2_third_launch/MINI0001.AVI"
        outputfolderstring = r"/Flight3_launch"

        cam1startframe = 16379  # roughly 9 minutes in
        cam2startframe = 15726

    elif flight == 'flight2':
        cam1filestring = r"Cam1_second_launch/MINI0010.AVI"
        cam2filestring = r"Cam2_second_launch/MINI0009.AVI"

        cam1startframe = 24271
        cam2startframe = 24422

    # navigate to video files
    fileDir = os.path.dirname(os.path.realpath('__file__'))
    cam1path = os.path.join(fileDir, cam1filestring)
    cam2path = os.path.join(fileDir, cam2filestring)

    # play video
    vs1 = cv2.VideoCapture(cam1path)
    vs2 = cv2.VideoCapture(cam2path)

    # fast forward to frame where you start
    vs1.set(1, cam1startframe-1)
    vs2.set(1, cam2startframe-1)

    # define VideoWriter
    fname = os.path.join(fileDir, 'output.avi')

    global width, height
    width = int(vs1.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(vs1.get(cv2.CAP_PROP_FRAME_HEIGHT))
    fps = 30
    writer = cv2.VideoWriter(
        filename=fname,
        fourcc=(ord('X') << 24) + (ord('V') << 16) +
        (ord('I') << 8) + ord('D'),
        fps=fps,
        frameSize=(width, height),
        isColor=1)

    # define horizon output file
    csvFile = open(fname.rpartition(".")[0] + "-horizon.txt", "w")
    csvWriter = csv.writer(csvFile)

    # Initialize windows
    windowNames = ['Cam1', 'Cam2']#, 'Fitted Line', 'Mask', 'M1', 'M2', 'M3', 'M4', 'Edges', 'Hcont']
    for wn in windowNames:
        cv2.namedWindow(wn, cv2.WINDOW_NORMAL)
        cv2.resizeWindow(wn, (600, 350))

    # Move relevant windows to forefront again
    cv2.moveWindow('Cam2', 300, 0)

    # Image processing loop
    fnum = 0
    while(vs1.isOpened() and vs2.isOpened()):
        fnum += 1

        # Print time and frame number
        time = fnum / 30.0
        c1frame = fnum + cam1startframe
        c2frame = fnum + cam2startframe
        #print(f'Time: {time: .4} \tCam1 Frame: {c1frame} \tCam2 Frame: {c2frame}')

        # Read frames
        ret, frame1 = vs1.read()
        ret, frame2 = vs2.read()
        
        #prepare processing parameters
        newwidth = int(width/4)
        newheight = int(height/4)

        # Define HSV colors of different scene items

        #clouds
        colorLower1=(0, 0, 179)
        colorHigher1=(82, 255, 255)

        #blue sky
        colorLower2=(2, 0, 153)
        colorHigher2=(19, 255, 252)

        # mountains
        colorLower3 = (9, 0, 115)
        colorHigher3 = (103, 78, 203)

        # sand
        colorLower4 = (77, 38, 163)
        colorHigher4 = (122, 100, 255)

        includeColors = [ (colorLower1, colorHigher1), (colorLower2, colorHigher2)]

        excludeColors = [ (colorLower3, colorHigher3), (colorLower4, colorHigher4)]

        # Identify and print horizons
        horiz1 = getHorizon(frame1, 'Cam1', includeColors, excludeColors, newwidth, newheight)

        horiz2 = getHorizon(frame2, 'Cam2', includeColors, excludeColors, newwidth, newheight)

        #pitch and yaw from:
        #https://www.grc.nasa.gov/WWW/K-12/rocket/rotations.html
        yaw = getYaw(horiz1, horiz2, 0)
        pitch = getPitch(horiz1, horiz2, 0)

        #write to csv
        row = [fnum, horiz1[0], horiz1[1], horiz1[2], horiz1[3], horiz2[0], horiz2[1], horiz2[2], horiz2[3]]

        csvWriter.writerow(row)

        key = cv2.waitKey(33) & 0xFF
        if key == ord('q'):
            break
        elif key == ord('p'):
            pdb.set_trace()
    

    vs1.release()
    vs2.release()
    csvFile.close()
    writer.release()
    cv2.destroyAllWindows()

def getYaw(horiz1, horiz2, expected_roll):
    """
    Returns the yaw in radians.
    
    Arguments:
        horiz1 {Tuple (x1, y1, x2, y2)} -- two points on the horizon of cam1
        horiz2 {Tuple (x1, y1, x2, y2)} -- two points on the horizon of cam2
        expected_roll {float} -- roll in radians (perhaps from magnetometer?)
    
    Returns:
        [float] -- yaw in radians as defined in https://www.grc.nasa.gov/WWW/K-12/rocket/rotations.html, where the pitch axis is out camera1.
    """
    # I have not implemented this yet.
    return None

def getPitch(horiz1, horiz2, expected_roll):
    """
    Returns the pitch in radians.
    
    Arguments:
        horiz1 {Tuple (x1, y1, x2, y2)} -- two points on the horizon of cam1
        horiz2 {Tuple (x1, y1, x2, y2)} -- two points on the horizon of cam2
        expected_roll {float} -- roll in radians (perhaps from magnetometer?)
    
    Returns:
        [float] -- pitch in radians as defined in https://www.grc.nasa.gov/WWW/K-12/rocket/rotations.html, where the pitch axis is out camera1.
    """
    # I have not implemented this yet.
    return None

def getMask(img, colorLower, colorHigher):
    """
    Finds the areas of an image falling between the color parameters colorLower and colorHigher. Pefforms erosion and dilation ("morphological restructuring") to make more uniform.
    
    Arguments:
        img {np.ndarray} -- the image of interest. Probably hsv.
        colorLower {Tuple (h, s, v)} -- hsv lower bound
        colorHigher {Tuple (h, s, v)} -- hsv higher bound
    
    Returns:
        [type] -- [description]
    """
    mask = cv2.inRange(img, colorLower, colorHigher)
    mask = cv2.erode(mask, None, iterations=2)
    mask = cv2.dilate(mask, None, iterations=2)
    return mask

def getHorizon(img, outputWindow, includeColors, excludeColors, newwidth, newheight, debug=False):
    """
    Finds the horizon, given two lists of colors that identify the sky and land. Returns two points on that horizon. Internally, the algorithm is roughly:

    1) make a mask that includes the sky but not the ground.
    2) find the edges of that mask
    3) pick the longest edge. Assume that contains the horizon.
    4) crop the image to only examine the central rectangle (and not the edges of the circular window)
    5) use Hough line transform to identify long straight line
    
    Arguments:
        img {np.ndarray} -- rgb image of interest
        outputWindow {String} -- Name of window to output to
        includeColors {List of Tuples } -- hsv colors to include in sky
        excludeColors {List of Tuples} -- hsv colors to exclude from sky
        newwidth {int} -- width to resize to
        newheight {int} -- height to resize to
    
    Keyword Arguments:
        debug {bool} -- Specifies whether to display intermediate steps of algorithm (default: {False})
    
    Returns:
        Tuple of Four Floats -- (x1, y1, x2, y2) coordinates of points that are on the identified horizon line
    """

    # Downsample image
    img = cv2.resize(img, (newwidth, newheight), 0, 0, cv2.INTER_NEAREST)

    # Create masks based on HSV values in image
    hsv=cv2.cvtColor(img, cv2.COLOR_RGB2HSV)

    #construct masks
    includeMasks = []
    for (lower, upper) in includeColors:
        includeMasks += [getMask(hsv, lower, upper)]
    
    excludeMasks = []
    for (lower, upper) in excludeColors:
        excludeMasks += [getMask(hsv, lower, upper)]
    
    # create a mask that includes cloud and blue sky but not sand or mountain
    mask = reduce(lambda x,y: x | y , includeMasks)
    mask = reduce(lambda x,y: x & ~y, excludeMasks, mask)

    # blur the mask for good measure
    mask=cv2.GaussianBlur(mask, (25, 25), 0)

    # Apply binary threshold for edgefinding
    ret, mask=cv2.threshold(mask, 127, 255, cv2.THRESH_BINARY)

    #extract the edges of the mask
    edges=auto_canny(mask)

    #retain only the largest contour
    cnts = cv2.findContours(edges, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    contours = imutils.grab_contours(cnts)

    # Create blank "horizon contour" image
    hcont = np.zeros((newheight,newwidth), np.uint8) 
    if len(contours) != 0:
        #find contour with max perimeter
        c= max(contours, key=lambda x: cv2.arcLength(x, False))

        #draw the longest contour to the horizon contour
        cv2.drawContours(hcont, c, -1, 255, 3, maxLevel=1)
    else:
        print("No contours found.")

    # Only look at the contour in the center edge
    centermask = np.zeros((newheight,newwidth), np.uint8)
    center = (int(newwidth/2), int(newheight/2))
    upperLeft = (int(center[0] - newwidth/4), int(center[1] - newheight/4))
    lowerRight = (int(center[0] + newwidth/4), int(center[1] + newheight/4))
    cv2.rectangle(centermask, upperLeft, lowerRight, 255, -1)
    hcont = hcont & centermask

    # Use Hough Line algorithm to find straight lines that are long enough
    minLineLength = 100
    maxLineGap = 200
    lines = cv2.HoughLinesP(hcont,rho=1,theta=np.pi/180,threshold=10,minLineLength=minLineLength,maxLineGap=maxLineGap)

    # Draw the fitted horizon on a copy of image
    imcopy=img.copy()
    if lines is not None:
        for x1, y1, x2, y2 in lines[0]:
            cv2.line(imcopy,(x1,y1),(x2,y2),(0,0,255),4)
            
    else:
        print("No Hough lines found")
        x1, y1, x2, y2 = (-1, -1, -1, -1)
    
    cv2.imshow(outputWindow, imcopy)
    #display all the images for debugging
    if debug:
        # Initialize windows
        windowNames = ['Mask', 'M1', 'M2', 'M3', 'M4', 'Edges', 'Hcont']
        for wn in windowNames:
            cv2.namedWindow(wn, cv2.WINDOW_NORMAL)
            cv2.resizeWindow(wn, (600, 350))

        cv2.imshow('Mask', mask)
        cv2.imshow('M1', includeMasks[0])
        cv2.imshow('M2', includeMasks[1])
        cv2.imshow('M3', excludeMasks[0])
        cv2.imshow('M4', excludeMasks[1])
        cv2.imshow('Edges', edges)
        cv2.imshow('Hcont', hcont)

    return (x1, y1, x2, y2)


# https://www.pyimagesearch.com/2015/04/06/zero-parameter-automatic-canny-edge-detection-with-python-and-opencv/
def auto_canny(image, sigma=0.33):
    """
    Performs Canny edge detection on an image, automatically selecting the upper and lower threshold from the median values.
    
    Arguments:
        image {np.ndarray} -- image of interest. RGB, for instance.
    
    Keyword Arguments: Not sure what this is (default: {0.33})
    
    Returns:
        np.ndarray -- binary image with edges drawn I believe, you might want to check the Canny specs
    """
    # compute the median of the single channel pixel intensities
    v=np.median(image)

    # apply automatic Canny edge detection using the computed median
    lower=int(max(0, (1.0 - sigma) * v))
    upper=int(min(255, (1.0 + sigma) * v))
    edged=cv2.Canny(image, lower, upper, apertureSize=3)

    # return the edged image
    return edged

# Deprecated
def trimCenter(img, width, height, denom):
    """
    Helper function to crop an image.
    """
    startX = int(width/denom)
    endX = width - startX
    startY = int(height/denom)
    endY = height - startY
    return img[startY:endY, startX:endX]

# Needed for Trackbar callback
def nothing(a):
    return

# Fiddle with values
def adjustValues(img, colorLower, colorHigher):
    """
    Simple GUI to pick HSV values for a given mask. Has sliders where you can adjust, in convenient alphabetical order, H, S, and V for lower and upper bounds. Continually loops through, updating and displaying the image mask. Indispensable for picking out colors.
    
    Arguments:
        img {np.ndarray} -- rgb image of interest
        colorLower {3-Tuple} -- initial hsv lower bound 
        colorHigher {3-Tuple} -- initial hsb higher bound
    """

    # Make and resize window
    cv2.namedWindow('TB', cv2.WINDOW_NORMAL)
    cv2.resizeWindow('TB', (350, 350))

    # Set up trackbars
    # Note: they will be alphabetized
    cv2.namedWindow('TB', cv2.WINDOW_NORMAL)
    cv2.createTrackbar('HLower', 'TB', 0, 255, nothing)
    cv2.createTrackbar('SLower', 'TB', 0, 255, nothing)
    cv2.createTrackbar('VLower', 'TB', 0, 255, nothing)
    cv2.createTrackbar('HHigher', 'TB', 0, 255, nothing)
    cv2.createTrackbar('SHigher', 'TB', 0, 255, nothing)
    cv2.createTrackbar('VHigher', 'TB', 0, 255, nothing)

    #initialize trackbar positions to function arguments
    hL=cv2.setTrackbarPos('HLower', 'TB', colorLower[0])
    sL=cv2.setTrackbarPos('SLower', 'TB', colorLower[1])
    vL=cv2.setTrackbarPos('VLower', 'TB', colorLower[2])
    hH=cv2.setTrackbarPos('HHigher', 'TB', colorHigher[0])
    sH=cv2.setTrackbarPos('SHigher', 'TB', colorHigher[1])
    vH=cv2.setTrackbarPos('VHigher', 'TB', colorHigher[2])

    #Continually monitor trackbar position and re-run mask creation
    while(True):
        hL=cv2.getTrackbarPos('HLower', 'TB')
        sL=cv2.getTrackbarPos('SLower', 'TB')
        vL=cv2.getTrackbarPos('VLower', 'TB')
        hH=cv2.getTrackbarPos('HHigher', 'TB')
        sH=cv2.getTrackbarPos('SHigher', 'TB')
        vH=cv2.getTrackbarPos('VHigher', 'TB')

        img=cv2.GaussianBlur(img, (3, 3), 0)

        # do hsv thresholding
        blurred=cv2.GaussianBlur(img, (3, 3), 0)
        hsv=cv2.cvtColor(blurred, cv2.COLOR_RGB2HSV)

        colorLower=(hL, sL, vL)
        colorHigher=(hH, sH, vH)

        print("colorLower: ", colorLower)
        print("colorHigher: ", colorHigher)

        mask = getMask(hsv, colorLower, colorHigher)
        cv2.imshow('Mask', mask & img)

        key=cv2.waitKey(33) & 0xFF
        if key == ord('q'):
            quit()
        elif key == ord('p'):
            pdb.set_trace()
        
main()
