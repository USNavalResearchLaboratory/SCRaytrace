"""
mathutil
========

Provide various useful functions and tools.

.. module:: mathutil
    :synopsis: Provides various useful functions and tools.

"""

import numpy as np
from scipy.ndimage import map_coordinates
import os
import glob
import copy
import ntpath
from scipy import ndimage as ndi
import collections
import matplotlib as mpl
import platform
if platform == "Darwin" :
    mpl.use('TkAgg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from astropy.wcs import WCS
from astropy.io import fits
import pathlib
import yaml
import pkgutil

# from PySide2.QtWidgets import (QWidget, QPushButton, 
#     QHBoxLayout, QVBoxLayout, QInputDialog)

if pkgutil.find_loader('PyQt5') is not None:
    from PyQt5.QtWidgets import (QWidget, QPushButton, 
                                 QHBoxLayout, QVBoxLayout, QInputDialog)
elif pkgutil.find_loader('PyQt6') is not None:
    from PyQt6.QtWidgets import (QWidget, QPushButton, 
                                 QHBoxLayout, QVBoxLayout, QInputDialog)
else:
    print("Could not find PyQt5 ot PyQt6")


def listFilesByDate(fileMask, extractFilename=True, return_pathlib=False):
    """List the files in a directory and order them by date."""
    lst = glob.glob(fileMask)
    lst.sort(key=os.path.getmtime)
    
    if extractFilename:
        lstFilename = []
        for i in lst:
            lstFilename.append(ntpath.split(i)[1])
        lst = lstFilename
    
    if return_pathlib :
        for i in range(len(lst)):
            lst[i] = pathlib.Path(lst[i])
    
    return lst


def readfits(b):
    hdu = fits.open(b)
    return hdu


def index_coords(data, origin=None):
    """Creates x & y coords for the indicies in a numpy array "data".
    "origin" defaults to the center of the image. Specify origin=(0,0)
    to set the origin to the lower left corner of the image."""
    ny, nx = data.shape[:2]
    if origin is None:
        origin_x, origin_y = nx // 2, ny // 2
    else:
        origin_x, origin_y = origin
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    x -= origin_x
    y -= origin_y
    return x, y

def cart2polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta

def polar2cart(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y


def reproject_image_into_polar(data, origin=None, nbpixRadial=None, nbpixAngle=None):
    """Reprojects a 3D numpy array ("data") into a polar coordinate system.
    "origin" is a tuple of (x0, y0) and defaults to the center of the image."""
    ny, nx = data.shape[:2]
    if origin is None:
        origin = (nx//2, ny//2)

    # Determine that the min and max r and theta coords will be...
    x, y = index_coords(data, origin=origin)
    r, theta = cart2polar(x, y)

    # Set the output image size
    if nbpixRadial is None:
        nbpixRadial = nx
    if nbpixAngle is None:
        nbpixAngle = ny


    # Make a regular (in polar space) grid based on the min and max r & theta
    r_i = np.linspace(r.min(), r.max(), nbpixRadial)
    theta_i = np.linspace(theta.min(), theta.max(), nbpixAngle)
    theta_grid, r_grid = np.meshgrid(theta_i, r_i)

    # Project the r and theta grid back into pixel coordinates
    xi, yi = polar2cart(r_grid, theta_grid)
    xi += origin[0]  # We need to shift the origin back to 
    yi += origin[1]  # back to the lower-left corner...
    xi, yi = xi.flatten(), yi.flatten()
    coords = np.vstack((xi, yi))  # (map_coordinates requires a 2xn array)

    # Reproject each band individually and the restack
    # (uses less memory than reprojection the 3-dimensional array in one step)
    
#    bands = []
#    zi = sp.ndimage.map_coordinates(data, coords, order=1)
    zi = map_coordinates(data, coords, order=1)
#    bands.append(zi.reshape((nx, ny)))
    output = zi.reshape((nbpixRadial, nbpixAngle))
    
#    for band in data.T:
#        zi = sp.ndimage.map_coordinates(band, coords, order=1)
#        bands.append(zi.reshape((nx, ny)))
#    output = np.dstack(bands)
    return output, r_i, theta_i






def rotmat(ang_rad, axis):
    """Compute a rotation matrix
    
    Args:
        ang_rad (float) : Angle in radian.
        axis : (1, 2, 3) or ('x', 'y', 'z')
        retArray (bool) : set to True to return an np.array(), else returns a np.matrix()
        
    Returns:
        3D rotation matrix
    
    """
    if axis == 1 or axis == 'x' or axis == 'X' :
        rmat = np.matrix([[1, 0, 0],
                         [0, np.cos(ang_rad), -np.sin(ang_rad)],
                         [0, np.sin(ang_rad),  np.cos(ang_rad)]])
    elif axis == 2 or axis =='y' or axis == 'Y' :
        rmat = np.matrix([[np.cos(ang_rad), 0, np.sin(ang_rad)],
                         [0, 1, 0],
                         [-np.sin(ang_rad), 0,  np.cos(ang_rad)]])
    elif axis == 3 or axis =='z' or axis == 'Z' :
        rmat = np.matrix([[np.cos(ang_rad), -np.sin(ang_rad), 0],
                         [np.sin(ang_rad), np.cos(ang_rad), 0],
                         [0, 0, 1]])
    else: rmat = -1

    return rmat


def rmat(ang_rad, axis):
    if axis == 1 or axis == 'x' or axis == 'X' :
        rmat = np.array([[1, 0, 0],
                         [0, np.cos(ang_rad), -np.sin(ang_rad)],
                         [0, np.sin(ang_rad),  np.cos(ang_rad)]])
    elif axis == 2 or axis =='y' or axis == 'Y' :
        rmat = np.array([[np.cos(ang_rad), 0, np.sin(ang_rad)],
                         [0, 1, 0],
                         [-np.sin(ang_rad), 0,  np.cos(ang_rad)]])
    elif axis == 3 or axis =='z' or axis == 'Z' :
        rmat = np.array([[np.cos(ang_rad), -np.sin(ang_rad), 0],
                         [np.sin(ang_rad), np.cos(ang_rad), 0],
                         [0, 0, 1]])
    else: rmat = -1

    return rmat



def dispim(im, title='', log=False, minmax=None, 
           cmap='gist_ncar', figsize=(12, 9), units='', 
           windowXYPos=(50, 50), origin='lower',
           extent=None, show=True,
           axin=None, axislabels=True):
    """Display an image with a color bar"""
    
    # palette = copy.copy(matplotlib.cm.get_cmap(cmap))
    palette = copy.copy(mpl.colormaps.get_cmap(cmap))
    palette.set_under(color='black')
    palette.set_bad(color='black')

    if axin is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    else:
        fig = None
        ax = axin
    
    if minmax is None:
        minmax = (np.min(im), np.max(im))
    
    logUnitPrefix = ''
    if log:
        norm = mpl.colors.LogNorm(vmin=minmax[0], vmax=minmax[1])
        logUnitPrefix = 'Log'
    else:
        norm = mpl.colors.Normalize(vmin=minmax[0], vmax=minmax[1])
        logUnitPrefix = ''

    pim = ax.imshow(im, origin=origin, cmap=palette, norm=norm, extent=extent)
    ax.set_title(title)
    if axin is None:
        cb = fig.colorbar(pim)
        cb.set_label("{0} {1}".format(logUnitPrefix, units))
    
    if axislabels:
        ax.set_xlabel('x [pix]')
        ax.set_ylabel('y [pix]')

    if(show == True):
        plt.show()
        
    plt.tight_layout()

    
    return fig, ax



def dispfits(im, header, log=True, title='', units='',
           minmax=(1, 16000), cmap='gist_ncar', 
           figsize=(12, 9), windowXYPos=(50, 50),
           key=None, wcsHead=None, plotCRPIX=False):
    """Display fits image with grid overlay based on WCS header"""
    
    c = copy.copy(mpl.colormaps.get_cmap(cmap))
    c.set_under(color='black')
    c.set_bad(color='black')

    if wcsHead is None:
        wcsHead = WCS(header, key=key)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=wcsHead)
    
#    logUnitPrefix = ''
    if log:
        norm = mpl.colors.LogNorm(vmin=minmax[0], vmax=minmax[1])
        logUnitPrefix = 'Log'
    else:
        norm = mpl.colors.Normalize(vmin=minmax[0], vmax=minmax[1])
        logUnitPrefix = ''

    pim = ax.imshow(im, origin='lower', cmap=c, norm=norm)
    ax.set_title(title)
    
    cb = fig.colorbar(pim)
    cb.set_label("{0} {1}".format(logUnitPrefix, units))

    ax.coords.grid(color='white')
    ax.coords[0].set_format_unit('deg')
    ax.coords[1].set_format_unit('deg')

    # ax.coords['HPLN'].set_axislabel('LN')
    # ax.coords['HPLT'].set_axislabel('LT')
    ax.set_title(title)
    if plotCRPIX:
        crpix = wcsHead.wcs.crpix
        ax.plot(crpix[0]-1, crpix[1]-1, marker='x')
    
       
    plt.show()

    thismanager = plt.get_current_fig_manager()
    
    geom = thismanager.window.geometry()
    x,y,dx,dy = geom.getRect()
    
    thismanager.window.setGeometry(windowXYPos[0], windowXYPos[1], dx, dy)

    return fig, ax


def dispFits(im, header, log=True, title='', units='',
           minmax=(1, 16000), cmap='gist_ncar', 
           figsize=(12, 9), windowXYPos=(50, 50),
           key=None, wcsHead=None, plotCRPIX=True):
    """Wraper for dispfits. No mo camel case."""
    fig, ax = dispfits(im, header, log=log, title=title, units=units,
                       minmax=minmax, cmap=cmap, 
                       figsize=figsize, windowXYPos=windowXYPos,
                       key=key, wcsHead=wcsHead, plotCRPIX=plotCRPIX)
    return fig, ax
    

def dispMosaic(imList, extent=None, 
               minmax=(1e1, 1e4), 
               log=True, 
               grid = (9, 10), 
               cmap='gray',
               titleLst=None):
    """Display the image mosaic"""

    fig = plt.figure(figsize=(17, 11))

    if log:
        norm = mpl.colors.LogNorm(vmin=minmax[0], vmax=minmax[1])
    else:
        norm = mpl.colors.Normalize(vmin=minmax[0], vmax=minmax[1])

    for i in range(len(imList)):
        ax = fig.add_subplot(grid[0], grid[1], i+1)
        ax.tick_params(axis='both', which='major', labelsize=8)
        
        ax.imshow(imList[i], origin='lower', 
                        cmap=cmap, norm=norm, extent=extent)
        if titleLst is None:
            ax.set_title("ID:{0}".format(i), size=8.)
        else:
            ax.set_title(titleLst[i], size=8.)
        
    fig.subplots_adjust(top=0.96, left=0.07, right=0.97, hspace=0.41, wspace=0.3, bottom=0.03)
    plt.show()

    return fig, ax



class Line2D():
    """Defines a 2D line"""
    def __init__(self, point, vector):
        self.point = point
        self.vector = vector
        # -- ensure that the vector is unit length
        self.vector = self.vector / np.linalg.norm(self.vector)
        
        


def lineIntersection(l1, l2):
    """Find the intersection point between two lines"""
    a = np.array([-l1.vector, l2.vector]).T
    b = l1.point - l2.point
    r = np.linalg.solve(a, b)
    
    p1 = l1.point + r[0] * l1.vector
    p2 = l2.point + r[1] * l2.vector
    
    return p1, p2




def imageStat(im=None, mask=None, getColNames=False):
    """Returns various statistics in a region of an image"""
    
    colnames = ['med', 'avg', 'std', 'pop', 'min', 'max','xcen', 'ycen']
    if getColNames: return colnames
    
    if mask is not None:
        m = mask > 0
    else:
        m = np.full_like(im, True, dtype=bool)
   
    med = np.median(im[m])
    avg = np.mean(im[m])
    std = np.std(im[m])
    pop = np.count_nonzero(mask)
    minimum = np.min(im[m])
    maximum = np.max(im[m])
    xy=np.mgrid[0:m.shape[0],0:m.shape[1]]
    x=xy[1,:]
    y=xy[0,:]
    posx = (np.max(x[m]) + np.min(x[m]))/2.     #np.mean(x[m])
    posy = (np.max(y[m]) + np.min(y[m]))/2.     #np.mean(y[m])
    
    Stat = collections.namedtuple('Stat', colnames)
    
    s = Stat(med, avg, std, pop, minimum, maximum, posx, posy)
    
    return s



class defineRoi():
    def drawZone(self, closed=False):
        
        if self.art is not None:
            self.art.remove()
        
        polyg = mpatches.Polygon(self.points, edgecolor='red', linewidth=1., fill=False, closed=closed)
        self.art = self.ax.add_artist(polyg)
        
        plt.draw()
    
    
    def onclick(self, event):
        """What to do when button is clicked
        
        button 1: new point
        button 2: errase previous point
        button 3: start new ROI definition
        """

        
#        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % ('double' if event.dblclick else 'single', event.button, event.x, event.y, event.xdata, event.ydata))
        
#        print(event.button, event.x, event.y, event.xdata, event.ydata)

        if event.inaxes!=self.ax: return

        if event.button == 1 :
            # -- add new point
            if event.xdata is not None and event.ydata is not None:
#                print("inaxes:  ", event.inaxes)
                print("event.x: ", event.x)
                self.points.append(np.array([event.xdata, event.ydata]))

        elif event.button == 2 :
            # -- errase previous point
            if self.points:
                print('Removing last point')
                self.points.pop()
        else:
            # -- this should be button 3 so it's done
            print('Done with this ROI, starting a new one')
            # self.disconnect()
            self.done = True
        
        if not self.done:
            print('Nb Points :', len(self.points))
        
        if len(self.points) > 1:
            self.drawZone(closed=self.done)

        if self.done:
            self.ROIs.append(self.points)
            self.points = []
            self.art = None
            self.done = False



    def disconnect(self):
        print('Disconnecting.')
        self.fig.canvas.mpl_disconnect(self.cid)


    def __init__(self, fig, ax):
        """Define a polygon region of interest in an image"""
    
#        im = np.random.randint(0,255,(255,511))
#    
#        fig, ax = dispim(im)
        self.fig = fig 
        self.ax = ax
        
        self.ROIs = []
        
        self.points = []
        
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        
        self.done = False
        self.art = None
        
    def __del__(self):
        self.disconnect()






class DoneButton(QWidget):
    """Implement a simple done widget button with QT5"""
    
    def __init__(self, callFunc=None, param=None):
        super().__init__()
        
        self.initUI()
        
        self.callFunc = callFunc
        self.param = param
        
        self.clicked = False
        
        
    def initUI(self):
        
        doneButton = QPushButton("Done", self)
        doneButton.clicked.connect(self.buttonClicked)
        
        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(doneButton)

        vbox = QVBoxLayout()
        vbox.addStretch(1)
        vbox.addLayout(hbox)
        
        self.setLayout(vbox)    
        
        self.setGeometry(300, 300, 300, 150)
        self.setWindowTitle('Press when done')    
        self.show()
        
    def buttonClicked(self):
        if self.callFunc is not None:
            if self.param is not None: 
                self.output = self.callFunc(self.param)
            else:
                self.output = self.callFunc()
        
        self.clicked = True
        self.close()



class InputBox(QWidget):
    """Simple input box to interactively enter text."""    
    def __init__(self, message='Enter Something: '):
        super().__init__()
        self.message = message
        self.initUI()
        
    def initUI(self):      
        text, ok = QInputDialog.getText(self, 'Input Dialog', self.message)
        
        self.text = text
        self.ok = ok




def ellipseMask(mask, cntr, xyRadius, rotAngDeg=0.):
    """Compute an elliptical mask in an 2d image"""
    a = np.arange(0, 2 * np.pi, np.deg2rad(0.1))
    x0 = xyRadius[0] * np.cos(a)
    y0 = xyRadius[1] * np.sin(a)
    
    x = np.array(np.round( cntr[0] + (x0 * np.cos(np.deg2rad(rotAngDeg)) + y0 * np.sin(np.deg2rad(rotAngDeg)))), dtype=int)
    y = np.array(np.round( cntr[1] + (-x0 * np.sin(np.deg2rad(rotAngDeg)) + y0 * np.cos(np.deg2rad(rotAngDeg)))), dtype=int)
    
    m = (x > 1535)
    if any(m) : x[m] = 1535
    m = (x < 0)
    if any(m) : x[m] = 0
    m = (y > 1023)
    if any(m) : y[m] = 1023
    m = (y < 0)
    if any(m) : y[m] = 0
    
    miny = np.min(y)
    maxy = np.max(y)

    for yy in np.arange(miny, maxy):
        
        m = (y == yy)
        minx = np.min(x[m])
        maxx = np.max(x[m])
        
        xx = np.arange(minx, maxx)
        x = np.append(x, xx)
        y = np.append(y, np.repeat(yy, np.shape(xx)[0]))
        
        
    mask[y, x] = 1





def profile_line(image, src, dst, linewidth=1,
                 order=1, mode='constant', cval=0.0):
    """Return the intensity profile of an image measured along a scan line.

    Parameters
    ----------
    image : numeric array, shape (M, N[, C])
        The image, either grayscale (2D array) or multichannel
        (3D array, where the final axis contains the channel
        information).
    src : 2-tuple of numeric scalar (float or int)
        The start point of the scan line.
    dst : 2-tuple of numeric scalar (float or int)
        The end point of the scan line. The destination point is *included*
        in the profile, in contrast to standard numpy indexing.
    linewidth : int, optional
        Width of the scan, perpendicular to the line
    order : int in {0, 1, 2, 3, 4, 5}, optional
        The order of the spline interpolation to compute image values at
        non-integer coordinates. 0 means nearest-neighbor interpolation.
    mode : {'constant', 'nearest', 'reflect', 'mirror', 'wrap'}, optional
        How to compute any values falling outside of the image.
    cval : float, optional
        If `mode` is 'constant', what constant value to use outside the image.

    Returns
    -------
    return_value : array
        The intensity profile along the scan line. The length of the profile
        is the ceil of the computed length of the scan line.

        Compared to the skimage version, this one returns the sum along the perpendicular direction, not the mean.

    Examples
    --------
    >>> x = np.array([[1, 1, 1, 2, 2, 2]])
    >>> img = np.vstack([np.zeros_like(x), x, x, x, np.zeros_like(x)])
    >>> img
    array([[0, 0, 0, 0, 0, 0],
           [1, 1, 1, 2, 2, 2],
           [1, 1, 1, 2, 2, 2],
           [1, 1, 1, 2, 2, 2],
           [0, 0, 0, 0, 0, 0]])
    >>> profile_line(img, (2, 1), (2, 4))
    array([ 1.,  1.,  2.,  2.])
    >>> profile_line(img, (1, 0), (1, 6), cval=4)
    array([ 1.,  1.,  1.,  2.,  2.,  2.,  4.])

    The destination point is included in the profile, in contrast to
    standard numpy indexing.
    For example:

    >>> profile_line(img, (1, 0), (1, 6))  # The final point is out of bounds
    array([ 1.,  1.,  1.,  2.,  2.,  2.,  0.])
    >>> profile_line(img, (1, 0), (1, 5))  # This accesses the full first row
    array([ 1.,  1.,  1.,  2.,  2.,  2.])
    """
    perp_lines = _line_profile_coordinates(src, dst, linewidth=linewidth)
    if image.ndim == 3:
        pixels = [ndi.map_coordinates(image[..., i], perp_lines,
                                      order=order, mode=mode, cval=cval)
                  for i in range(image.shape[2])]
        pixels = np.transpose(np.asarray(pixels), (1, 2, 0))
    else:
        pixels = ndi.map_coordinates(image, perp_lines,
                                     order=order, mode=mode, cval=cval)
        
#    intensities = pixels.mean(axis=1)
    intensities = pixels.sum(axis=1)

    return intensities

def profile_line_pos(image, src, dst, cval=0.0):
    src_row, src_col = src = np.asarray(src, dtype=float)
    dst_row, dst_col = dst = np.asarray(dst, dtype=float)
    d_row, d_col = dst - src
    theta = np.arctan2(d_row, d_col)

    length = int(np.ceil(np.hypot(d_row, d_col) + 1))
    # we add one above because we include the last point in the profile
    # (in contrast to standard numpy indexing)
    line_col = np.linspace(src_col, dst_col, length)
    line_row = np.linspace(src_row, dst_row, length)


    return np.array()

    

def _line_profile_coordinates(src, dst, linewidth=1):
    """Return the coordinates of the profile of an image along a scan line.

    Parameters
    ----------
    src : 2-tuple of numeric scalar (float or int)
        The start point of the scan line.
    dst : 2-tuple of numeric scalar (float or int)
        The end point of the scan line.
    linewidth : int, optional
        Width of the scan, perpendicular to the line

    Returns
    -------
    coords : array, shape (2, N, C), float
        The coordinates of the profile along the scan line. The length of the
        profile is the ceil of the computed length of the scan line.

    Notes
    -----
    This is a utility method meant to be used internally by skimage functions.
    The destination point is included in the profile, in contrast to
    standard numpy indexing.
    """
    src_row, src_col = src = np.asarray(src, dtype=float)
    dst_row, dst_col = dst = np.asarray(dst, dtype=float)
    d_row, d_col = dst - src
    theta = np.arctan2(d_row, d_col)

    length = int(np.ceil(np.hypot(d_row, d_col) + 1))
    # we add one above because we include the last point in the profile
    # (in contrast to standard numpy indexing)
    line_col = np.linspace(src_col, dst_col, length)
    line_row = np.linspace(src_row, dst_row, length)

    # we subtract 1 from linewidth to change from pixel-counting
    # (make this line 3 pixels wide) to point distances (the
    # distance between pixel centers)
    col_width = (linewidth - 1) * np.sin(-theta) / 2
    row_width = (linewidth - 1) * np.cos(theta) / 2
    perp_rows = np.array([np.linspace(row_i - row_width, row_i + row_width,
                                      linewidth) for row_i in line_row])
    perp_cols = np.array([np.linspace(col_i - col_width, col_i + col_width,
                                      linewidth) for col_i in line_col])
    return np.array([perp_rows, perp_cols])




def pixccd2pixim(pixccd, pixsidesize):
    """Change coordinate system from pixel image to pixel CCD
    

; CATEGORY:
;  image processing
;
; DESCRIPTION:
;  This program is useful when measuring features in images of
;   different resolution. It permits to deal with the shift induced by
;   the rebining of images. The method is to convert image pixel
;   position into physical CCD pixel position.
;
;  pixel image : pixel position in an image, centered on the center 
;                of the pixel. It's given in units of pixels. 
;                The pixel 0,0 is at the center of the
;                bottom left pixel.
;
;  pixel CCD : pixel position on the CCD. It's centered on the bottom
;              left corner of each pixel. It's given in units of 
;              distance. The position 0,0 is at the bottom left
;              corner of the bottom left pixel of the CCD.
;  pixsidesize : size of the pixel size in distance units 
;               (mm for example)
;
; INPUTS:
;  pixccd: pixel position on the CCD
;  pixsidesize : size of the pixel size in distance units (mm for example)
;
; OUTPUTS:
;  return : the pixel position in image coordinate
;
;-    
    
    """

    return pixccd / pixsidesize - 0.5




def pixim2pixccd(pixim,pixsidesize):
    """Change coordinate system from pixel image to pixel CCD
    
; CATEGORY:
;  image processing
;
; DESCRIPTION:
;  This program is useful when measuring features in images of
;   different resolution. It permits to deal with the shift induced by
;   the rebining of images. The method is to convert image pixel
;   position into physical CCD pixel position.
;
;  pixel image : pixel position in an image, centered on the center 
;                of the pixel. It's given in units of pixels. 
;                The pixel 0,0 is at the center of the
;                bottom left pixel.
;
;  pixel CCD : pixel position on the CCD. It's centered on the bottom
;              left corner of each pixel. It's given in units of 
;              distance. The position 0,0 is at the bottom left
;              corner of the bottom left pixel of the CCD.
;  pixsidesize : size of the pixel size in distance units 
;               (mm for example)
;
; INPUTS:
;  pixim: pixel position on the image
;  pixsidesize : size of the pixel size in distance units (mm for example)
;
; OUTPUTS:
;  return : the pixel position in physical CCD coordinate
;
    
    """

    return pixim * pixsidesize + pixsidesize / 2.





def piximchangereso(pixim, reso):
    """Change resolution of a pixel image coordinate
 
#; CATEGORY:
#;  image processing
#;
#; DESCRIPTION:
#;  This program is useful when measuring features in images of
#;   different resolution. It permits to deal with the shift induced by
#;   the rebining of images. The method is to convert image pixel
#;   position into physical CCD pixel position.
#;
#;  pixel image : pixel position in an image, centered on the center 
#;                of the pixel. It's given in units of pixels. 
#;                The pixel 0,0 is at the center of the
#;                bottom left pixel.
#;
#;  pixel CCD : pixel position on the CCD. It's centered on the bottom
#;              left corner of each pixel. It's given in units of 
#;              distance. The position 0,0 is at the bottom left
#;              corner of the bottom left pixel of the CCD.
#;  pixsidesize : size of the pixel size in distance units 
#;               (mm for example)
#;
#; INPUTS:
#;  pixim: pixel position on the image
#;  reso : resolution:  < 0 if want to decrease resolution
#;                      = 0 does nothing
#;                      > 0 increase resolution
#;
#; OUTPUTS:
#;  return : the pixel position in image coordinate in the new resolution    
    
    """    

    return pixccd2pixim(pixim2pixccd(pixim, 2.**reso), 1.)



def loadLocalConf(fnConfig, defaultLocalConf=None):
    """Load the yaml configuration file"""
        
    # -- check if local config file exist
    if (not fnConfig.exists()) and (defaultLocalConf is not None):
        with open(fnConfig, 'w') as file:
            file.write(defaultLocalConf)
            raise Exception(f'Default file {fnConfig} Created. Please edit before rerunning.')
    
    with open(fnConfig) as f:    
        localConf = yaml.load(f, Loader=yaml.FullLoader)
        # To access the yaml variables, use the following syntax
        # localConf['ccor1_datapath']

    return localConf


def progress_bar(current, total, bar_length=20, message=''):
    """Simple but useful text based progress bar"""
    if total == 0:
        total = 1
    percent = float(current) * 100 / total
    progress = int(bar_length * current / total)
    bar = '=' * progress + '-' * (bar_length - progress)
    if len(message) > 0 :
        message = ' : ' + message
    print(f'\r[{bar}] {percent:.1f}%{message}', end='')
    if int(percent) >= 100: print()


if __name__ == '__main__':
    pass