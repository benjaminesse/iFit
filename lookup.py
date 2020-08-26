# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 16:35:14 2020

Script to interpolate a single SO2 spectrometer data point onto light dilution
curves to create a best estimate of SO2 and the light dilution factor

@author: Matthew Varnam - The University of Manchester
@email : matthew.varnam(-at-)manchester.ac.uk
"""

# =============================================================================
# Import additional python libraries
# =============================================================================

# Numpy is a useful mathematics library that is good at handling arrays
import numpy as np

# Matplotlib produces publication quality figures
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Use pandas for the Dataframe class that makes exporting results table easy
import pandas as pd

# Glob allows listing of all files in a particular directory
from glob import glob

#Import shapely to create easily queryable objects
from shapely.geometry import Point, Polygon
from shapely.strtree  import STRtree

# =============================================================================
# Create new functions
# =============================================================================

def create_polygons(px,py):
    '''Create triangles from a 2D data grid'''

    # Record shape of the data
    shape    = np.shape(px)
    
    # Create empty array to store answer
    polygons = np.empty((shape[0]-1,shape[1]-1,2,3,2))
    indices  = np.empty((shape[0]-1,shape[1]-1,2,3),dtype = int)
    
    # Create and store the polygons
    for j in np.arange(shape[1]-1):
        for i in np.arange(shape[0]-1):
            for k in np.arange(2):
                if k == 0:
                    # Create single polygon corners
                    polygon = np.array([[px[i,j]    ,py[i,j]],
                                        [px[i+1,j]  ,py[i+1,j]],
                                        [px[i,j+1]  ,py[i,j+1]]])
                    
                elif k == 1:
                    # Create second polygon corners
                    polygon = np.array([[px[i+1,j]  ,py[i+1,j]],
                                        [px[i+1,j+1],py[i+1,j+1]],
                                        [px[i,j+1]  ,py[i,j+1]]])
                                
                # Store answer in correct location
                indices[i,j,k]  = np.array([[i,j,k]])  
                polygons[i,j,k] = polygon

    return indices,polygons

def check_polygon(point,poly):
    '''Find barycentric coordinates a,b and c of point in triangle'''
    
    # Extract x and y coordinates of point and polygon
    x0,y0 = point
    x1,y1 = poly[0]
    x2,y2 = poly[1]
    x3,y3 = poly[2]
    
    # Find a and b denominator
    denominator = ((y2-y3)*(x1-x3) + (x3-x2)*(y1-y3))
    
    # Find a and b numerators
    a_numerator = ((y2-y3)*(x0-x3) + (x3-x2)*(y0-y3))
    b_numerator = ((y3-y1)*(x0-x3) + (x1-x3)*(y0-y3))
    
    # Find barycentric coordinates
    a = np.divide(a_numerator,denominator)
    b = np.divide(b_numerator,denominator)
    c = 1 - a - b

    return(a,b,c) 

def ellipse(point,error,grid306,grid312):
    '''
    Calculate where points in a grid are relative to an ellipse centred 
    on x,y with a half width of error
    '''
    # Calculate x component of ellipse
    xnum = np.power(np.subtract(grid306,point[0]),2)
    xdenom = np.power(error[0],2)
    xval = np.divide(xnum,xdenom)
    
    # Calculate y component of ellipse
    ynum = np.power(np.subtract(grid312,point[1]),2)
    ydenom = np.power(error[1],2)  
    yval = np.divide(ynum,ydenom)
    
    return np.add(xval,yval)

def calc_uncertainty(point,error,curves0,curves1):
    '''
    Fine the maximum and minimum lookup table SO2 and LDF that lies within
    error of the uncorrected SO2 SCD data point
    '''
    
    # Extract error of the current set of points
    point_err = points_err[i]
    
    # Calculate shape of curves object
    shape = np.shape(curves0)
    
    # Find any points that lie inside the error ellipse
    point_args = np.where(ellipse(point,point_err,curves0,curves1)<1)
    
    # Raise exception if error ellipse lies entirely outside lookup table
    if len(point_args[0]) == 0:
        raise ValueError('No lookup points inside error ellipse')
    
    # Determine the maximum and minimum so2 and ldf indices of points inside 
    # the ellipse
    x_idx_min, x_idx_max = np.min(point_args[0]), np.max(point_args[0])
    y_idx_min, y_idx_max = np.min(point_args[1]), np.max(point_args[1])
    
    # If the ellipse intersects the outer edges of the lookup table then 
    # choose the max and min values modelled
    if x_idx_max != shape[0] - 1:
        x_idx_max += 1
        
    elif x_idx_min != 0:   
        x_idx_min -= 1
         
    if y_idx_max != shape[1] - 1: 
        y_idx_max += 1
          
    elif y_idx_min != 0:
        y_idx_min -= 1
    
    # Group outputs
    x_idx = [x_idx_min,x_idx_max]
    y_idx = [y_idx_min,y_idx_max]
    
    return point_args,x_idx,y_idx

def calc_value(point,tree,polygons,indices,index_by_id):
    '''
    Interpolate within the lookup table to find the best guess of SO2 and LDF.
    Using STR tree here dramatically improves performance as most invalid 
    polygons are removed.
    '''    
    # Find id of possible polygons. Using STR tree here dramatically improves
    # performance as most invalid polygons are removed.
    valid_ids = [(index_by_id[id(poly)]) for poly in tree.query(point)]
    valid_indices = indices[valid_ids]
    
    # Find possible polygons
    polygons_poss = np.array([polygons[i,j,k] for i,j,k in valid_indices])
    
    # Create empty arrays to store correct answer
    index_answer    = []
    bary_cords_list = []  
    
    # Iterate over all polygons that could contain the point as determined by
    # the STR tree
    for i,poly in enumerate(polygons_poss):
        
        # Calculate barycentric coordinates
        bary_cords = check_polygon(np.array(point),poly)
        
        # Check that the barycentric coordinates are all between 
        # 0 and 1
        if bary_cords[0] >= 0 and bary_cords[0] <= 1 and \
            bary_cords[1] >= 0 and bary_cords[1] <= 1 and \
                bary_cords[2] >= 0 and bary_cords[2] <= 1:
                    index_answer.append(valid_indices[i])
                    bary_cords_list.append(bary_cords)                   
                    
    # Check number of answers is 1 (point lies inside a polygon)                 
    if len(index_answer) == 1:
        
        # Create flag to show single answer
        answer_flag = 'Single'
        
    # Check if there are no answers (point lies outside all polygons)    
    elif len(index_answer) == 0:
        
        # Create flag to show no answer
        answer_flag = 'None'
    
    # Raise an error if there is more than one answer (very unlikely due to
    # precision of calculation). This should be handled better, but will only 
    # arise if a point lies exactly on a polygon edge. Edit code to handle
    # in future
    elif len(index_answer) > 1:
        print('Using single solution of ' + str(len(index_answer)))
        
        index_answer = [index_answer[0]]
        bary_cords_list = [bary_cords_list[0]]
        
        # Create flag to show multiple answers
        answer_flag = 'Single'
    
    return index_answer,bary_cords_list,answer_flag

def load_points(fpath):
    '''
    Load in uncorrected SO2 SCDs from a normal ifit analysis
    '''
    # Check file type to determine if in spec format or scan format
    if file_type == 'spec':
        
        # Select file name for points for both wavelength
        points_0_fname = fpath + ('ifit_out_306_316.csv')
        points_1_fname = fpath + ('ifit_out_312_322.csv')
        
        # Load SO2 at first wavelenth
        df0 = pd.read_csv(points_0_fname)
        
        # Load SO2 at second wavelength
        df1 = pd.read_csv(points_1_fname)
                        
    # Extract SO2 and SO2 error info
    points_0 = np.array(df0['SO2'])
    error_0  = np.array(df0['SO2_err'])
    points_1 = np.array(df1['SO2'])
    error_1  = np.array(df1['SO2_err'])

    # Reshape imported data
    points     = np.vstack((points_0,points_1)).T
    points_err = np.vstack((error_0 ,error_1 )).T
        
    return points,points_err

# =============================================================================
# Setup interpolation
# =============================================================================

# Load light dilution curves
mst_dir = 'F:/Scripts/iFit3.4/Dilution_Curves/'
curves306_raw = np.load(mst_dir + 'Telica_sum/306-316.npy')
curves312_raw = np.load(mst_dir + 'Telica_sum/312-322.npy')

# Define so2 and light dilution grid used in light dilution curve generation
so2_grid = np.arange(0,7510,20)
ldf_grid = np.arange(0,1.00,0.002)

plot_flag = True

# Define location of ifit results
point_dir = ('F:/Black_Hole/Data/Burton_Data/SO2Cam Data/Telica Summit 090513/')

# Load points via function
points,points_err = load_points(point_dir)

# =============================================================================
# Pre-processing and loading calculations
# =============================================================================

# Seperate error from curves
curves306,error306 = curves306_raw
curves312,error312 = curves312_raw

# Make polygons out of curves
indices,polygons = create_polygons(curves306,curves312)

# Reshape polygon object to enable STRtrees 
poly_shape = polygons.shape
shaped_polygons = polygons.reshape(poly_shape[0]*poly_shape[1]*2,3,2)
shaped_indices  = indices.reshape(poly_shape[0]*poly_shape[1]*2,3)

# Record dimensions of curve array 
shape = np.shape(curves306)

# Convert test point from molecules/cm2 to ppm.m
points     = np.divide(points     , 2.652E15)
points_err = np.divide(points_err , 2.652E15)

# List column names
col_names = ('Number','SO2_0','SO2_1',
             'SO2','SO2_min','SO2_max',
             'LDF','LDF_min','LDF_max')

n_spec = np.arange(len(points))

# Create dataframe
results_df = pd.DataFrame(index = n_spec, columns = col_names)

# Use shapely objects to improve speed
points_shapely   = [Point(point) for point in points]

# Ignore IllegalArgumentException on this line
poly_shapely = [Polygon(poly) for poly in shaped_polygons]

# Create dictionary to index point list for faster querying
index_by_id = dict((id(poly), i) for i, poly in enumerate(poly_shapely))

# Create STRtree
tree = STRtree(poly_shapely)

for i,point in enumerate(points):
    if i%50 == 0:
        print(i)
    
    # Extract error and shapely point on current loop
    point_shapely = points_shapely[i]
    point_err = points_err[i]
    
    # =========================================================================
    # Calculate uncertainty
    # =========================================================================
    
    try:
        point_args,x_idx,y_idx = calc_uncertainty(point,point_err,
                                                  curves306,curves312)
        
        so2_min,so2_max = so2_grid[x_idx]
        ldf_min,ldf_max = ldf_grid[y_idx]
        
    except:
        so2_min,so2_max = np.nan,np.nan
        ldf_min,ldf_max = np.nan,np.nan

    # =========================================================================
    # Calculate best guess
    # =========================================================================

    answers,bary_coords,answer_flag = calc_value(point_shapely,tree,
                                                 polygons,
                                                 shaped_indices,index_by_id)
    
    # Make sure that only a single answer was found.
    if answer_flag == 'Single':
        
        answer = answers[0]
        bary_coord = bary_coords[0]
        
        # Find vertices of triangle containing point
        vertices_so2 = polygons[(answer[0],answer[1],answer[2])]
        vertices_so2 = np.vstack((vertices_so2,vertices_so2[0]))
        
        # Extract triangle index to find vertices
        j,k = answer[0],answer[1] 
        if answer[2] == 0:
            # Create polygon corners for triangle type a
            vertices_model = np.array([[so2_grid[j]  ,ldf_grid[k]],
                                       [so2_grid[j+1],ldf_grid[k]],
                                       [so2_grid[j]  ,ldf_grid[k+1]]])
        else:       
            # Create polygon corners for triangle type b
            vertices_model = np.array([[so2_grid[j+1],ldf_grid[k]  ],
                                       [so2_grid[j+1],ldf_grid[k+1]],
                                       [so2_grid[j]  ,ldf_grid[k+1]]])
            
        # Create best guess using the barycentric coordinates
        so2_best,ldf_best = np.sum(vertices_model.T*bary_coords[0],
                                   axis = 1)
        
        # Check that for both ldf and so2, min <= best <= max
        if so2_max < so2_best:so2_max = so2_best
        if so2_min > so2_best:so2_min = so2_best
        if ldf_max < ldf_best:ldf_max = ldf_best
        if ldf_min > ldf_best:ldf_min = ldf_best
     
    # If no anaswer is found then the uncorrected SO2 SCDs lie outside graph
    # Examine location of the point to determine best way forward
    elif answer_flag == 'None':
        
        # Check if point is below error zero on either uncorrected SO2 SCD
        if (point[0] < -(np.mean(2*points_err[0])/2.652e15) or
            point[1] < -(np.mean(2*points_err[1])/2.652e15)):
            so2_best = np.nan
            ldf_best = np.nan
        
        else:
            # Give undefined best estimate if error range is entire dataset
            if so2_min == np.nan or (so2_min == 0 and so2_max == so2_grid[-1]):
                so2_best = np.nan
                    
            # Set SO2 to average of the maximum and minimum error value
            else:
                so2_best = np.mean((so2_min,so2_max))
        
            # If 306 is bigger than 312, then set best ldf guess to 0
            if point[0] > point[1]:
                ldf_best = 0.0
               
            # If 306 is smaller than 312 and no answer, LDF must be large or
            # outside modelled data range
            else:
                ldf_best = 1.0
    
    # If data point is within error of zero, give an undefined SO2 
    if np.logical_or((point[0] - 2*point_err[0] < 0), 
                     (point[1] - 2*point_err[1] < 0)):
        so2_best = np.nan
        ldf_best = np.nan
        
    # Combine answers into a single row, which will then go to an output .csv
    row = [i,point[0],point[1],
           so2_best,so2_min,so2_max,
           ldf_best,ldf_min,ldf_max]
    results_df.loc[i] = row 

# Create a plot of the point, error ellipse and light dilution curves
import matplotlib
matplotlib.rcParams.update({'font.size': 14}) 

if plot_flag == True:
    
    # Create figure axes
    fig,ax = plt.subplots(1,2,figsize = [18,8*(9/7)])
    
    # Plot lines of equal model light dilution and so2
    ax[1].plot(curves306[0],curves312[0], zorder = 1,c = 'black',alpha = 0.3,label = 'Lookup grid')
    ax[1].plot(curves306,curves312, zorder = 1,c = 'black',alpha = 0.3)
    ax[1].plot(curves306.T,curves312.T, zorder = 1,c = 'black',alpha = 0.3) 
    
    for ii,x in enumerate(x_idx):
        if ii == 0:
            ax[1].plot(curves306[x],curves312[x],zorder = 2,c = 'dimgrey',
                    linewidth = 5,label = 'SO$_2$ Limit')
        else:
            ax[1].plot(curves306[x],curves312[x],zorder = 2,c = 'dimgrey',
                    linewidth = 5)
    for jj,y in enumerate(y_idx):
        if jj == 0:
            ax[1].plot(curves306.T[y],curves312.T[y],zorder = 2,c = 'darkgray',
                    linewidth = 5,label = 'LDF Limit')
        else:
            ax[1].plot(curves306.T[y],curves312.T[y],zorder = 2,c = 'darkgray',
                    linewidth = 5)

    ax[1].annotate('LDF Max: ' + str(ldf_max), 
                xy=(421, 657), xycoords='data', annotation_clip=False)
    ax[1].annotate('LDF Min: ' + str(ldf_min), 
                xy=(465, 657), xycoords='data', annotation_clip=False)
    ax[1].annotate('SO$_2$ Max: ' + str(so2_max), 
                xy=(440, 657), xycoords='data', annotation_clip=False)
    ax[1].annotate('SO$_2$ Min: ' + str(so2_min), 
                xy=(481.5, 636), xycoords='data', annotation_clip=False) 
    
    letter0 = ax[0].annotate('A)', 
                             xy=(10, 715), xycoords='data') 
    letter0.set_bbox(dict(facecolor='white',edgecolor = 'white'))
    letter1 = ax[1].annotate('B)', 
                             xy=(422, 652), xycoords='data') 
    letter1.set_bbox(dict(facecolor='white',edgecolor = 'white'))
    
    for i in np.arange(len(curves306.T)):
        if i % 50 == 0:
            ax[0].plot(curves306.T[i],curves312.T[i],#alpha = 0.5,
                     label = 'LDF: '+'{0:.1f}'.format(ldf_grid[i]),
                     linewidth = 2)
            ax[1].plot(curves306.T[i],curves312.T[i],#alpha = 0.5,
                     #label = '{0:.1f}'.format(ldf_grid[i]),
                     linewidth = 3)
            
    #plt.legend(title = 'LDF',prop={'size': 11})
    
    # Add datum point to the graph
    #ax[0].scatter(points.T[0][5:750],points.T[1][5:750],zorder = 4,c='C1',
    #           s=15,alpha = 0.5,edgecolor = 'black',label = 'SO$_2$ data')
    ax[0].scatter(point[0],point[1],zorder = 11,c='C0',label = 'Measurement',
               s = 100,marker = 's')
    ax[1].scatter(point[0],point[1],zorder = 11,c='C0',
               s = 100,marker = 's')
    
    # If there is a best guess at the solution, plot triangle it lies inside
    #if answer_flag == 'Single':
    #    ax.plot(vertices_so2.T[0],vertices_so2.T[1],'black',zorder = 2)
      
    # Plot all points that lie within the error ellipse
    ax[1].scatter(curves306[point_args],curves312[point_args],
               c = 'black',s = 15,zorder = 3,label = 'Lookup points \nin error ellipse')

    rect = patches.Rectangle((420,580),60,70,linewidth=1,edgecolor='black',
                             facecolor='none',label = 'Fig 4B) bounds')

    # Add the patch to the Axes
    ax[0].add_patch(rect)

    # Create error ellipse around the datum point
    plot_ellipse0 = patches.Ellipse(point,
                                   point_err[0]*2,point_err[1]*2,fill=False,
                                   zorder = 10,edgecolor='C0',lw = 2,
                                   label = 'Error ellipse')
    plot_ellipse1 = patches.Ellipse(point,
                                   point_err[0]*2,point_err[1]*2,fill=False,
                                   zorder = 10,edgecolor='C0',lw = 2)
    # Add ellipse to graph
    ax[0].add_patch(plot_ellipse0)
    ax[1].add_artist(plot_ellipse1)

    
    ax[0].legend()
    ax[1].legend(loc = 4,framealpha=1)
    
    ax[0].set_xlabel('SO$_2$ at W1 (ppm.m)',fontsize = 16)
    ax[0].set_ylabel('SO$_2$ at W2 (ppm.m)',fontsize = 16)
    ax[0].set_xlim(0,750)
    ax[0].set_ylim(0,750)
    
    ax[1].axis('equal')
    ax[1].set_xlim(420,480)
    ax[1].set_ylim(585,645)
    ax[1].set_xlabel('SO$_2$ at W1 (ppm.m)',fontsize = 16)
    
    fig.subplots_adjust(top=0.925,bottom=0.105,left=0.075,right=0.885,
                        hspace=0.20,wspace=0.10) 

plt.figure()
# Plot lines of equal model light dilution and so2
for i,curve in enumerate(curves306.T):    
    if i % 50 == 0:
        plt.plot(curves306.T[i],curves312.T[i], zorder = 1,c = 'black',alpha = 0.3) 
        
plt.scatter(points.T[0],points.T[1])

results_df.to_csv('lookup_results_bad.csv')