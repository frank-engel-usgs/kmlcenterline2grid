Introduction
==================
 KMLCENTERLINE2GRID creates a curvilinear grid based on an input KML centerline file.
 When deployed, this function behaves differently (and includes a seperate
 help).
 
Summary
--------
   This function fits continuously differentiable piece-wise cubic splines
   to the input centerline, determining the normal vectors at a spacing
   (ds) specified by the user. Then, the lateral spacing (dn) is used to
   contruct points on either side of the PCS centerline up to the
   specified half width (beta).
 
   This function is especially handy for constructing 2-D computational
   and/or sampling grids along curved channels with roughly parallel
   banks.

Installation
------------
1. To run from the source code, you will need Matlab, and the Mapping Toolbox.
2. Also required is the KML Toolbox, found on the [Matlab FileExcange](http://www.mathworks.com/matlabcentral/fileexchange/34694-kml-toolbox-v2-7 "KML Toolbox")
3. Clone this repository using git:
        -if you have a key associated with your github account
        `git clone git@github.com:frank-engel-usgs/kmlcenterline2grid.git`

        -otherwise
        `git clone https://github.com/frank-engel-usgs/kmlcenterline2grid.git`

Usage
-----

###Inputs
   ```matlab
   inputfile: KML file of the grid centerline
   outputfile: name of the output file
   ds: Desired spacing between point in the curvilinear grid in the along
       centerline direction (streamwise)
   dn: Desired spacing between point in the curvilinear grid perpendicular
       to the centerline (cross-stream)
   beta: Half width of the desired grid (half channel width)```

###Outputs
   ```matlab
   centerline: Structure array containing the following fields
       s:   Streamwise distance along curvilinear grid path
       xn:   UTM East coordinates for each grid node
       yn:   UTM North coordinates for each grid node
       ds:   Streamwise grid node spacing
       dn:   Cross-stream grid node spacing
       beta: Cross-stream grid half width
   KML file with the name outputfile containing the computational grid
   CSV file (no header) with the coordinates of the the grid in UTM WGS84```
 
###Examples
   Create a sampling grid for the Missuori River near St. Charles, saved
   as a KML file to be viewed in Google Earth
   The following code:
     ```matlab
     centerline = kmlcenterline2grid('Example.kml','TestOut',100,50,300)```
   Results in:
       ```matlab
       centerline = 
     
            s: [1x62 double]
           xn: [62x13 double]
           yn: [62x13 double]
           ds: 100
           dn: 50
         beta: 300```
 
   **Tip**: To create orthogonal transects, set `dn = beta`

Attribution
----------- 
Original PCS curvature code by: Inci Guneralp (Texas A&M)
Grid Production code by: Frank L. Engel, USGS
 