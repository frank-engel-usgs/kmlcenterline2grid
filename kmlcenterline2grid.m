function varargout = kmlcenterline2grid(varargin)  %inputfile,ds,dn,beta);
% Creates a curvilinear grid based on an input KML centerline file.
% When deployed, this function behaves differently (and includes a seperate
% help).
% 
% SUMMARY
%   This function fits continuously differentiable piece-wise cubic splines
%   to the input centerline, determining the normal vectors at a spacing
%   (ds) specified by the user. Then, the lateral spacing (dn) is used to
%   contruct points on either side of the PCS centerline up to the
%   specified half width (beta).
% 
%   This function is especially handy for constructing 2-D computational
%   and/or sampling grids along curved channels with roughly parallel
%   banks.
%
% INPUTS
%   inputfile: KML file of the grid centerline
%   outputfile: name of the output file
%   ds: Desired spacing between point in the curvilinear grid in the along
%       centerline direction (streamwise)
%   dn: Desired spacing between point in the curvilinear grid perpendicular
%       to the centerline (cross-stream)
%   beta: Half width of the desired grid (half channel width)
%
% OUTPUTS
%   centerline: Structure array containing the following fields
%       s:   Streamwise distance along curvilinear grid path
%       xn:   UTM East coordinates for each grid node
%       yn:   UTM North coordinates for each grid node
%       ds:   Streamwise grid node spacing
%       dn:   Cross-stream grid node spacing
%       beta: Cross-stream grid half width
%   KML file with the name outputfile containing the computational grid
%   CSV file (no header) with the coordinates of the the grid in UTM WGS84
% 
% EXAMPLES
%   Create a sampling grid for the Missuori River near St. Charles, saved
%   as a KML file to be viewed in Google Earth
%   The following code:
%     centerline = kmlcenterline2grid('Example.kml','TestOut',100,50,300)
%   Results in:
%       centerline = 
%     
%            s: [1x62 double]
%           xn: [62x13 double]
%           yn: [62x13 double]
%           ds: 100
%           dn: 50
%         beta: 300
% 
%   Tip: To create orthogonal transects, set dn = beta
%
% Original PCS curvature code by: Inci Guneralp (Texas A&M)
% Grid Production code by: Frank L. Engel, USGS
% Last modified: 8/28/2013

% Parse the inputs
if isdeployed
    if nargin==0 || (nargin==1 && strcmp(varargin{1},'-help'))
    disp({...
        'Creates a curvilinear grid based on an input KML centerline file.';...
        '';...
        'INPUTS';...
        '  inputfile: KML file of the grid centerline';...
        '  outputfile: name of the output file';...
        '  ds: Desired spacing between point in the curvilinear grid in the along';...
        '      centerline direction (streamwise)';...
        '  dn: Desired spacing between point in the curvilinear grid perpendicular';...
        '      to the centerline (cross-stream)';...
        '  beta: Half width of the desired grid (half channel width)';...
        '';...
        'OUTPUTS';...
        '  KML file with the name outputfile containing the computational grid';...
        '  CSV file (no header) with the coordinates of the the grid in UTM WGS84';...
        '';...
        'SUMMARY';...
        '  This function fits continuously differentiable piece-wise cubic splines';...
        '  to the input centerline, determining the normal vectors at a spacing';...
        '  (ds) specified by the user. Then, the lateral spacing (dn) is used to';...
        '  contruct points on either side of the PCS centerline up to the';...
        '  specified half width (beta).';...
        '';...
        '  This function is especially handy for constructing 2-D computational';...
        '  and/or sampling grids along curved channels with roughly parallel';...
        '  banks.';...
        '';...
        'EXAMPLES';...
        '  Create a sampling grid for the Missuori River near St. Charles, saved';...
        '  as a KML file to be viewed in Google Earth';...
        '          text';...
        '    kmlcenterline2grid Example.kml TestOut 100 50 300';...
        '  Tip: To create orthogonal transects, set dn = beta:';...
        '    kmlcenterline2grid Example.kml TestOut 100 300 300)';...
        '';...
        'Original PCS curvature code by: Inci Guneralp (Texas A&M)';...
        'Grid Production code by: Frank L. Engel, USGS'})
    return
    elseif nargin<5 && nargin>1
        error('Too few input arguments. Try the ''-help'' flag')
    elseif nargin>5
        error('Too many input arguments. Try the ''-help'' flag')
    else
        % If deployed, convert required arguments. This is required because inputs
        % passed to the compile function from the command line are always
        % interprested as strings.
        
        inputfile      = varargin{1};
        outputfilename = varargin{2};
        ds             = varargin{3};
        dn             = varargin{4};
        beta           = varargin{5};
        
        ds   = str2double(ds);
        dn   = str2double(dn);
        beta = str2double(beta);
    end
    
else % not deployed
    if nargin==4 % no output file requested
        inputfile      = varargin{1};
        outputfilename = [];
        ds             = varargin{2};
        dn             = varargin{3};
        beta           = varargin{4};
    elseif nargin==5
        inputfile      = varargin{1};
        outputfilename = varargin{2};
        ds             = varargin{3};
        dn             = varargin{4};
        beta           = varargin{5};
    elseif nargin<4
        error('Too few input arguments.')
    elseif nargin>5
        error('Too many input arguments.')
        
    end
end

% Read the input KML file to a MapStruct
indata = kml2struct(inputfile);

% Extract centerline Lat/Lon and convert to UTM
lat = indata.Lat;
lon = indata.Lon;
[x,y,utmzone] = deg2utm(lat,lon);

% Fit PCS and create a regular centerline with spacing ~ds
[pcs_out,x_der1,y_der1]=pcscurvature(x,y,ds);
x=pcs_out(:,1);y=pcs_out(:,2);curv = pcs_out(:,3); % Reassign x,y, and curvature

% Compute theta, tangent, and normal vectors
dx = x_der1; %diff(x);
dy = y_der1; %diff(y);
for i = 2:length(x)
    theta(i) = atan(dx(i-1)/dy(i-1)); % angle of centerline to valley
    dss(i) = sqrt(dx(i-1)^2+dy(i-1)^2); % actual ds (function of tolerance)
    xts(i) = dx(i)/sqrt(dx(i)^2+dy(i)^2);
    yts(i) = dy(i)/sqrt(dx(i)^2+dy(i)^2);
    xns(i) = -dy(i)/sqrt(dx(i)^2+dy(i)^2);
    yns(i) = dx(i)/sqrt(dx(i)^2+dy(i)^2);
end
Ts = [xts', yts']; % tangent vectors
Ns = [xns', yns']; % normal vectors


% Compute and assign (S,N) coords to a grid
s = 0:ds:ds*(numel(x)-1);
n = -beta:dn:beta;
for j = 1:length(x)
    for i = 1:numel(n)
        xn(j,i)= Ns(j,1)*n(i)+x(j);
        yn(j,i) = Ns(j,2)*n(i)+y(j);
    end
end
[S,N] = meshgrid(s,n);

% Screen the first point, which is just a scrap bit of the input centerline
xn = xn(2:end,:);
yn = yn(2:end,:);
s = s(2:end)-ds;

% Convert resulting grid back into DD (lat/lon)
[I,J] = ind2sub(size(xn),1:numel(xn));
[lat_out,lon_out] = utm2deg(xn(:),yn(:),repmat(utmzone(1,:),length(xn(:)),1));

% If required, write the KML and CSV output files
if ~isempty(outputfilename)
    % Write the CSV table
    csvoutmat = [xn(:) yn(:)];
    dlmwrite([outputfilename '.csv'],csvoutmat,'precision', '%12.3f')
    
    % Write a KML, and open in Google Earth
    % The KML will show as a graticule. Also, add markers at the Endpoints
    k = kml(outputfilename);
    % Streamwise lines
    for JJ = 1:numel(xn(1,:))
        k.plot(lon_out(J==JJ),lat_out(J==JJ));
    end
    % Crossstream lines
    for II = 1:numel(xn(:,1))
        k.plot(lon_out(I==II),lat_out(I==II));
    end
    % End points
    lat_eps = [lat_out(J==1); lat_out(J==numel(xn(1,:)))];
    lon_eps = [lon_out(J==1); lon_out(J==numel(xn(1,:)))];
    k.scatter(lon_eps,lat_eps,'iconColor','501400FF','iconScale',0.5);
    k.run
end


% Assign outputs
if nargout>0
    centerline.s = s;
    centerline.xn = xn;
    centerline.yn = yn;
    centerline.ds = ds;
    centerline.dn = dn;
    centerline.beta = beta;
    varargout{1} = centerline;
%     varargout{2} = N;
%     varargout{3} = xn;
%     varargout{4} = yn;
end



function mypostcallback_zoom(obj,evd)
ticks_format('%6.0f','%8.0f'); %formats the ticks for UTM (when zooming)

function mypostcallback_pan(obj,evd)
ticks_format('%6.0f','%8.0f'); %formats the ticks for UTM (when panning)

