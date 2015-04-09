function [S,N,Xout,Yout,Ix] = xy2sn(centerline,X,Y)
% Transforms Cartesian coordinates into channel-fitted coordinates
% 
% SUMMARY
%   This function determines the channel-fitted (S,N) coordinates for any
%   input Cartesian coordinate pairs (X,Y) falling within the supplied
%   curvilinear grid. It works by creating local search polygons along the
%   grid which are resampled with a higher resolution piecewise cubic
%   spline centerline. The minimum distance of any point within a search
%   polygon to the refined centerline is taken as its N coordinate, whereas
%   the arc-length position of the centerline vertex is taken as its S
%   coordinate. Note that while there is discretization error in this
%   approach, the sub-sampling and refinement of the centerline within each
%   search polygon reduces this error somewhat. 
% 
%   This method employs Guneralp and Rhoads (2007) PCS curvature method to
%   create the continuously differentiable centerline. It then uses the
%   forward transformation approach described by Legleiter and Kyriakidis
%   to convert to S,N space.
% 
% INPUTS
%   This function is designed to take the output of kmlcenterline2grid
%   centerline: Structure array containing the following fields
%       s:   Streamwise distance along curvilinear grid path
%       xn:   UTM East coordinates for each grid node
%       yn:   UTM North coordinates for each grid node
%       ds:   Streamwise grid node spacing
%       dn:   Cross-stream grid node spacing
%       beta: Cross-stream grid half width
%   X: vector of Cartesian X coordinates of data to be transformed into S,N
%      coordinate space
%   Y: vector of Cartesian Y coordinates of data to be transformed into S,N
%      coordinate space
% 
% OUTPUTS
%   S: streamwise coordinate of input X,Y Cartesian data 
%   N: cross-stream coordinate of input X,Y Cartesian data
%   X: sorted X coordinates
%   Y: sorted Y coordinates
% 
% EXAMPLE
%   This:
%       [S,N,X,Y] = xy2sn(centerline,Xcart,Ycart);
%   Yields 4 vectors containing the original X,Y coordinate pairs, and
%   their corresponding S,N coordinates within the supplied grid. Points
%   which lie outside of the grid are ommited.
% 
% By: Frank L. Engel, USGS
% Last modified: 20150409

% Create local search polygons
for i = 1:numel(centerline.s)-1
    poly_left_idx(:,:,i) = [...
        i   2;...
        i+1 2;...
        i+1 3;...
        i   3;...
        i   2];
    poly_right_idx(:,:,i) = [...
        i   2;...
        i+1 2;...
        i+1 1;...
        i   1;...
        i   2];
end

% Subsample the centerline between coarse local search polygons
nsubs = 10;
ds = centerline.ds/nsubs;
disp(['Streamwise resolution (ds): ' num2str(ds) ' (units are same as inputs)'])
S = []; N = [];
Xout = []; Yout = [];
Ix = []; % added by KMK
ref_vec = 1:1:(size(X,1)); %added by KMK
% figure(1);clf
% plot(centerline.xn,centerline.yn); hold all
% plot(centerline.xn(:,2),centerline.yn(:,2),'ko')
% axis equal

%%%%%%%%%%%%%%
% RIGHT SIDE %
%%%%%%%%%%%%%%
for j = 1:size(poly_left_idx,3)
    idxl = poly_left_idx(:,:,j);
    %idxr = poly_right_idx(:,:,j);
    %ds  = (max(max(centerline.s(idxl))) - min(min(centerline.s(idxl))))/nsubs;
    
    % Redo the PCS computation on the local polygon
    [pcs_out,x_der1,y_der1] = pcscurvature([centerline.xn(j,2) centerline.xn(j+1,2)],[centerline.yn(j,2) centerline.yn(j+1,2)],ds);
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
        
    % Compute a subseted s corresponding to the working local polygon
    s = centerline.s(j):ds:centerline.s(j+1);
    n = -centerline.beta:centerline.dn:centerline.beta;
    for jj = 1:length(x)
        for ii = 1:numel(n)
            xn(jj,ii)= Ns(jj,1)*n(ii)+x(jj);
            yn(jj,ii) = Ns(jj,2)*n(ii)+y(jj);
        end
    end
    % Discard the 1st point (always will be a scrap)
    xn = xn(2:end,:); yn = yn(2:end,:);
    
    % Construct more precise polygon, and determine points to evaluate
    xv = [centerline.xn(j+1,1);centerline.xn(j+1,2);xn(:,2);centerline.xn(j,2);centerline.xn(j,1)];
    yv = [centerline.yn(j+1,1);centerline.yn(j+1,2);yn(:,2);centerline.yn(j,2);centerline.yn(j,1)];
    [IN,ON] = inpolygon(X,Y,xv,yv);
    Ixx = ref_vec(IN)'; %line added by Kory Konsoer for keeping track of original point order
    %figure(1); hold on
    %patch(xv,yv,'b','EdgeColor','none')
    
    % Find the nearest vertex on the centerline to each X,Y point contained
    % within the current polygon
    px = X(IN|ON)'; py = Y(IN|ON)';
    sx = xn(:,2)'; sy = yn(:,2)';
    [dist,idx,~,~] = closest_points(px,py,sx,sy);
    
    % Handle the case where there are no data points found in the local
    % polygon
    if ~isnan(dist)
        ss = s(idx+1); ss = ss(:);
        S = [S; ss];
        N = [N; dist(:)];
        Xout = [Xout; px(:)];
        Yout = [Yout; py(:)];
        Ix = [Ix; Ixx(:)]; % line added by Kory Konsoer
    else
        S = [S; nan(size(px))];
        N = [N; nan(size(px))];
    end

end


%%%%%%%%%%%%%%
% LEFT SIDE %
%%%%%%%%%%%%%%
for k = 1:size(poly_left_idx,3)
    idxl = poly_left_idx(:,:,j);
    %idxr = poly_right_idx(:,:,j);
    %ds  = (max(max(centerline.s(idxl))) - min(min(centerline.s(idxl))))/nsubs;
    
    % Redo the PCS computation on the local polygon
    [pcs_out,x_der1,y_der1] = pcscurvature([centerline.xn(k,2) centerline.xn(k+1,2)],[centerline.yn(k,2) centerline.yn(k+1,2)],ds);
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
        
    % Compute a subseted s corresponding to the working local polygon
    s = centerline.s(k):ds:centerline.s(k+1);
    n = -centerline.beta:centerline.dn:centerline.beta;
    for jj = 1:length(x)
        for ii = 1:numel(n)
            xn(jj,ii)= Ns(jj,1)*n(ii)+x(jj);
            yn(jj,ii) = Ns(jj,2)*n(ii)+y(jj);
        end
    end
    % Discard the 1st point (always will be a scrap)
    xn = xn(2:end,:); yn = yn(2:end,:);
    
    % Construct more precise polygon, and determine points to evaluate
    xv = [centerline.xn(k+1,3);centerline.xn(k+1,2);xn(:,2);centerline.xn(k,2);centerline.xn(k,3)];
    yv = [centerline.yn(k+1,3);centerline.yn(k+1,2);yn(:,2);centerline.yn(k,2);centerline.yn(k,3)];
    [IN,ON] = inpolygon(X,Y,xv,yv);
    Ixx = ref_vec(IN)'; %line added by Kory Konsoer for keeping track of original point order
    %figure(1); hold on
    %patch(xv,yv,'r','EdgeColor','none')
    
    % Find the nearest vertex on the centerline to each X,Y point contained
    % within the current polygon
    px = X(IN|ON)'; py = Y(IN|ON)';
    sx = xn(:,2)'; sy = yn(:,2)';
    [dist,idx,~,~] = closest_points(px,py,sx,sy);
    
    % Handle the case where there are no data points found in the local
    % polygon
    if ~isnan(dist)
        ss = s(idx+1); ss = ss(:);
        S = [S; ss];
        N = [N; -dist(:)];
        Xout = [Xout; px(:)];
        Yout = [Yout; py(:)];
        Ix = [Ix; Ixx(:)];
    else
        S = [S; nan(size(px))];
        N = [N; nan(size(px))];
    end

end


%%%%%%%%%%%%%%%%%
% SUB FUNCTIONS %
%%%%%%%%%%%%%%%%%

function [out,x_der1,y_der1]=pcscurvature(X,Y,in2)
% PARAMETRIC CUBIC SPLINE INTERPOLATION AND ARC-LENGTH PARAMETERIZATION OF
% DIGITIZED DATA POINTS OF MEANDERING RIVERS - CURVATURE CALCULATION
%
%                           I-B Güneralp
%                           March 24, 2007
%                           Prolific Oven, Palo Alto, CA
%
% in1m;
% ins1m;               % reads the X and Y coodinates of the digitized points
ins1 = [X', Y']; % Smoothed centerline coordinates
in1 = ins1;
rB1 = ins1;
rB1=reorient(ins1);  % moves the first data point to origin (0,0) 
                    % and recalculates the new coords of the rest
                    % of the data points accordingly

% does the piecewise cubic spline (pcs) interpolation:
[ xt, yt, tt ] = ParametricSpline(rB1(:,1),rB1(:,2));
[breaks_xt,coefs_xt,npolys_xt,ncoefs_xt,dim_xt]=unmkpp(xt);
[breaks_yt,coefs_yt,npolys_yt,ncoefs_yt,dim_yt]=unmkpp(yt);
coefs = [coefs_xt; coefs_yt];
npolys = size(coefs,1)/2; ncoefs = 4;
breaks = breaks_xt;
% [breaks,coefs,npolys,ncoefs,dim]=unmkpp(cscvn(rB1'));

% initializations
% in2 = input('Enter the arc-length, delta s: '); % specifies the arc-length (constant)        
spln_number = npolys;       % takes as input from pcs interpolation
splngt = 0;
cntr_t = 1;
ctotlngt = 0;
cintlngt = 0;               % initial value of the current arc-length (containing the remaining part from the previous spline + the increments on the current spline)
delta_t_default = 2;
delta_t = delta_t_default;  % initial cord-increment value to be used in the calculation of arc-lengths
cut_off = 1e-9;             % cut-off value for the arc-length determination
spln_no(cntr_t) = 0;
eq_t(cntr_t) = 0;
rem_intlngt = 0;            % initial value of the length of the remaining part of the previous spline
cintlngt = 0;

% rearrange the coefficient matrix as (ax bx cx dx, ay by cy dy)
coefs_pcs = [coefs(1:npolys,:) coefs(npolys+1:end,:)];
% for i=1:npolys
% 	coefs_pcs(i,:)=[coefs(2*i-1,:) coefs(2*i,:)];
% end;

% rearrange the break points matrix
breaks=breaks';

for counter1 = 1:spln_number    % for each spline returned by cscvnt at line 14
%     fprintf('Spline = %d\n', counter1);
    dax = coefs_pcs(counter1,1); dbx = coefs_pcs(counter1,2); dcx = coefs_pcs(counter1,3); ddx = coefs_pcs(counter1,4); % assigns the pcs' X function coefs
    day = coefs_pcs(counter1,5); dby = coefs_pcs(counter1,6); dcy = coefs_pcs(counter1,7); ddy = coefs_pcs(counter1,8); % assigns the pcs' Y function coefs
    F = @(q)sqrt((3*dax*q.^2+(2*dbx)*q+(dcx)).^2+(3*day*q.^2+(2*dby)*q+(dcy)).^2);
    splnend = breaks(counter1 + 1) - breaks(counter1); % calculates the cord-length of the current spline
    splngt = quad(F,0,splnend); % calculates the spline length for current spline
%     fprintf('Slength %f\n', splngt);
    nextspline = 1;     % resets to 1 
    leftover = splngt;  % initilizes the remaining portion to the whole length of the current spline
    while ((rem_intlngt + leftover) >= in2)     % while the length of the remaining portion of the current spline is equal or longer than the specified arc-length
%         fprintf('Leftover = %d\n', leftover);
        if (nextspline == 1)    % if this is a new spline
            cur_pnt = 0;        % shows the location on the cords (used in the arc-length calculation)
            cintlngt = 0;           % resets the current interval length to zero
            cintlngt = rem_intlngt + quad(F, cur_pnt, cur_pnt + delta_t); % the length of the remaining part of the privious spline + the length of the arc for delta t increment on the current spline
        else
            cur_pnt = eq_t(cntr_t); % assigns the length (location) of the current point to the sum of the cords + the increments on the current cord
            rem_intlngt = 0;        % resets the remaining length to 0 (in the current spline)
            cintlngt = 0;           % resets the current interval length to zero
            cintlngt = quad(F, cur_pnt, cur_pnt + delta_t); % calculates the current interval length for delta t
        end % for -if
        while (abs(cintlngt - in2) > cut_off)   % while the absolute value of the difference between cintlngt and the specified arc-length is larger than the cut-off
            while (cintlngt < in2)              % while cintlngt is smaller than the specified arc-length size
                
                cur_pnt = cur_pnt + delta_t;    % updates current point 
                cintlngt = cintlngt + quad(F, cur_pnt, cur_pnt + delta_t); % updates cintlngt
            end % for -inner while
            if (abs(cintlngt - in2) > cut_off)  % if the absolute value of the difference between cintlngt and the specified arc-length is larger than the cut-off
                cintlngt = cintlngt - quad(F, cur_pnt, cur_pnt + delta_t);  % takes back the most recent addition to cintlngt            
%                 cur_pnt = cur_pnt - delta_t;
                delta_t = delta_t / 10;                                     % updates delta t
                cintlngt = cintlngt + quad(F, cur_pnt, cur_pnt + delta_t);  % calculates the cintlngt using new delta t           
            end % for -if
        end % for -outer while
        
%         fprintf('currentinterval %f\n', cintlngt);
%         fprintf('equalt %f\n',cur_pnt + delta_t);
        cntr_t = cntr_t + 1;   
        spln_no(cntr_t) = counter1;         % keeps track of the spline numbers that correspond to the arc-lengths
        eq_t(cntr_t) = cur_pnt + delta_t;   % keeps track of the total cord-lengths that correspond to the coordinates of the coordinates of the arcs
        x (cntr_t)= dax * eq_t(cntr_t)^3 + dbx * eq_t(cntr_t)^2 + dcx * eq_t(cntr_t) +ddx;  % calculates the x coordinates on the arc-length parameterized composite curve
        y (cntr_t)= day * eq_t(cntr_t)^3 + dby * eq_t(cntr_t)^2 + dcy * eq_t(cntr_t) +ddy;  % calculates the y coordinates on the arc-length parameterized composite curve
        x_der1 (cntr_t)= 3* dax * eq_t(cntr_t)^2 + 2* dbx * eq_t(cntr_t) + dcx;             % calculates the first-order derivative of x
        y_der1 (cntr_t)= 3* day * eq_t(cntr_t)^2 + 2* dby * eq_t(cntr_t) + dcy;             % calculates the first-order derivative of y
        x_der2 (cntr_t)= 6* dax * eq_t(cntr_t) + 2* dbx;                                    % calculates the second-order derivative of x
        y_der2 (cntr_t)= 6* day * eq_t(cntr_t) + 2* dby;                                    % calculates the second-order derivative of y
        curv (cntr_t)= (((x_der1(cntr_t) * y_der2(cntr_t)) - (y_der1(cntr_t) * x_der2(cntr_t))) / (((x_der1(cntr_t))^2) + ((y_der1(cntr_t))^2))^1.5);   % calculates the parametric-curvature of (x,y)
        leftover = leftover - (cintlngt - rem_intlngt);     % the remaining length of the current arc
%         fprintf('Leftover = %f\n', leftover);
%         fprintf('cintlngt = %f\n', cintlngt);
%         fprintf('c eksi rem %f\n', cintlngt - rem_intlngt);
        rem_intlngt = 0;
        nextspline = 0;
        delta_t = delta_t_default;
    end
    rem_intlngt = rem_intlngt + abs(leftover);
%     fprintf('rem_intlngt = %f\n', rem_intlngt);
end % for -counter

    eq_t=eq_t';         % corresponding lengths on the cord-lengths for XY coordinates
    spln_no=spln_no';   % corresponding spline numbers of the XY coordinates
    
    x=x';               % X coordinates of the arc-length parameterized composite curve 
    y=y';               % Y coordinates of the arc-length parameterized composite curve
    x_der1=x_der1';     % 1st order derivatives of the X coordinates of the arc-length parameterized composite curve 
    y_der1=y_der1';     % 1st order derivatives of the Y coordinates of the arc-length parameterized composite curve 
    x_der2=x_der2';     % 2nd order derivatives of the X coordinates of the arc-length parameterized composite curve 
    y_der2=y_der2';     % 2nd order derivatives of the Y coordinates of the arc-length parameterized composite curve 
    curv=curv';         % curvatures of the arc-length parameterized composite curve 

    XY_coord= [x y];    % XY coordinates of arc-length parameterized composite curve
    XY_coord_org=[x+ins1(1,1) y+ins1(1,2)];
    
    for space=1:size(XY_coord,1),
        cikti(space,1) = XY_coord_org(space,1);
        cikti(space,2) = XY_coord_org(space,2);
        cikti(space,3) = curv(space,1);        
    end;
    
    Xplot=XY_coord_org(:,1);
    Yplot=XY_coord_org(:,2);
%     subplot(2,1,1), plot(in1(:,1),in1(:,2),'bo'), hold on, % plots the pcs-interpolated and arc-length parameterized smoothed XY coordinates obtained with arc-length interval given in _conf2  
%     subplot(2,1,1), plot(Xplot,Yplot,'.','MarkerSize',5), hold off; % plots the pcs-interpolated and arc-length parameterized smoothed XY coordinates obtained with arc-length interval given in _conf2  
%     subplot(2,1,2), plot(curv,'.','MarkerSize',5); % plots the curvature series of this planform
    
out=cikti;

% [EOF] pcscurvature

function y=reorient(b)
%   Recalculates the coordinates by moving the first coordinates
%   (i.e. {b(1,1),b(1,2)}) to the origin {0,0}
%
%                           ?nci-Burak Güneralp
%                           March 24, 2007
%                           Prolific Oven, Palo Alto, CA

H=b;                    % reassign the input matrix b to local variable H
dummy(1,1)=H(1,1);      % create dummy matrix populated with
dummy(1,2)=H(1,2);      % first coordinates of H
for r=1:size(H,1),      % for all rows of H
    for c=1:size(H,2),  % for all columns of H
        H(r,c)=H(r,c)-dummy(1,c);   % subtract the first coordinates from all
    end                 % end for-c
end                 % end for-r

y=H;
% [EOF] reorient

function [dist,idx,dist2,idx2] = closest_points(px,py,sx,sy)
% Finds the two closest points in set S to each point in set P
if isempty(px) % no points to find
    dist  = nan;
    idx   = nan;
    dist2 = nan;
    idx2  = nan;
    return
end
for p1 = 1:(size(px,2))
    for p2 = 1:size(sx,2)
        d(p1,p2) = hypot(px(p1)-sx(p2),py(p1)-sy(p2));
    end
    [dist(p1),idx(p1)] = nanmin(d(p1,:)); % Closest break
    d2(p1,:) = d(p1,:); d2(p1,idx(p1))=nan;
    [dist2(p1),idx2(p1)] = nanmin(d2(p1,:)); % Next closest break
end

function [xt,yt,tt] = ParametricSpline(x,y)
%ParametricSpline Generates ppform parametric splines
% Uses matlab's built-in spline function to generate parametric cubic
% splines that are are a f(x,y). This is an alternative to the Curve
% Fitting Toolbox's cscvn function. 
% 
% Taken from: http://www.physicsforums.com/showthread.php?t=468582

arc_length = 0;
n = length(x);
t = zeros(n, 1);
for i=2:n
    arc_length = sqrt((x(i)-x(i-1))^2 + (y(i)-y(i-1))^2);
    t(i) = t(i-1) + arc_length;
end
t=t./t(length(t));
xt = spline(t, x);
yt = spline(t, y);
tt = linspace(0,1,1000);
% [EOF] ParametricSpline