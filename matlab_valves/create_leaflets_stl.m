clear all; close all; clc;

% plotting this function x^2/a^2 + y^2/b^2 = c

% parameters for implicitly defined surface
a = 17;
b = 13;

% discretization parameters for the point cloud
Nz = 100; % number of points in the z direction for the valve
Nz_top = 5; % number of points for the top addition to the valve
Ntheta = 100; 

% geometric parameters for implicitly defined surface
ztop = 0;
zmin = 2;
zmax = 27;
NDIM = 3;

% vectors for making point cloud
zvec = linspace(zmin,zmax,Nz);
%zvec_top = linspace(ztop,zmin - zmin/Nz_top,Nz_top);
zvec_top = linspace(ztop,zmin,Nz_top);
theta = linspace(0,2*pi,Ntheta);

% read in leaflet data at (x,y) tuples
leaflet_data = xlsread('mitral_valve_points.xlsx','coordinates');
leafletXval = leaflet_data(1:end,1);
leafletYval = leaflet_data(1:end,2);
leafletYval(1) = leafletYval(end); % make sure function is periodic

% plot leaflet data
figure(2)
plot(leafletXval, leafletYval,'b');

% find bounding box for leaflet data
Xmin = min(leafletXval);
Xmax = max(leafletXval);
Ymin = min(leafletYval);
Ymax = max(leafletYval);

% scale leaflet data to a shifted unit cube
leafletXval = leafletXval - Xmin; leafletXval = leafletXval./(Xmax - Xmin);
leafletYval = leafletYval - Ymax; leafletYval = leafletYval./(Ymax - Ymin);

% create temporary list of points and tangent vectors for the leaflet data
temp_leafletMatrix = [leafletXval leafletYval];
temp_leafletTangent = zeros(length(leafletXval), 2);
temp_leafletTangent(1:end-1,:) = temp_leafletMatrix(2:end,:) - temp_leafletMatrix(1:end-1,:);
temp_leafletTangent(end,:) = temp_leafletTangent(end-1,:);

% normalize tangent vectors
for ii = 1:length(temp_leafletTangent); temp_leafletTangent(ii,:) = temp_leafletTangent(ii,:)./norm(temp_leafletTangent(ii,:)); end;

% identify problematic points.
% these problematic data points are identified by looking at the dot product between tangent vectors
% which are close together.  if this dot product is negative, then these data points
% are probably low quality and need to be removed.
problem_indicator = zeros(length(temp_leafletTangent)-1, 1);
threshold = 1e-1;
for ii = 1:length(temp_leafletTangent)-1; problem_indicator(ii) = temp_leafletTangent(ii,:)*temp_leafletTangent(ii+1,:)'; end;
problem_indices = find(problem_indicator < threshold);
remove_indices = problem_indices + 1;
temp_leafletTangent(remove_indices,:) = [];
temp_leafletMatrix(remove_indices,:) = [];

% plot scaled leaflet data
figure(3);
for ii = 1:length(leafletXval)
  plot(leafletXval(ii), leafletYval(ii),'o'); hold on
  text(leafletXval(ii),leafletYval(ii),int2str(ii));    
end

% create smooth spline with leaflet data
figure(3)
leafletSpline = cscvn(temp_leafletMatrix');
fnplt(leafletSpline);

% uniformly sample spline, plot it, and store points on spline
% in matrices.  also compute tangent vectors.
figure(11)
Nsample = 1000;
leafletMatrix = zeros(Nsample,2);
leafletTangent = zeros(Nsample,2);
blah = linspace(leafletSpline.breaks(1),leafletSpline.breaks(end),Nsample);
for ii = 1:length(blah)
  spline_val = fnval(leafletSpline,blah(ii));
  plot(spline_val(1), spline_val(2),'x'); hold on
  leafletMatrix(ii,1) = spline_val(1);
  leafletMatrix(ii,2) = spline_val(2);
end
leafletTangent(1:end-1,:) = leafletMatrix(2:end,:) - leafletMatrix(1:end-1,:);
leafletTangent(end,:) = leafletTangent(end-1,:);

% normalize tangent vectors
for ii = 1:length(leafletTangent); leafletTangent(ii,:) = leafletTangent(ii,:)./norm(leafletTangent(ii,:)); end;

% plotting point cloud and populating a matrix
point_cloud = cell(Ntheta, Nz);
for zz = 1:length(zvec)
    c = (1 - (zvec(zz)/93)^(2/3))^3;
    xvec = a*sqrt(c)*cos(theta);
    yvec = b*sqrt(c)*sin(theta);
    for tt = 1:length(theta)
        point_cloud{tt,zz} = [xvec(tt) yvec(tt) zvec(zz)];
    end
end

%creating elements on surface
facets = cell(2*(Nz-1)*(Ntheta-1),1);

count = 1;
for tt = 1:Ntheta-1
   for zz = 1:Nz-1
    facets{count,1} = zeros(NDIM,3);
    facets{count+1,1} = zeros(NDIM,3);
    
    facets{count,1}(1,1) = point_cloud{tt,zz}(1); 
    facets{count,1}(2,1) = point_cloud{tt,zz}(2); 
    facets{count,1}(3,1) = point_cloud{tt,zz}(3); 
    
    facets{count,1}(1,2) = point_cloud{tt,zz+1}(1);
    facets{count,1}(2,2) = point_cloud{tt,zz+1}(2); 
    facets{count,1}(3,2) = point_cloud{tt,zz+1}(3);
    
    facets{count,1}(1,3) = point_cloud{tt+1,zz+1}(1);
    facets{count,1}(2,3) = point_cloud{tt+1,zz+1}(2); 
    facets{count,1}(3,3) = point_cloud{tt+1,zz+1}(3);
    
    facets{count+1,1}(1,1) = point_cloud{tt+1,zz+1}(1); 
    facets{count+1,1}(2,1) = point_cloud{tt+1,zz+1}(2); 
    facets{count+1,1}(3,1) = point_cloud{tt+1,zz+1}(3); 
    
    facets{count+1,1}(1,2) = point_cloud{tt+1,zz}(1); 
    facets{count+1,1}(2,2) = point_cloud{tt+1,zz}(2); 
    facets{count+1,1}(3,2) = point_cloud{tt+1,zz}(3);
    
    facets{count+1,1}(1,3) = point_cloud{tt,zz}(1); 
    facets{count+1,1}(2,3) = point_cloud{tt,zz}(2); 
    facets{count+1,1}(3,3) = point_cloud{tt,zz}(3);
    count = count + 2;
   end
end

%adding facets for top to final_facets array 
for zz = 1:length(zvec_top)-1
    c_zz = (1 - (zvec_top(zz)/93)^(2/3))^3;
    xvec_top_zz = a*sqrt(c_zz)*cos(theta);
    yvec_top_zz = b*sqrt(c_zz)*sin(theta);
    c_zzp1 = (1 - (zvec_top(zz+1)/93)^(2/3))^3;
    xvec_top_zzp1 = a*sqrt(c_zzp1)*cos(theta);
    yvec_top_zzp1 = b*sqrt(c_zzp1)*sin(theta);
    for tt = 1:length(theta)-1
          facets{count,1} = zeros(NDIM,3);
          facets{count+1,1} = zeros(NDIM,3);
    
          facets{count,1}(1,1) = xvec_top_zz(tt); 
          facets{count,1}(2,1) = yvec_top_zz(tt); 
          facets{count,1}(3,1) = zvec_top(zz); 
    
          facets{count,1}(1,2) = xvec_top_zzp1(tt);
          facets{count,1}(2,2) = yvec_top_zzp1(tt); 
          facets{count,1}(3,2) = zvec_top(zz+1);
    
          facets{count,1}(1,3) = xvec_top_zzp1(tt+1);
          facets{count,1}(2,3) = yvec_top_zzp1(tt+1); 
          facets{count,1}(3,3) = zvec_top(zz+1);
    
          facets{count+1,1}(1,1) = xvec_top_zzp1(tt+1); 
          facets{count+1,1}(2,1) = yvec_top_zzp1(tt+1); 
          facets{count+1,1}(3,1) = zvec_top(zz+1); 
    
          facets{count+1,1}(1,2) = xvec_top_zz(tt+1); 
          facets{count+1,1}(2,2) = yvec_top_zz(tt+1); 
          facets{count+1,1}(3,2) = zvec_top(zz);
    
          facets{count+1,1}(1,3) = xvec_top_zz(tt); 
          facets{count+1,1}(2,3) = yvec_top_zz(tt); 
          facets{count+1,1}(3,3) = zvec_top(zz);
          count = count + 2;
    end
end

% determine indicator function on the point cloud
point_cloud_indicator = zeros(Ntheta, Nz);
zvec_shifted = zvec - zmin;
for zz = 1:length(zvec)
    indicator = 0;
    for tt = 1:length(theta)
        z_scaled = -zvec_shifted(zz)./max(abs(zvec_shifted));
        theta_scaled = theta(tt)./max(abs(theta));
        % find the closest (x,y) to (theta,z)
        diff = zeros(length(leafletMatrix),1);
        for ii = 1:length(leafletMatrix)
	  diff(ii) = (leafletMatrix(ii,1) - theta_scaled)^2 + (leafletMatrix(ii,2) - z_scaled)^2; 
	end
        [val ind] = min(diff);
	    normal = [theta_scaled - leafletMatrix(ind,1) z_scaled - leafletMatrix(ind,2)];
	    normal = normal./norm(normal);
	    tangent = [leafletTangent(ind,1) leafletTangent(ind,2)];
        if(normal(1)*tangent(2) - normal(2)*tangent(1) <= 0)
 	      indicator = 1;
        else
	    indicator = 0;
        end
        point_cloud_indicator(tt,zz) = indicator;
    end
end

% plotting indicator function
figure(4);
count_points = 0;
count_squares = 0;
for zz = 1:length(zvec)
   for tt = 1:length(theta)
      if (point_cloud_indicator(tt,zz) == 0)
        plot3(point_cloud{tt,zz}(1), point_cloud{tt,zz}(2), point_cloud{tt,zz}(3), 'r.'); hold on
      else
        plot3(point_cloud{tt,zz}(1), point_cloud{tt,zz}(2), point_cloud{tt,zz}(3), 'b.'); hold on
	    count_points = count_points + 1;
        if(zz <= length(zvec)-1 && tt <= length(theta)-1)
            % count_squares = twice the number of facets in the trimmed stl
            count_squares = count_squares + 1;
        end
      end
   end
end

% create final point cloud with additional top part
% in this loop we just add the points for the top
final_point_cloud = zeros(Nz_top*Ntheta + count_points, 3);
for zz = 1:length(zvec_top)
    c = (1 - (zvec_top(zz)/93)^(2/3))^3;
    xvec_top = a*sqrt(c)*cos(theta);
    yvec_top = b*sqrt(c)*sin(theta);
    for tt = 1:length(theta)
        final_point_cloud(Ntheta*(zz-1) + tt, 1) = xvec_top(tt);
	    final_point_cloud(Ntheta*(zz-1) + tt, 2) = yvec_top(tt);
	    final_point_cloud(Ntheta*(zz-1) + tt, 3) = zvec_top(zz);
    end
end

% add trimmed part of valve to the final point cloud
jj = 1;
for zz = 1:length(zvec)
   for tt = 1:length(theta)
      if (point_cloud_indicator(tt,zz) ~= 0)
          final_point_cloud(Nz_top*Ntheta + jj, 1) = point_cloud{tt,zz}(1);				
          final_point_cloud(Nz_top*Ntheta + jj, 2) = point_cloud{tt,zz}(2);
          final_point_cloud(Nz_top*Ntheta + jj, 3) = point_cloud{tt,zz}(3);
	      jj = jj + 1;
      end
   end
end

% plot final valve
for ii = 1:Nz_top*Ntheta + count_points
    figure(5); plot3(final_point_cloud(ii,1), final_point_cloud(ii,2), final_point_cloud(ii,3), 'k.'); hold on
end


% populate the final facet array
final_facets = cell(2*count_squares + 2*((Nz_top-1)*(Ntheta-1)) , 1);
count_final_facets = 1;
for zz = 1:length(zvec)-1
    for tt = 1:length(theta)-1
      if (point_cloud_indicator(tt,zz) ~= 0)
          final_facets{count_final_facets,1} = zeros(NDIM,3);
          final_facets{count_final_facets+1,1} = zeros(NDIM,3);
    
          final_facets{count_final_facets,1}(1,1) = point_cloud{tt,zz}(1); 
          final_facets{count_final_facets,1}(2,1) = point_cloud{tt,zz}(2); 
          final_facets{count_final_facets,1}(3,1) = point_cloud{tt,zz}(3); 
    
          final_facets{count_final_facets,1}(1,2) = point_cloud{tt,zz+1}(1);
          final_facets{count_final_facets,1}(2,2) = point_cloud{tt,zz+1}(2); 
          final_facets{count_final_facets,1}(3,2) = point_cloud{tt,zz+1}(3);
    
          final_facets{count_final_facets,1}(1,3) = point_cloud{tt+1,zz+1}(1);
          final_facets{count_final_facets,1}(2,3) = point_cloud{tt+1,zz+1}(2); 
          final_facets{count_final_facets,1}(3,3) = point_cloud{tt+1,zz+1}(3);
    
          final_facets{count_final_facets+1,1}(1,1) = point_cloud{tt+1,zz+1}(1); 
          final_facets{count_final_facets+1,1}(2,1) = point_cloud{tt+1,zz+1}(2); 
          final_facets{count_final_facets+1,1}(3,1) = point_cloud{tt+1,zz+1}(3); 
    
          final_facets{count_final_facets+1,1}(1,2) = point_cloud{tt+1,zz}(1); 
          final_facets{count_final_facets+1,1}(2,2) = point_cloud{tt+1,zz}(2); 
          final_facets{count_final_facets+1,1}(3,2) = point_cloud{tt+1,zz}(3);
    
          final_facets{count_final_facets+1,1}(1,3) = point_cloud{tt,zz}(1); 
          final_facets{count_final_facets+1,1}(2,3) = point_cloud{tt,zz}(2); 
          final_facets{count_final_facets+1,1}(3,3) = point_cloud{tt,zz}(3);
          count_final_facets = count_final_facets + 2;
      end   
    end
end

disp('here 1');
count_squares
count_final_facets

%adding facets for top to final_facets array 
for zz = 1:length(zvec_top)-1
    c_zz = (1 - (zvec_top(zz)/93)^(2/3))^3;
    xvec_top_zz = a*sqrt(c_zz)*cos(theta);
    yvec_top_zz = b*sqrt(c_zz)*sin(theta);
    c_zzp1 = (1 - (zvec_top(zz+1)/93)^(2/3))^3;
    xvec_top_zzp1 = a*sqrt(c_zzp1)*cos(theta);
    yvec_top_zzp1 = b*sqrt(c_zzp1)*sin(theta);
    for tt = 1:length(theta)-1
          final_facets{count_final_facets,1} = zeros(NDIM,3);
          final_facets{count_final_facets+1,1} = zeros(NDIM,3);
    
          final_facets{count_final_facets,1}(1,1) = xvec_top_zz(tt); 
          final_facets{count_final_facets,1}(2,1) = yvec_top_zz(tt); 
          final_facets{count_final_facets,1}(3,1) = zvec_top(zz); 
    
          final_facets{count_final_facets,1}(1,2) = xvec_top_zzp1(tt);
          final_facets{count_final_facets,1}(2,2) = yvec_top_zzp1(tt); 
          final_facets{count_final_facets,1}(3,2) = zvec_top(zz+1);
    
          final_facets{count_final_facets,1}(1,3) = xvec_top_zzp1(tt+1);
          final_facets{count_final_facets,1}(2,3) = yvec_top_zzp1(tt+1); 
          final_facets{count_final_facets,1}(3,3) = zvec_top(zz+1);
    
          final_facets{count_final_facets+1,1}(1,1) = xvec_top_zzp1(tt+1); 
          final_facets{count_final_facets+1,1}(2,1) = yvec_top_zzp1(tt+1); 
          final_facets{count_final_facets+1,1}(3,1) = zvec_top(zz+1); 
    
          final_facets{count_final_facets+1,1}(1,2) = xvec_top_zz(tt+1); 
          final_facets{count_final_facets+1,1}(2,2) = yvec_top_zz(tt+1); 
          final_facets{count_final_facets+1,1}(3,2) = zvec_top(zz);
    
          final_facets{count_final_facets+1,1}(1,3) = xvec_top_zz(tt); 
          final_facets{count_final_facets+1,1}(2,3) = yvec_top_zz(tt); 
          final_facets{count_final_facets+1,1}(3,3) = zvec_top(zz);
          count_final_facets = count_final_facets + 2;
    end
end

disp('here 2');
count_final_facets

2*count_squares + 2*((Nz_top-1)*(Ntheta-1))

% outputting implicit surface that is not trimmed
fid = fopen('implicit_surface_test.stl','w');
fprintf(fid,'%s\n', 'solid Volume 1');
for ff = 1:length(facets)
    n = cross(facets{ff,1}(:,3) - facets{ff,1}(:,2), facets{ff,1}(:,1) - facets{ff,1}(:,2));
    n = n./norm(n);
    fprintf(fid,'%s %2.3f %2.3f %2.3f\n', '  facet normal', n(1), n(2), n(3));
    fprintf(fid,'%s\n', '    outer loop');
    fprintf(fid,'%s %2.3f %2.3f %2.3f\n', '      vertex', facets{ff,1}(1,1), facets{ff,1}(2,1), facets{ff,1}(3,1));
    fprintf(fid,'%s %2.3f %2.3f %2.3f\n', '      vertex', facets{ff,1}(1,2), facets{ff,1}(2,2), facets{ff,1}(3,2));
    fprintf(fid,'%s %2.3f %2.3f %2.3f\n', '      vertex', facets{ff,1}(1,3), facets{ff,1}(2,3), facets{ff,1}(3,3));
    fprintf(fid,'%s\n', '    endloop');
    fprintf(fid,'%s\n', '  endfacet');
end
fprintf(fid,'%s\n', 'endsolid Volume 1');
fclose(fid);

% outputting implicit surface that is not trimmed
fid2 = fopen('leaflet_test.stl','w');
fprintf(fid2,'%s\n', 'solid Volume 2');
for ff = 1:length(final_facets)
    n = cross(final_facets{ff,1}(:,3) - final_facets{ff,1}(:,2), final_facets{ff,1}(:,1) - final_facets{ff,1}(:,2));
    n = n./norm(n);
    fprintf(fid2,'%s %2.3f %2.3f %2.3f\n', '  facet normal', n(1), n(2), n(3));
    fprintf(fid2,'%s\n', '    outer loop');
    fprintf(fid2,'%s %2.3f %2.3f %2.3f\n', '      vertex', final_facets{ff,1}(1,1), final_facets{ff,1}(2,1), final_facets{ff,1}(3,1));
    fprintf(fid2,'%s %2.3f %2.3f %2.3f\n', '      vertex', final_facets{ff,1}(1,2), final_facets{ff,1}(2,2), final_facets{ff,1}(3,2));
    fprintf(fid2,'%s %2.3f %2.3f %2.3f\n', '      vertex', final_facets{ff,1}(1,3), final_facets{ff,1}(2,3), final_facets{ff,1}(3,3));
    fprintf(fid2,'%s\n', '    endloop');
    fprintf(fid2,'%s\n', '  endfacet');
end
fprintf(fid2,'%s\n', 'endsolid Volume 2');
fclose(fid2);

