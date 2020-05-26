clc
clear all
%% Reading .stl file

file = fopen('D638_TypeI.stl','r');
read_file = fread(file,inf,'uint8=>uint8');
%A binary STL file has an 80-character header(ignored)
num_triangles = typecast(read_file(81:84),'uint32');
%print no. of traingles
num_triangles
%Following the header is a 4-byte integer indicating the number of triangular facets in the file. 
triangles_data = zeros(num_triangles,12);
%Each triangle is described by twelve 32-bit floating-point numbers: 
%three for the normal and then three for the X/Y/Z coordinate of each vertex
%After these follows a 2-byte ("short") unsigned integer that is the "attribute byte count"
% 12*4+2 = 50
tt1 = reshape(read_file(85:end),50,num_triangles);
ttt2 = reshape(typecast(reshape(tt1(1:48,1:num_triangles),1,48*num_triangles),'single'),12,num_triangles)';
triangles_data(:,1:9) = ttt2(:,4:12);
triangles_data(:,10:12) = ttt2(:,1:3);

%% 
    axis='x'; 
if axis == 'x'
    rotated = triangles_data;
elseif axis == 'y'
    rotated = [triangles_data(:,2),triangles_data(:,1),triangles_data(:,3), triangles_data(:,5),triangles_data(:,4),triangles_data(:,6),triangles_data(:,8),triangles_data(:,7),triangles_data(:,9),triangles_data(:,10:12)];
elseif axis == 'z'
    rotated = [triangles_data(:,3:-1:1),triangles_data(:,6:-1:4),triangles_data(:,9:-1:7),triangles_data(:,10:12)];
else
    error('axis needs to be x, y or z');
end
%% 
    axis='z'; 
    theta =0;
 if axis == 'x'
        rotmat = [1, 0, 0; 0, cosd(theta),-sind(theta); 0, sind(theta),cosd(theta)];
    elseif axis == 'y'
        rotmat = [cosd(theta), 0 , sind(theta); 0,1,0; -sind(theta) 0, cosd(theta)];
    elseif axis == 'z'
        rotmat = [cosd(theta), -sind(theta), 0; sind(theta), cosd(theta),0; 0, 0,1];
    else
        error('axis should be x y or z')
    end
    triangles_data(:,1:3) = triangles_data(:,1:3) * rotmat;
    triangles_data(:,4:6) = triangles_data(:,4:6) * rotmat;
    triangles_data(:,7:9) = triangles_data(:,7:9) * rotmat;
    triangles_data(:,1:3:end) = triangles_data(:,1:3:end) - min(min(triangles_data(:,1:3:end)));
    triangles_data(:,2:3:end) = triangles_data(:,2:3:end) - min(min(triangles_data(:,2:3:end)));
    triangles_data(:,3:3:end) = triangles_data(:,3:3:end) - min(min(triangles_data(:,3:3:end)));
    
    %% 
    
slice_height = input('enter the slice height');

min_z = min([triangles_data(:,3); triangles_data(:,6);triangles_data(:,9)])-1e-5;
max_z = max([triangles_data(:,3); triangles_data(:,6);triangles_data(:,9)])+1e-5;

z_slices = min_z: slice_height :max_z;

triangles_data_new = [triangles_data(:,1:12),min(triangles_data(:,[3 6 9]),[],2), max(triangles_data(:,[ 3 6 9]),[],2)];
%gives z min of each  row ( min of z coordinte of each trisngle) in the form of a column vector

%find intersecting triangles

slices = z_slices;

z_triangles = zeros(size(z_slices,2),400);
%triangle with rows equal to no. of slicing planes and columns as a big
%number

z_triangles_size=zeros(size(z_slices,2),1);
%column matrix with no. of rows equal to no. of slicing planes

for i = 1:size(triangles_data_new,1)        
    node = triangles_data_new(i,13);
    % stores minimum of z_coordinate, see line 61.
    high= size(slices,2);
    % its basically no. of slicing planes
    low = 1;
    
    not_match = true;
    a = true;
    b = true;
    
    while not_match
        middle = low + floor((high - low)/2);
        %middle indicates the mid plane where floow gives the greatest
        %integer value of that function
        
        % In the following code slices(middle) represents the z coordinate
        % of the plane and node represents the z_min value of the traingles
        
        if middle == 1 && slices(middle) >= node
            check = 2;
        elseif slices(middle) <= node && middle == size(slices,2)
            check = 2;
            
 %In first condition middle plane goes below the z_mincoordinate  and in
 %the second condition mddle plane goes below the z_max coordinate in
 %these cases we will assign them a value of 2 and later they are forced to
 %leave the llop by changing the a value from true to false
            
        elseif slices(middle)>node && slices(middle-1)<node
            check = 0;
%if the z_min lies in between the two adjacent slicing planes             

        elseif slices(middle)>node
            check = -1;
%If z of plane is greater than z_min of traingle
            
        elseif slices(middle) < node
            check = 1;
%If z of plane is less than z_min of traingle
        end
        
%     check

      if check == -1
          high = middle - 1;
%if this is the case highest node is shifted from high to the plane below
%the middle plane
      elseif check == 1
          low = middle + 1;
%if this is the case lowest node is shifted from low to a plane above the
%middle plane
      elseif check == 0
          node = middle;
          not_match = false;
% in this we change the node to the middle plane and the loop exits with not_match =false,
      elseif high > low || check == 2
          a = false;
          not_match = false;
      end
%if this is the case then we will exit the loop with setting the conditions to false 
    end
    
  
    z_low_index = middle;

 %binary check high
 
    node = triangles_data_new(i,14);
    % stores maximum of z_coordinate, see line 61.
    
    high= size(slices,2);
    low = 1;
    
    not_match = true;
    
    while not_match
        middle = low + floor((high - low)/2);

        if middle == 1 && slices(1) <= node
            check = 2;
        elseif middle == size(slices,2) && slices(middle) <=node
            check = 2;
        elseif slices(middle)>node && slices(middle-1)<node
            check = 0;
        elseif slices(middle)>node
            check = -1;
        elseif slices(middle) < node
            check = 1;
        end

      if check == -1
          high = middle - 1;
      elseif check == 1
          low = middle + 1;
      elseif check == 0
          node = middle;
          not_match = false;
      elseif high > low || check == 2
          b = false;
          not_match = false;
      end

    end
    z_high_index = middle;
    if z_high_index > z_low_index 
        for j = z_low_index:z_high_index-1
            z_triangles_size(j) = z_triangles_size(j) + 1;
            z_triangles(j,z_triangles_size(j)) = i;
        end
    end
end

c=4
%list formed
'list formed'
triangle_checklist2 = z_triangles;
for  k = 1:size(z_slices,2)
    
  triangle_checklist = triangle_checklist2(k,1:z_triangles_size(k));

    [lines,linesize] = triangle_plane_intersection(triangles_data_new(triangle_checklist,:), z_slices(k));

     if linesize ~= 0
            %find all the points assign nodes and remove duplicates
            start_nodes = lines(1:linesize,1:2);
            end_nodes = lines(1:linesize,4:5);
            nodes = [start_nodes; end_nodes];
%             hold on
%             for i=1:size(start_nodes,1)
%             plot3([start_nodes(i,1) end_nodes(i,1)], [start_nodes(i,2) end_nodes(i,2)], [z_slices(k) z_slices(k)])
%             end
            a=1;
            min_x = min(nodes(:,1));
            x_start = min_x+.001;
            max_x = max(nodes(:,1));
            max_y = max(nodes(:,2))+5;
            min_y = min(nodes(:,2))+5;
            while x_start < max_x+0.001
            j=1;
            point_all=zeros(2,2);
            for i = 1:size(start_nodes,1)
            lines = [x_start min_y; x_start max_y; start_nodes(i,1) start_nodes(i,2);end_nodes(i,1) end_nodes(i,2)];
            x = lines(:,1);
            y = lines(:,2);
            % Calculation
            denominator = (x(1)-x(2))*(y(3)-y(4))-(y(1)-y(2))*(x(3)-x(4));
            point = [((x(1)*y(2)-y(1)*x(2))*(x(3)-x(4))-(x(1)-x(2))*(x(3)*y(4)-y(3)*x(4)))/denominator 
            ,((x(1)*y(2)-y(1)*x(2))*(y(3)-y(4))-(y(1)-y(2))*(x(3)*y(4)-y(3)*x(4)))/denominator];
            if point(1)>min(start_nodes(i,1),end_nodes(i,1))&& point(1)<max(start_nodes(i,1),end_nodes(i,1))
            point_all(j,1) = point(1);
            point_all(j,2) = point(2);
            j= j+1;
            end   
            end
            point_all;
%             while b<3
%                fprintf('N%d G00 X%dY%dZ%d \n',c,round(point_all(1,1)),round(point_all(1,2)),round(z_slices(k)))
%                b=5;
%                c=c+1;
%             end
            fprintf('N%d G00 X%dY%dZ%d \n',c,round(point_all(1,1)),round(point_all(1,2)),round(z_slices(k)))
            c=c+1;
            fprintf('N%d G01 X%dY%dZ%d \n',c,round(point_all(2,1),2),round(point_all(2,2),2),round(z_slices(k),2))
            c=c+1;
            x_start = x_start+a;
            end
            
    end
% hold off
end