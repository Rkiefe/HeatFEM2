function model = geometry(choice)
if nargin < 1
	choice = 1;
end
model = createpde;

if choice == 1
	% box with Circle
	container_width = 15;
	circle_diameter = 1;

	box_size = [container_width,container_width]; % width and hight
	position = [0,0];	

	rect1 = [3;
		4;
		position(1)-box_size(1)/2;
		position(1)+box_size(1)/2;
		position(1)+box_size(1)/2;
		position(1)-box_size(1)/2;
		position(2)-box_size(2)/2;
		position(2)-box_size(2)/2;
		position(2)+box_size(2)/2;
		position(2)+box_size(2)/2];

	% Magnet
	% inner bottom circle
	position = [0,0];
	radious = circle_diameter;
	circ2 = [1;
			position(1);
			position(2);
			radious];
	circ2 = [circ2;zeros(6,1)];

	% Combine the geometries
	gd = [rect1,circ2];

	% Create names for the shapes
	ns = char('rect1','circ2');
	ns = ns';

	% Set logic for combining shapes
	% sf = 'rect1 + (circ1-rect2) + (circ2-rect2) + (circ3-rect2) + (circ4-rect2) + rect3';
	sf = 'rect1 + circ2'; 

	% Create geometry
	[dl,bt] = decsg(gd,sf,ns);

	geometryFromEdges(model,dl);

else
	% Just two boxes
	
	% Container
	container_width = 15;

	box_size = [container_width,container_width]; % width and hight
	position = [0,0];	

	rect1 = [3;
		4;
		position(1)-box_size(1)/2;
		position(1)+box_size(1)/2;
		position(1)+box_size(1)/2;
		position(1)-box_size(1)/2;
		position(2)-box_size(2)/2;
		position(2)-box_size(2)/2;
		position(2)+box_size(2)/2;
		position(2)+box_size(2)/2];

	% Magnet
	rectangle_dimensions = [2,4]; % width and hight
	position = [0,0];

	rect2 = [3;
		4;
		position(1)-rectangle_dimensions(1)/2;
		position(1)+rectangle_dimensions(1)/2;
		position(1)+rectangle_dimensions(1)/2;
		position(1)-rectangle_dimensions(1)/2;
		position(2)-rectangle_dimensions(2)/2;
		position(2)-rectangle_dimensions(2)/2;
		position(2)+rectangle_dimensions(2)/2;
		position(2)+rectangle_dimensions(2)/2];

	% Combine the geometries
	gd = [rect1,rect2];

	% Create names for the shapes
	ns = char('rect1','rect2');
	ns = ns';

	% Set logic for combining shapes
	% sf = 'rect1 + (circ1-rect2) + (circ2-rect2) + (circ3-rect2) + (circ4-rect2) + rect3';
	sf = 'rect1 + rect2'; 

	% Create geometry
	[dl,bt] = decsg(gd,sf,ns);

	geometryFromEdges(model,dl);

end


end
