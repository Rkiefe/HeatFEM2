%{
	Solves the heat equation with linear thermal conductivity
	
	The boundary condition: k grad T . n = h(Text-T)
	h is the heat transfer coefficient (reciprocal to the thermal insulance)
%}

close all
clear
clc

% ------------ User defined variables and constants ------------
time = 10;	% time of experiment, s

rho_in = 7.9;			% density, g/cm^3
rho_out = 1.225e-3; 	% g/cm^3 wikipedia

% k_in = 8.8; 		% conductivity, W/m/K (https://www.matweb.com/search/datasheet_print.aspx?matguid=750a3dd8d69b44b79468fbaf72a2beef)
k_in = 4e-2;
k_out = 2.262e-2;	% wikipedia

Cp_in = 0.3; 	% J/g/K
Cp_out = 1.012;	% J/g/K wikipedia

Ti_f = @(x,y) (293);	% Initial temperature function
Text = 280;				% Outside temperature, K
hCoef = 1e-2; 			% Heat transfer coefficient, W/m^2/K

heatSource = 0;	% volumetric heat source (W/cm^3)

% Solver settings
dt = 1e-2; 	% time step, s

% Mesh parameters
hmax = 0; % 0 -> let matlab choose
hmin = 0; % ...
% --------------------------------------------------------------

% PDE model
% model = geometry(2); % 1 for circle, 2 for rectangle
model = createpde();

% Add a rectangle
container_width = 2;
container_height = 1;
position = [0,0];

box_size = [container_width,container_height]; % width and hight
rect = [3;
	4;
	position(1)-box_size(1)/2;
	position(1)+box_size(1)/2;
	position(1)+box_size(1)/2;
	position(1)-box_size(1)/2;
	position(2)-box_size(2)/2;
	position(2)-box_size(2)/2;
	position(2)+box_size(2)/2;
	position(2)+box_size(2)/2];
gd = rect;

% Names for the shapes
ns = char('rect1');
ns = ns';

% Set logic for combining shapes
sf = 'rect1'; 

% Create geometry
[dl,bt] = decsg(gd,sf,ns);
geometryFromEdges(model,dl);

% >> Plot model
pdegplot(model,"EdgeLabels","on","FaceLabels","on");

% >> Mesh
mesh = processMesh(model,Hmax=hmax,Hmin=hmin); % 0 for 'let matlab choose'

% >> Boundary conditions
g = zeros(mesh.numFaces,1);
g(1:4) = hCoef; % [3,4,5,8]

% >> Initial temperature
T = zeros(mesh.nv,1);
for i = 1:mesh.nv
	r = mesh.p(1:2,i);

	T(i) = Ti_f(r(1),r(2));
end

% >> Stiffness matrix
% kTh = zeros(mesh.nt,1) + k_out;		% element wise thermal conductiviy
% kTh(mesh.InsideElements) = k_in;
kTh = zeros(mesh.nt,1) + k_in;
A = stiffnessMatrix(mesh,kTh);

% >> Mass matrix
% c = zeros(mesh.nv,1) + rho_out*Cp_out;
% c(mesh.InsideElements) = rho_in*Cp_in;
c = zeros(mesh.nt,1) + rho_in*Cp_in;
M = massMatrix(mesh,c);

% Load Vector | skiped if heatSource == 0
if heatSource ~= 0
	% q_V = zeros(mesh.nt,1); q_V(mesh.InsideElements) = heatSource;
	q_V = zeros(mesh.nt,1) + heatSource;
	b = loadVector(mesh,q_V);
else
	b = zeros(mesh.nv,1);
end

% Boundary matrix
R = boundaryMatrix(mesh,g);

% Boundary vector
r = boundaryVector(mesh,Text*g);

t = 0;
while t < time
	% Solve matrix equation
	T_new = (M+dt*(A+R))\(dt*(b+r) + M*T);

	% Update
	t = t + dt;
	T = T_new;
end


% >> Plots
figure
% pdegplot(model,"EdgeLabels","on"); hold on
pdemesh(model); hold on; fig = gca;
fig.Children(2).Color = [fig.Children(2).Color,0.1];

% Plot T map
scatter(mesh.p(1,:),mesh.p(2,:),10,T,'filled')
cbar = colorbar;
cbar.Label.String = "T (K)";


% Functions
function mesh = processMesh(model,options)

	arguments
		model;
		options.Hmax = 0;
		options.Hmin = 0;
	end

	generateMesh(model,"GeometricOrder","linear","Hmax",options.Hmax,"Hmin",options.Hmin);

	[~,AE] = area(model.Mesh); % area of each element

	% Process mesh to get a list of edges
	[p,edgeList,t] = meshToPet(model.Mesh);

	% edge list: nd1, nd2, edge index
	edgeList = edgeList([1,2,5],:);

	nv = length(p); % Number of nodes
	nt = length(t); % number of elements
	nedg = length(edgeList); % number of edges

	p = [p;zeros(1,nv)];
	for e = 1:nedg
		p(end,edgeList(1:2,e)) = edgeList(end,e);
	end

	mesh.InsideElements = []; mesh.nInside = 0;
	try % there might not be an inside geometry defined
		mesh.InsideElements = findElements(model.Mesh,"region",Face=2);
		mesh.nInside = numel(mesh.InsideElements);
	end

	mesh = struct;
	mesh.p = p; clear p
	mesh.t = t; clear t
	mesh.nv = nv; clear nv
	mesh.nt = nt; clear nt
	mesh.ne = nedg; clear nedg
	mesh.edgeList = edgeList; clear edgeList
	mesh.AE = AE; clear AE
	mesh.numFaces = model.Geometry.NumEdges;
end

function M = massMatrix(mesh,c)
	M_k = 1/12 * [2,1,1;1,2,1;1,1,2];

	M = zeros(mesh.nv);
	for k = 1:mesh.nt
		nds = mesh.t(1:3,k);

		M(nds,nds) = M(nds,nds) + mesh.AE(k)*M_k*c(k);
	end
end

function A = stiffnessMatrix(mesh,mu)
	A = zeros(mesh.nv,mesh.nv);
	for k = 1:mesh.nt
		
		% Nodes of the element
		nds = mesh.t(1:3,k);

		% For each node
		for i = 1:length(nds)
			[~,bi,ci] = abc(mesh.p,nds,nds(i)); % basis function parameters

			for j = i:length(nds) % matrix is symetric
				[~,bj,cj] = abc(mesh.p,nds,nds(j));

				A(nds(j),nds(i)) = A(nds(j),nds(i)) + mu(k)*(bi*bj + ci*cj)*mesh.AE(k);
				A(nds(i),nds(j)) = A(nds(j),nds(i)); % A is symetric
			end
		end
	end
end

function b = loadVector(mesh,F)
	b = zeros(mesh.nv,1);
	for s = 1:mesh.ne
		nds = mesh.edgeList(1:2,s);
		f = mesh.edgeList(end,s);

		center = mean(mesh.p(1:2,nds),2);
		for i = 1:length(nds)
			nd = nds(i);
			[ai,bi,ci] = abc(mesh.p,nds,nd);

			b(nd) = b(nd) + mesh.AE(k)*F(k)*(ai + bi*center(1) + ci*center(2));
		end
	end
end

function R = boundaryMatrix(mesh,h)
	R = zeros(mesh.nv);

	mat = [2,1;1,2];

	for s = 1:mesh.ne
		nds = mesh.edgeList(1:2,s);
		r = norm(mesh.p(1:2,nds(2))-mesh.p(1:2,nds(1)));
		
		R(nds,nds) = R(nds,nds) + h(mesh.edgeList(end,s))*mat*r/6;
	end
end

function b = boundaryVector(mesh,F)
	b = zeros(mesh.nv,1);
	for s = 1:mesh.ne
		nds = mesh.edgeList(1:2,s);
		r = norm(mesh.p(1:2,nds(2))-mesh.p(1:2,nds(1)));

		b(nds) = b(nds) + 0.5*F(mesh.edgeList(end,s))*r;
	end
end

function [a,b,c] = abc(p,nds,nd) % [a,b,c]
	nd1 = 0;
	nd2 = 0;
	% Find opposing nodes to nd
	for i = 1:length(nds)
		if nds(i) ~= nd && nd1 == 0
			nd1 = nds(i);
		end

		if nd1 ~= nds(i) && nds(i) ~= nd && nd2 == 0
			nd2 = nds(i);
		end
	end
	
	% Vectors in-plane
	v1 = [p(1,nd)-p(1,nd1), p(2,nd)-p(2,nd1),1];
	v2 = [p(1,nd)-p(1,nd2), p(2,nd)-p(2,nd2),1];

	n = cross(v1,v2);

	a = 1 + (n(1)*p(1,nd) + n(2)*p(2,nd))/n(3);
	b = -n(1)/n(3);
	c = -n(2)/n(3);
end