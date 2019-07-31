% All data is imported from the excel file - "test.xlsx". It needs to be
% located in the same folder as this matlab script. The X and Y coordinates
% are vectors and imported from the sheets Xdata and Ydata respectively.
% The corresponding Z-data which is a matrix is imported from the sheet,
% Zdata. Note that X and Y data vectors need to have the same size if you
% want to puncture holes(create a mesh) in the surface. The data that is
% being imported describes the singlet surface of CO+O. 

% The ranges are as follows:
% X: -2.5 to 2.5 angstroms, size - 200
% Y: 0 to 3 angstroms, size - 200
% Z: f(X,Y), size - 200x200


% surf_print code written by Phalgun Lolur, Missouri S&T.
% Acknowledgements: Sven Holcombe and Paul Kassebaum. 



function [] = surf_print ()
tic
clc;
clear all;
% Import the data from the excel sheet
% Do not add additional information/comments to the excel sheet. It may
% result in errors in data importation.

% 1cm = 20 grid
% 1m = 2000 grid
% X=linspace(-0.25, 0.25, 1000);  % 0.5m, 19.68503937 inch
% Y=linspace(-0.25, 0.25, 1000);  % 0.5m, 19.68503937 inch
% [xgrid,ygrid]=meshgrid(X,Y);
% Z=0.10725*((2*xgrid/0.18).^4+(2*ygrid/0.2).^2);
% Z(Z(:,:)>0.1177) = 0.1177;  % 0.1177 m, 4.6338582677 inch 

X=linspace(0, 0.198, 1000);  % 0.5m, 19.68503937 inch
Y=linspace(0, 0.1782, 1000);  % 0.5m, 19.68503937 inch
[xgrid,ygrid]=meshgrid(X,Y);
Z=ones(1000,1000);

% X=linspace(-0.155, 0.155, 621);  % 0.31m, 12.204724409 inch
% Y=linspace(-0.155, 0.155, 621);  % 0.31m, 12.204724409 inch
% [xgrid,ygrid]=meshgrid(X,Y);
% Z=0.195*((2*xgrid/0.3).^4+(2*ygrid/0.3).^2);
% Z(Z(:,:)>0.192) = 1;  % 1 m, 39.37007874 inch 
% Z(Z(:,:)<=0.192) = 2;  % 2 m, 78.74015748 inch 
% Z(Z(:,:))

if (numel(X)*numel(Y)~= numel(Z))
    fprintf('There is data missing in the excel sheet. \nCorrect the error and run the script again.\n\n');
    return;
end

% Plot a surface from the imported data
surf(X,Y,Z);


% Get the desired dimensions of the final object from the user)
fprintf('STL files dont have units.\nWe set the units in the printer.\nFor convenience let us work in inches.\n');
fprintf('The CO+O singlet example PES in this sheet is: 5" x 3" x 4"\nIt is important to assign consistent dimensions to the object otherwise they result in a stretched/distorted final object.\n\n');
ol=input('Enter the desired length of the object (X dimension, 5 for CO+O example)\n');
ow=input('Enter the desired width of the object (Y dimension, 3 for CO+O example)\n');
oh=input('Enter the desired height of the object (Z dimension, 4 for CO+O example)\n');
xf=max(X(:))-min(X(:));
yf=max(Y(:))-min(Y(:));
zf=max(Z(:))-min(Z(:));

% Scaling the dimensions of the model so that they can be printed as required by the user.

X=X*ol/xf;
Y=Y*ow/yf;
Z=Z*oh/zf;

fprintf('\nThe choice is made to print as a smooth surface, or as a mesh by \nintroducing an array of rectanglar holes.\n');
fprintf('\nNote, using a mesh can save material, but can take longer to generate the STL file \nand may cause problems for the printer if for example the mesh is too fine.\n');
choice=input('\nDo you want to puncture holes (create a mesh) in the surface?\nEnter 1 for yes, 2 for no.\n\n');
fprintf('\nThe outputted STL file is created in the same folder with the name "thick_surf.stl"\n\n');

if (choice==1)
    if (numel(X)~=numel(Y))
        fprintf('X and Y data vectors are not of the same size. \nCorrect the error and run the script again.\n\n');
     return;
    end
        
    n= numel(X);
    h = surf(X,Y,Z);
    set(h,'visible','off');
    p = surf2patch(h);
    sub = [1:n-1,repmat(n,1,n-1),n:-1:2,ones(1,n-1)];
    boundVert = sub2ind(size(Z),...
    sub,...
    circshift(sub,[0 n-1]));
    t=input('Enter the desired thickness of the surface in inches (typically between 0.04 and 0.125)\n');
%     Note that stl files don't deal with units. So, if we specify the
%     thickness as 0.1 units, its considers the thickness in the units we
%     are working in. If we set the units as inches in our printer, the
%     thickness of the surface will be 0.1 inches.

    holesize=input('Enter the desired hole size in each grid from a scale of 0 to 1\n For example, if holesize=0.5, it means that 50 % of the surface will be removed from each grid');

    
    if (holesize >=1 || holesize <0)
%         if holesize =1 means that 100% of the surface will be removed
%         if holesize =0 means that none of the surface is removed
%         if holesize <0 implies a negative holesize
        fprintf('Invalid choice\n');
        return;
    end
%     holesize refers to the fraction of surface that will be
%     removed/punctured from each grid. If we specify holesize as 0.5, half
%     of each grid will be removed.
    boundVert = [boundVert,boundVert(1)];
    s = punctureSurface(p, boundVert, holesize, -t);
    stlwrite('thick_surf.stl',s);
end
if (choice==2)
    t=input('Enter the desired thickness of the surface in inches (typically between 0.04 and 0.125)\n');
    s=surf2solid(X,Y,Z,'thickness',-t);
    stlwrite('thick_surf.stl',s);

end

o=input('Do you want to open the generated STL file? \nEnter 1 for yes, 2 for no.\n');
if (o==1)
    winopen thick_surf.stl
end
if (o~=1)
    return;
end

toc
end


function result = punctureSurface(manifold, boundVert, holeSize, shellThickness)
% PUNCTURESURFACE. Puncture a surface defined by faces and vertices.
% 
% INPUT.
%  manifold - structure with fields describing a surface.
%    manifold.faces
%    manifold.vertices
%
%  boundVert - vector of indices into manifold.vertices that
%              specify the boundary vertices. The first and last components
%              of the vector must be equal, representing a closed circuit.
%
%  holeSize - scalar in the range [0,1] determining the size of holes.
%
%  shellThickness - thickness of the shell created from the original
%                   manifold.
%
%
% OUTPUT.
%  result - structure with fields describing the punctured surface. 
%           Every face is a triangle to facilitate its STL-file 
%           description for computer aided manufacturing.
%    result.faces
%    result.vertices
%
%
% EXAMPLE USE CASE.
% n = 15;
% x = linspace(-1,1,n);
% thickness = 0.05*max(x);
% holeSize = 0.5;
% [x,y] = meshgrid(x);
% z = abs(cos(complex(x,y)));
% figure;
% h = surf(x,y,z);
% set(h,'visible','off');
% p = surf2patch(h);
% sub = [1:n-1,repmat(n,1,n-1),n:-1:2,ones(1,n-1)];
% boundVert = sub2ind(size(z),...
%   sub,...
%   circshift(sub,[0 n-1]));
% boundVert = [boundVert,boundVert(1)];
% s = punctureSurface(p, boundVert, holeSize, thickness);
% zeta = complex(s.vertices(:,1),s.vertices(:,2));
% w = sin(zeta);
% u = real(w);
% v = imag(w);
% s.vertices = [u,v,s.vertices(:,3)];
% ps = patch(s);
% axis vis3d;
% axis equal;
% set(gca,'visible','off');
% set(ps, ...
%     'EdgeColor','none', ...
%     'FaceColor',[0.9 0.2 0.2], ...
%     'FaceLighting','phong', ...
%     'AmbientStrength',0.3, ...
%     'DiffuseStrength',0.6, ... 
%     'Clipping','off',...
%     'BackFaceLighting','lit', ...
%     'SpecularStrength',1, ...
%     'SpecularColorReflectance',1, ...
%     'SpecularExponent',7);
% l1 = light('Position',[40 100 20], ...
%     'Style','local', ...
%     'Color',[0 0.8 0.8]);
% l2 = light('Position',[.5 -1 .4], ...
%     'Color',[0.8 0.8 0]);

% Copyright 2013 The MathWorks, Inc.


% Initialization
holeSize = holeSize*(-0.5) + 0.5; % rescale for convenient input.
puncturedManifold.faces = [];
puncturedManifold.vertices = [];
hole.vertices = [];
hole.faces = {};
verticesAddedSoFar = 0;

% Determine how many edges are on each face
uniformNumberOfEdges = ~iscell(manifold.faces);
if ~uniformNumberOfEdges
  numEdges = zeros(size(manifold.faces,1),1);
  for j = 1:size(numEdges,1)
    numEdges(j) = length(manifold.faces{j});
  end
else
  numEdges = repmat(size(manifold.faces,2),size(manifold.faces,1),1);
  manifold.faces = num2cell(manifold.faces,2);
end

uniqueNumEdges = unique(numEdges);
% for each family of ngons with uniqueNumEdges(i) edges.
% A family of ngons all have the same number of edges.
for i = 1:length(uniqueNumEdges) 
  cellBoundary.faces = reshape(...
    [manifold.faces{numEdges == uniqueNumEdges(i)}],...
    uniqueNumEdges(i), [])';
  
  cellBoundaryFacet = manifold.vertices';
  cellBoundaryFacet = reshape(cellBoundaryFacet(:,cellBoundary.faces'),...
    3, uniqueNumEdges(i), []);
  
  % Create vertices that define each hole for faces
  if mod(uniqueNumEdges(i),2) % ngons with odd number of edges
    % Hole vertex lies on line connecting a cellBoundary vertex and the
    % center of its opposite edge.
    intermediateHoleVertices = cellBoundaryFacet ...
    +(...
      (...
         circshift(cellBoundaryFacet,[0 floor(uniqueNumEdges(i)/2) 0]) ...
       - circshift(cellBoundaryFacet,[0 ceil( uniqueNumEdges(i)/2) 0]) ...
      )./2 ...
      + circshift(cellBoundaryFacet,[0 ceil( uniqueNumEdges(i)/2) 0]) ...
      - cellBoundaryFacet ...
     ).*holeSize;
  else % ngons with even number of edges
    % Hole vertex lies on line connecting a cellBoundary vertex and its
    % opposite cellBoundary vertex.
    intermediateHoleVertices  = cellBoundaryFacet ...
      +(circshift(cellBoundaryFacet,[0 uniqueNumEdges(i)/2 0])...
      -cellBoundaryFacet).*holeSize;
  end
  
  intermediateHoleVertices = reshape(intermediateHoleVertices,3,[])';
  % Collect ngon hole vertices for all n
  hole.vertices = [hole.vertices; intermediateHoleVertices];

  
  intermediateHoleFaces = reshape(1:size(intermediateHoleVertices,1),...
    uniqueNumEdges(i),[])';
  % Account for all hole vertices defined so far from previous ngon
  % families and the original manifold.
  intermediateHoleFaces = intermediateHoleFaces ...
    + verticesAddedSoFar + size(manifold.vertices,1);
  
  % Create triangular faces that connect the cellBoundary vertices to the
  % hole vertices.
  for j = 1:uniqueNumEdges(i)-1
    puncturedManifold.faces = [puncturedManifold.faces;
      cellBoundary.faces(:,j),...
      cellBoundary.faces(:,j+1),...
      intermediateHoleFaces(:,j);...
      ...
      intermediateHoleFaces(:,j),...
      cellBoundary.faces(:,j+1),...
      intermediateHoleFaces(:,j+1)...
      ];
  end
  puncturedManifold.faces = [puncturedManifold.faces;
    cellBoundary.faces(:,uniqueNumEdges(i)), ...
    cellBoundary.faces(:,1), ...
    intermediateHoleFaces(:,uniqueNumEdges(i));...
    ...
    intermediateHoleFaces(:,uniqueNumEdges(i)),...
    cellBoundary.faces(:,1),...
    intermediateHoleFaces(:,1)...
    ];
  
  % Collect ngon faces for all families
  hole.faces = [hole.faces; num2cell(intermediateHoleFaces,2)];
  
  % Keep track of how many vertices were created for this ngon family
  verticesAddedSoFar = verticesAddedSoFar + size(intermediateHoleVertices,1);
end

% Collect all of the vertices defined for each family of ngon.
puncturedManifold.vertices = [...
    manifold.vertices;...
    hole.vertices];

% Bottom layer to add thickness
facets = puncturedManifold.vertices';
facets = reshape(facets(:,puncturedManifold.faces'),...
  3, 3, []);
V1 = squeeze(facets(:,2,:) - facets(:,1,:));
V2 = squeeze(facets(:,3,:) - facets(:,1,:));
normals = V1([2 3 1],:) .* V2([3 1 2],:) - V2([2 3 1],:) .* V1([3 1 2],:);
normals = bsxfun(@times, normals, 1 ./ sqrt(sum(normals .* normals, 1)));
% Compute mean normals of faces shared by each vertex
meanNormal = zeros(3,length(puncturedManifold.vertices));
for k = 1:length(puncturedManifold.vertices)
  % Find all faces shared by a vertex
  [sharedFaces,~] = find(puncturedManifold.faces == k);
  % Compute the mean normal of all faces shared by a vertex
  meanNormal(:,k) = mean(normals(:,sharedFaces),2);
end
meanNormal = bsxfun(@times, meanNormal,...
  shellThickness ./ sqrt(sum(meanNormal .* meanNormal, 1)));
bottom.vertices = puncturedManifold.vertices - meanNormal';
bottom.faces = puncturedManifold.faces;

% Walls of holes
numEdges = zeros(size(hole.faces,1),1);
for j = 1:size(numEdges,1)
  numEdges(j) = length(hole.faces{j});
end
uniqueNumEdges = unique(numEdges);
% for each family of ngons with uniqueNumEdges(i) edges.
% A family of ngons all have the same number of edges.
wall.faces = [];
for i = 1:length(uniqueNumEdges) 
  cellBoundary.faces = reshape(...
    [hole.faces{numEdges == uniqueNumEdges(i)}],...
    uniqueNumEdges(i), [])';
  numVerticesPerLayer = length(puncturedManifold.vertices);
  for j = 1:uniqueNumEdges(i)-1
    wall.faces = [wall.faces;...
      cellBoundary.faces(:,j),...
      cellBoundary.faces(:,j+1),...
      cellBoundary.faces(:,j)+numVerticesPerLayer;...
      ...
      cellBoundary.faces(:,j)+numVerticesPerLayer,...
      cellBoundary.faces(:,j+1),...
      cellBoundary.faces(:,j+1)+numVerticesPerLayer...
      ];
  end
  wall.faces = [wall.faces;...
    cellBoundary.faces(:,uniqueNumEdges(i)), ...
    cellBoundary.faces(:,1), ...
    cellBoundary.faces(:,uniqueNumEdges(i))+numVerticesPerLayer;...
    ...
    cellBoundary.faces( :,uniqueNumEdges(i))+numVerticesPerLayer,...
    cellBoundary.faces(:,1),...
    cellBoundary.faces(:,1)+numVerticesPerLayer...
    ];
end % walls of holes

boundary.vertices = [...
  manifold.vertices(boundVert,:);
  bottom.vertices(boundVert,:)];
% Number of wall vertices on each surface (nwv).
nwv = length(boundary.vertices)/2;
% Allocate memory for wallFaces.
boundary.faces = zeros(2*(nwv-1),3);
% Define the faces.
for k = 1:nwv-1
    boundary.faces(k      ,:) = [k+1  ,k      ,k+nwv];
    boundary.faces(k+nwv-1,:) = [k+nwv,k+1+nwv,k+1];
end

% All together now!
result.vertices = [...
  puncturedManifold.vertices; 
  bottom.vertices; 
  boundary.vertices];
result.faces = [...
  fliplr(puncturedManifold.faces); 
  bottom.faces+numVerticesPerLayer; 
  wall.faces; 
  boundary.faces + 2*numVerticesPerLayer];

end


function stlwrite(filename, varargin)
%STLWRITE   Write STL file from patch or surface data.
%
%   STLWRITE(FILE, FV) writes a stereolithography (STL) file to FILE for a
%   triangulated patch defined by FV (a structure with fields 'vertices'
%   and 'faces').
%
%   STLWRITE(FILE, FACES, VERTICES) takes faces and vertices separately,
%   rather than in an FV struct
%
%   STLWRITE(FILE, X, Y, Z) creates an STL file from surface data in X, Y,
%   and Z. STLWRITE triangulates this gridded data into a triangulated
%   surface using triangulation options specified below. X, Y and Z can be
%   two-dimensional arrays with the same size. If X and Y are vectors with
%   length equal to SIZE(Z,2) and SIZE(Z,1), respectively, they are passed
%   through MESHGRID to create gridded data. If X or Y are scalar values,
%   they are used to specify the X and Y spacing between grid points.
%
%   STLWRITE(...,'PropertyName',VALUE,'PropertyName',VALUE,...) writes an
%   STL file using the following property values:
%
%   MODE          - File is written using 'binary' (default) or 'ascii'.
%
%   TITLE         - Header text (max 80 chars) written to the STL file.
%
%   TRIANGULATION - When used with gridded data, TRIANGULATION is either:
%                       'delaunay'  - (default) Delaunay triangulation of X, Y
%                       'f'         - Forward slash division of grid quads
%                       'b'         - Back slash division of quadrilaterals
%                       'x'         - Cross division of quadrilaterals
%                   Note that 'f', 'b', or 't' triangulations now use an
%                   inbuilt version of FEX entry 28327, "mesh2tri".
%
%   FACECOLOR     - Single colour (1-by-3) or one-colour-per-face (N-by-3) 
%                   vector of RGB colours, for face/vertex input. RGB range
%                   is 5 bits (0:31), stored in VisCAM/SolidView format
%                   (http://en.wikipedia.org/wiki/STL_(file_format)#Color_in_binary_STL)
%
%   Example 1:
%       % Write binary STL from face/vertex data
%       tmpvol = zeros(20,20,20);       % Empty voxel volume
%       tmpvol(8:12,8:12,5:15) = 1;     % Turn some voxels on
%       fv = isosurface(tmpvol, 0.99);  % Create the patch object
%       stlwrite('test.stl',fv)         % Save to binary .stl
%
%   Example 2:
%       % Write ascii STL from gridded data
%       [X,Y] = deal(1:40);             % Create grid reference
%       Z = peaks(40);                  % Create grid height
%       stlwrite('test.stl',X,Y,Z,'mode','ascii')

%   Original idea adapted from surf2stl by Bill McDonald. Huge speed
%   improvements implemented by Oliver Woodford. Non-Delaunay triangulation
%   of quadrilateral surface courtesy of Kevin Moerman. FaceColor
%   implementation by Grant Lohsen.
%
%   Author: Sven Holcombe, 11-24-11


% Check valid filename path
path = fileparts(filename);
if ~isempty(path) && ~exist(path,'dir')
    error('Directory "%s" does not exist.',path);
end

% Get faces, vertices, and user-defined options for writing
[faces, vertices, options] = parseInputs(varargin{:});
asciiMode = strcmp( options.mode ,'ascii');

% Create the facets
facets = single(vertices');
facets = reshape(facets(:,faces'), 3, 3, []);

% Compute their normals
V1 = squeeze(facets(:,2,:) - facets(:,1,:));
V2 = squeeze(facets(:,3,:) - facets(:,1,:));
normals = V1([2 3 1],:) .* V2([3 1 2],:) - V2([2 3 1],:) .* V1([3 1 2],:);
clear V1 V2
normals = bsxfun(@times, normals, 1 ./ sqrt(sum(normals .* normals, 1)));
facets = cat(2, reshape(normals, 3, 1, []), facets);
clear normals

% Open the file for writing
permissions = {'w','wb+'};
fid = fopen(filename, permissions{asciiMode+1});
if (fid == -1)
    error('stlwrite:cannotWriteFile', 'Unable to write to %s', filename);
end

% Write the file contents
if asciiMode
    % Write HEADER
    fprintf(fid,'solid %s\r\n',options.title);
    % Write DATA
    fprintf(fid,[...
        'facet normal %.7E %.7E %.7E\r\n' ...
        'outer loop\r\n' ...
        'vertex %.7E %.7E %.7E\r\n' ...
        'vertex %.7E %.7E %.7E\r\n' ...
        'vertex %.7E %.7E %.7E\r\n' ...
        'endloop\r\n' ...
        'endfacet\r\n'], facets);
    % Write FOOTER
    fprintf(fid,'endsolid %s\r\n',options.title);
    
else % BINARY
    % Write HEADER
    fprintf(fid, '%-80s', options.title);             % Title
    fwrite(fid, size(facets, 3), 'uint32');           % Number of facets
    % Write DATA
    % Add one uint16(0) to the end of each facet using a typecasting trick
    facets = reshape(typecast(facets(:), 'uint16'), 12*2, []);
    % Set the last bit to 0 (default) or supplied RGB
    facets(end+1,:) = options.facecolor;
    fwrite(fid, facets, 'uint16');
end

% Close the file
fclose(fid);
fprintf('Wrote %d facets\n',size(facets, 3));


%% Input handling subfunctions
function [faces, vertices, options] = parseInputs(varargin)
% Determine input type
if isstruct(varargin{1}) % stlwrite('file', FVstruct, ...)
    if ~all(isfield(varargin{1},{'vertices','faces'}))
        error( 'Variable p must be a faces/vertices structure' );
    end
    faces = varargin{1}.faces;
    vertices = varargin{1}.vertices;
    options = parseOptions(varargin{2:end});
    
elseif isnumeric(varargin{1})
    firstNumInput = cellfun(@isnumeric,varargin);
    firstNumInput(find(~firstNumInput,1):end) = 0; % Only consider numerical input PRIOR to the first non-numeric
    numericInputCnt = nnz(firstNumInput);
    
    options = parseOptions(varargin{numericInputCnt+1:end});
    switch numericInputCnt
        case 3 % stlwrite('file', X, Y, Z, ...)
            % Extract the matrix Z
            Z = varargin{3};
            
            % Convert scalar XY to vectors
            ZsizeXY = fliplr(size(Z));
            for i = 1:2
                if isscalar(varargin{i})
                    varargin{i} = (0:ZsizeXY(i)-1) * varargin{i};
                end                    
            end
            
            % Extract X and Y
            if isequal(size(Z), size(varargin{1}), size(varargin{2}))
                % X,Y,Z were all provided as matrices
                [X,Y] = varargin{1:2};
            elseif numel(varargin{1})==ZsizeXY(1) && numel(varargin{2})==ZsizeXY(2)
                % Convert vector XY to meshgrid
                [X,Y] = meshgrid(varargin{1}, varargin{2});
            else
                error('stlwrite:badinput', 'Unable to resolve X and Y variables');
            end
            
            % Convert to faces/vertices
            if strcmp(options.triangulation,'delaunay')
                faces = delaunay(X,Y);
                vertices = [X(:) Y(:) Z(:)];
            else
                if ~exist('mesh2tri','file')
                    error('stlwrite:missing', '"mesh2tri" is required to convert X,Y,Z matrices to STL. It can be downloaded from:\n%s\n',...
                        'http://www.mathworks.com/matlabcentral/fileexchange/28327')
                end
                [faces, vertices] = mesh2tri(X, Y, Z, options.triangulation);
            end
            
        case 2 % stlwrite('file', FACES, VERTICES, ...)
            faces = varargin{1};
            vertices = varargin{2};
            
        otherwise
            error('stlwrite:badinput', 'Unable to resolve input types.');
    end
end

if ~isempty(options.facecolor) % Handle colour preparation
    facecolor = uint16(options.facecolor);
    %Set the Valid Color bit (bit 15)
    c0 = bitshift(ones(size(faces,1),1,'uint16'),15);
    %Red color (10:15), Blue color (5:9), Green color (0:4)
    c0 = bitor(bitshift(bitand(2^6-1, facecolor(:,1)),10),c0);
    c0 = bitor(bitshift(bitand(2^11-1, facecolor(:,2)),5),c0);
    c0 = bitor(bitand(2^6-1, facecolor(:,3)),c0);
    options.facecolor = c0;    
else
    options.facecolor = 0;
end

function options = parseOptions(varargin)
IP = inputParser;
IP.addParamValue('mode', 'binary', @ischar)
IP.addParamValue('title', sprintf('Created by stlwrite.m %s',datestr(now)), @ischar);
IP.addParamValue('triangulation', 'delaunay', @ischar);
IP.addParamValue('facecolor',[], @isnumeric)
IP.parse(varargin{:});
options = IP.Results;

function [F,V]=mesh2tri(X,Y,Z,tri_type)
% function [F,V]=mesh2tri(X,Y,Z,tri_type)
% 
% Available from http://www.mathworks.com/matlabcentral/fileexchange/28327
% Included here for convenience. Many thanks to Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 15/07/2010
%------------------------------------------------------------------------

[J,I]=meshgrid(1:1:size(X,2)-1,1:1:size(X,1)-1);

switch tri_type
    case 'f'%Forward slash
        TRI_I=[I(:),I(:)+1,I(:)+1;  I(:),I(:),I(:)+1];
        TRI_J=[J(:),J(:)+1,J(:);   J(:),J(:)+1,J(:)+1];
        F = sub2ind(size(X),TRI_I,TRI_J);
    case 'b'%Back slash
        TRI_I=[I(:),I(:)+1,I(:);  I(:)+1,I(:)+1,I(:)];
        TRI_J=[J(:)+1,J(:),J(:);   J(:)+1,J(:),J(:)+1];
        F = sub2ind(size(X),TRI_I,TRI_J);
    case 'x'%Cross
        TRI_I=[I(:)+1,I(:);  I(:)+1,I(:)+1;  I(:),I(:)+1;    I(:),I(:)];
        TRI_J=[J(:),J(:);    J(:)+1,J(:);    J(:)+1,J(:)+1;  J(:),J(:)+1];
        IND=((numel(X)+1):numel(X)+prod(size(X)-1))';
        F = sub2ind(size(X),TRI_I,TRI_J);
        F(:,3)=repmat(IND,[4,1]);
        Fe_I=[I(:),I(:)+1,I(:)+1,I(:)]; Fe_J=[J(:),J(:),J(:)+1,J(:)+1];
        Fe = sub2ind(size(X),Fe_I,Fe_J);
        Xe=mean(X(Fe),2); Ye=mean(Y(Fe),2);  Ze=mean(Z(Fe),2);
        X=[X(:);Xe(:)]; Y=[Y(:);Ye(:)]; Z=[Z(:);Ze(:)];
end

V=[X(:),Y(:),Z(:)];
end
end
end
end


function varargout = surf2solid(varargin)
%SURF2SOLID  Convert a thin surface to a closed triangulated solid volume.
%
%   SOLID_FV = SURF2SOLID(FV,...) takes in a triangulated patch defined by
%   FV (a structure with fields 'vertices' and 'faces'), and returns a
%   solid patch SOLID_FV closed by options (described below).
%
%   SOLID_FV = SURF2SOLID(F, V,...) takes faces and vertices separately.
%
%   [F,V] = SURF2SOLID(...) returns solid faces and vertices separately.
%
%   SURF2SOLID(...) with no output argument plots the 3 components
%   (orig-surface, side-walls, under-surface) to a new figure.
%
%   SURF2SOLID(X, Y, Z, ...) reads in surface data in X, Y, and Z matrices,
%   and triangulates this gridded data into a surface using triangulation
%   options specified below. Z must be a 2D matrix. X and Y can be 2D
%   matrices of the same size as Z, or vectors of length equal to SIZE(Z,2)
%   and SIZE(Z,1), respectively. If X or Y are scalar values, they are used
%   to specify the X and Y spacing between grid points.
%
%   SURF2SOLID(...,'PropertyName',VALUE,...) makes a solid volume from thin
%   surface using any of the following property/value options:
%
%   ELEVATION     - Extends the surface down to a flat base at the given
%                   (Z) elevation value. Useful for turning a thin
%                   elevation map into a solid block with a flat base. The
%                   ELEVATION value should be below the lowest (or above
%                   the highest) data point. If no other options are given,
%                   ELEVATION defaults to MIN(Z)-0.1*(MAX(Z)-MIN(Z)).
%                   Variable ELEVATION may also be given per-point, via a
%                   2D matrix (the same size as Z for X,Y,Z style input) or
%                   a 1D array (with length equal to the number of vertices
%                   given in face/vertex input).
%
%   THICKNESS     - Value to offset the given thin surface to make a
%                   thickened solid slab. Each node on the surface will be
%                   projected along its normal direction by thickness. When
%                   negative thickness is given, offset will be away from
%                   face normal direction. Variable thickness can also be
%                   specified via a 2D matrix (of same size as Z, for X,Y,Z
%                   input) or an N-by-1 array of thicknesses (where N is
%                   the number of vertices in the thin surface)
%
%   TRIANGULATION - When used with gridded data, TRIANGULATION is either:
%                    'delaunay'  - (default) Delaunay triangulation of X, Y
%                    'f'         - Forward slash division of grid quads
%                    'b'         - Back slash division of quadrilaterals
%                    'x'         - Cross division of quadrilaterals
%                   Note that 'f', 'b', or 'x' triangulations use an
%                   inbuilt version of FEX entry 28327, "mesh2tri". 'x'
%                   style triangulation cannot be used with variable
%                   ELEVATION or THICKNESS parameters.
%
% Note 1: Currently surf2solid will return a closed surface with face
% normals pointing "out". With user feedback, I'd be happy to change this
% behaviour to either "in" or "unchanged from input direction".
% Note 2: If a single ELEVATION value is specified (i.e., flat base), the
% resulting patch will have minimal triangles on the flat base to reduce
% patch/file size.
%
%   Example (shows both THICKNESS and ELEVATION forms):
%     n = 30;
%     [X,Y] = meshgrid(linspace(0,1,2*n+1));
%     L = (40/51/0.9)*membrane(1,n);
%     figure, subplot(2,2,[1 3]), title 'Thin surface'
%     surf(X,Y,L,'EdgeColor','none'); colormap pink; axis image; camlight
%     subplot(2,2,2), title 'Block elevation'
%     surf2solid(X,Y,L,'elevation',min(L(:))-0.05); axis image; camlight; camlight
%     subplot(2,2,4), title 'Thickness'
%     surf2solid(X,Y,L,'thickness',-0.1); axis image; camlight;
%
%   Original idea adapted from Paul Kassebaum's blog post
%   http://blogs.mathworks.com/community/2013/06/20/paul-prints-the-l-shaped-membrane/
%   Many thanks to Paul for his further input and improvements.
%
%   Author: Sven Holcombe, 07-20-2013

% 1.0 (2013-07) Original
% 1.01(2013-07) Set other-faces as opposite norm direction to input faces
%               Ensured optional per-node thickness input is allowed
% 1.1 (2013-07) Added per-node elevation and thickness, minimal-flat-base
%               file size, and default face orientation.

% Get faces, vertices, and user-defined options for writing
[F, V, options] = parseInputs(varargin{:});

% Get the latest triangulation class. 2013a brought in "triangulation" with
% some nice options. Earlier versions must use "TriRep".
if exist('triangulation','class')
    T = triangulation(F,V);
    options.oldVersion = false;
else
    T = TriRep(F,V); %#ok
    options.oldVersion = true;
end
% Extract boundary edges from input surface. These will connect to walls.
boundEdges = T.freeBoundary;
boundVerts = boundEdges([1:end 1],1);

% Define "other" faces opposite to input faces. These will either sit at
% given Z-elevations, or they will be offset from input faces by thickness.
if ~isempty(options.elevation)  % ELEVATION was specified
    % If scalar elevation was given, it will be assigned here. If variable
    % elevation was given, it will also be assigned.
    V_extrude = V;
    V_extrude(:,3) = options.elevation(:);
else                            % THICKNESS was specified
    % Get vertex normals the hard or the easy way
    if options.oldVersion
        facets = V';
        facets = permute(reshape(facets(:,F'), 3, 3, []),[2 1 3]);
        allEdgeVecs = facets([2 3 1],:,:) - facets(:,:,:);
        allFacetNormals =  bsxfun(@times, allEdgeVecs(1,[2 3 1],:), allEdgeVecs(2,[3 1 2],:)) - ...
            bsxfun(@times, allEdgeVecs(2,[2 3 1],:), allEdgeVecs(1,[3 1 2],:));
        allFacetNormals = bsxfun(@rdivide, allFacetNormals, sqrt(sum(allFacetNormals.^2,2)));
        facesByVertex = T.vertexAttachments;
        Vnormals = zeros(size(V));
        for i = 1:length(facesByVertex)
            Vnormals(i,:) =  mean(allFacetNormals(:,:,facesByVertex{i}),3);
        end
        Vnormals = bsxfun(@rdivide, Vnormals, sqrt(sum(Vnormals.^2,2)));
    else
        Vnormals = T.vertexNormal;
    end
    % Extrudes by thickness in each normal direction.
    % bsxfun is used in case the user wants to supply variables offsets by
    % each vertex.
    V_extrude = V + bsxfun(@times, Vnormals, options.thickness(:));
end

% If a scalar elevation was supplied, we can save file size (or triangle
% count) and define the base by only its boundary vertices.
if isscalar(options.elevation)
    V_wall = [V(boundVerts,:); V_extrude(boundVerts,:)];
else % thickness was specified
    V_wall = [V(boundVerts,:); V_extrude(boundVerts,:)];
end

% Number of wall vertices on each surface (nwv).
nwv = length(V_wall)/2;
% Allocate memory for wallFaces.
F_wall = zeros(2*(nwv-1),3);
% Define the faces.
for k = 1:nwv-1
    F_wall(k      ,:) = [k+1  ,k      ,k+nwv];
    F_wall(k+nwv-1,:) = [k+nwv,k+1+nwv,k+1];
end

% Let's use the first vertex to test if faces are pointed in/out
testNormal = cross(...
    V(F(1,2),:) - V(F(1,1),:),...
    V(F(1,3),:) - V(F(1,1),:));
if ~isempty(options.elevation)
    firstVertZoffset = options.elevation(1) - V(1,3);
else
    firstVertZoffset = options.thickness(1);
end
% If the first face is pointing in the same direction as the direction of
% offset, we will have all faces pointing "in". We want them pointing "out"
if sign(testNormal(3)) == sign(firstVertZoffset)
    F = fliplr(F);
end

if isscalar(options.elevation)     % SCALAR ELEVATION was specified
    % Each boundary vertex forms a triangle with its right-hand neighbor
    % and a newly formed centroid of the extruded face.
    n = length(boundVerts);
    V_extrude = [mean(V_extrude,1); V_extrude(boundVerts,:)];% prepend centroid.
    % Ensure extruded faces are properly oriented.
    testNormal = cross(...
        V_extrude(2,:)-V_extrude(1,:),...
        V_extrude(3,:)-V_extrude(1,:));
    if sign(testNormal(3)) == sign(options.elevation)
        F_extrude = [ones(n-1,1),(2:n)',[(3:n)';2]];
    else
        F_extrude = [[(3:n)';2],(2:n)',ones(n-1,1)];
    end
else % Thickness or variable elevation was specified
    % Simply ensure the extruded faces are oriented opposite the originals
    F_extrude = fliplr(F);
end

% Compile 3 sets of faces together
allVertices = [V; V_wall; V_extrude];
allFaces = [F;                           % Use original faces
    F_wall+size(V,1);                    % Add wall faces
    F_extrude+size(V,1)+size(V_wall,1)]; % Add opposite faces (flipped)

% Ouput based on requested variables.
varargout = {};
if nargout == 0
  %figure;
  hold on;
  view(3);
  axis vis3d;
  patch('Faces'    ,F        ,...
        'Vertices' ,V        ,...
        'FaceColor','r'      );
  patch('Faces'    ,F_wall   ,...
        'Vertices' ,V_wall   ,...
        'FaceColor','g'      );
  patch('Faces'    ,F_extrude,...
        'Vertices' ,V_extrude,...
        'FaceColor','b'      );
  hold off;
%   set(gca,'zlim',[min(allVertices(:,3)),max(allVertices(:,3))]);
elseif nargout == 1
  varargout = {struct('faces',allFaces,'vertices',allVertices)};
elseif nargout >= 2
  varargout = {allFaces, allVertices};
end


%% Input handling subfunctions
function [faces, vertices, options] = parseInputs(varargin)
% Determine input type
if isstruct(varargin{1}) % surf2solid(FVstruct, ...)
    if ~all(isfield(varargin{1},{'vertices','faces'}))
        error( 'Variable p must be a faces/vertices structure' );
    end
    faces = varargin{1}.faces;
    vertices = varargin{1}.vertices;
    options = parseOptions(varargin{2:end});
    
elseif isnumeric(varargin{1})
    firstNumInput = cellfun(@isnumeric,varargin);
    firstNumInput(find(~firstNumInput,1):end) = 0; % Only consider numerical input PRIOR to the first non-numeric
    numericInputCnt = nnz(firstNumInput);
    
    options = parseOptions(varargin{numericInputCnt+1:end});
    switch numericInputCnt
        case 3 % surf2solid(X, Y, Z, ...)
            % Extract the matrix Z
            Z = varargin{3};
            
            % Convert scalar XY to vectors
            ZsizeXY = fliplr(size(Z));
            for i = 1:2
                if isscalar(varargin{i})
                    varargin{i} = (0:ZsizeXY(i)-1) * varargin{i};
                end                    
            end
            
            % Extract X and Y
            if isequal(size(Z), size(varargin{1}), size(varargin{2}))
                % X,Y,Z were all provided as matrices
                [X,Y] = varargin{1:2};
            elseif numel(varargin{1})==ZsizeXY(1) && numel(varargin{2})==ZsizeXY(2)
                % Convert vector XY to meshgrid
                [X,Y] = meshgrid(varargin{1}, varargin{2});
            else
                error('surf2solid:badinput', 'Unable to resolve X and Y variables');
            end
            
            % Convert to faces/vertices
            if strcmp(options.triangulation,'delaunay')
                faces = delaunay(X,Y);
                vertices = [X(:) Y(:) Z(:)];
            else
                if ~exist('mesh2tri','file')
                    error('surf2solid:missing', '"mesh2tri" is required to convert X,Y,Z matrices to STL. It can be downloaded from:\n%s\n',...
                        'http://www.mathworks.com/matlabcentral/fileexchange/28327')
                end
                [faces, vertices] = mesh2tri(X, Y, Z, options.triangulation);
            end
            
        case 2 % surf2solid(FACES, VERTICES, ...)
            faces = varargin{1};
            vertices = varargin{2};
            
        otherwise
            error('surf2solid:badinput', 'Unable to resolve input types.');
    end
end
% Ensure *some* information is there to make a thickness
if isempty(options.thickness) && isempty(options.elevation)
    options.elevation = min(vertices(:,3)) - 0.1 * (max(vertices(:,3))-min(vertices(:,3)));
end

function options = parseOptions(varargin)
IP = inputParser;
IP.addParamValue('triangulation', 'delaunay', @ischar);
IP.addParamValue('elevation',[],@isnumeric)
IP.addParamValue('thickness',[],@isnumeric)
IP.parse(varargin{:});
options = IP.Results;

%% INCLUDED FUNCTIONS %%
function [F,V]=mesh2tri(X,Y,Z,tri_type)
% function [F,V]=mesh2tri(X,Y,Z,tri_type)
% 
% Available from http://www.mathworks.com/matlabcentral/fileexchange/28327
% Included here for convenience. Many thanks to Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 15/07/2010
%------------------------------------------------------------------------
[J,I]=meshgrid(1:1:size(X,2)-1,1:1:size(X,1)-1);
switch tri_type
    case 'f'%Forward slash
        TRI_I=[I(:),I(:)+1,I(:)+1;  I(:),I(:),I(:)+1];
        TRI_J=[J(:),J(:)+1,J(:);   J(:),J(:)+1,J(:)+1];
        F = sub2ind(size(X),TRI_I,TRI_J);
    case 'b'%Back slash
        TRI_I=[I(:),I(:)+1,I(:);  I(:)+1,I(:)+1,I(:)];
        TRI_J=[J(:)+1,J(:),J(:);   J(:)+1,J(:),J(:)+1];
        F = sub2ind(size(X),TRI_I,TRI_J);
    case 'x'%Cross
        TRI_I=[I(:)+1,I(:);  I(:)+1,I(:)+1;  I(:),I(:)+1;    I(:),I(:)];
        TRI_J=[J(:),J(:);    J(:)+1,J(:);    J(:)+1,J(:)+1;  J(:),J(:)+1];
        IND=((numel(X)+1):numel(X)+prod(size(X)-1))';
        F = sub2ind(size(X),TRI_I,TRI_J);
        F(:,3)=repmat(IND,[4,1]);
        Fe_I=[I(:),I(:)+1,I(:)+1,I(:)]; Fe_J=[J(:),J(:),J(:)+1,J(:)+1];
        Fe = sub2ind(size(X),Fe_I,Fe_J);
        Xe=mean(X(Fe),2); Ye=mean(Y(Fe),2);  Ze=mean(Z(Fe),2);
        X=[X(:);Xe(:)]; Y=[Y(:);Ye(:)]; Z=[Z(:);Ze(:)];
end
V=[X(:),Y(:),Z(:)];
end
end
end
end
