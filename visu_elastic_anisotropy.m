%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is supplementary material to visualize elastic anisotropy
% based on the following article.
%
% Manuscript title: Visualising elastic anisotropy: theoretical background 
% and computational implementation
% Manuscript Authors: Joachim Nordmann, Marcus Aßmus, Holm Altenbach
% DOI: 10.1007/s00161-018-0635-9
%
% Inputs: Cijkl (material parameters of the material chosen)
%
% Work address :    Chair of Engineering Mechanics
%                   Institute of Mechanics
%                   Faculty of Mechanical Engineering
%                   Otto von Guericke University
%                   Universitaetsplatz 2                    
%                   39108 Magdeburg 
%                   Germany
% Website:          https://www.ifme.ovgu.de/ltm
% Contact:          joachim.nordmann@ovgu.de
%                   marcus.assmus@0vgu.de
%                   holm.altenbach@ovgu.de
%
% December 2017; Last revision: 28-November-2019
%
% When refering to the script in publications please cite as follows:
%
% @article{noasal2018,
%  author = {J. Nordmann and M. A{\ss}mus and H. Altenbach},
%  title = {Visualising Elastic Anisotropy: Theoretical Background and 
%           Computational Implementation},
%  journal = {Continuum Mechanics and Thermodynamics},
%  volume = {30},
%  number = {4},
%  year = {2018},
%  pages = {689--708},
%  doi = {10.1007/s00161-018-0635-9}
% }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
%% Input
% As input only the stiffness matrix in the correct Voigt Notation and the 
% division of azimuth angle are required. At the end this code produces 4 
% plots for the different material parameters.  
% The azimuth angle is in the range of 0 to 2*pi, the elongation angle goes
% from 0 to pi and the rotation angle psi is in the range of 0 to 2*pi too.
%
% Scheme for the stiffness matrix
%         
% C = [         C11         C12         C13 sqrt(2)*C14 sqrt(2)*C15 sqrt(2)*C16;
%               C12         C22         C23 sqrt(2)*C24 sqrt(2)*C25 sqrt(2)*C26; 
%               C13         C23         C33 sqrt(2)*C34 sqrt(2)*C35 sqrt(2)*C36;
%       sqrt(2)*C14 sqrt(2)*C24 sqrt(2)*C34       2*C44       2*C45       2*C46;
%       sqrt(2)*C15 sqrt(2)*C25 sqrt(2)*C35       2*C45       2*C55       2*C56;
%       sqrt(2)*C16 sqrt(2)*C26 sqrt(2)*C36       2*C46       2*C56       2*C66 ];


% Division / Elements of the azimuth angle
n         = 70; 
% Division / Elements of the elongation angle
m         = n/2;
% Division / Elements of the rotation angle
o         = 360;
% Silicon in GPa (Values are from [18])
C     = [ 165.7  63.9  63.9  0    0    0;
           63.9 165.7  63.9  0    0    0;
           63.9  63.9 165.7  0    0    0;
            0     0     0  2*79.6  0    0;
            0     0     0    0   2*79.6  0;
            0     0     0    0    0   2*79.6 ];

% Test materials in Voigt-Notation

% EVA+Si in GPa (Values are from [32])
% C = [ 495.29  18.2  13.2  0    0        0;
%        18.2  495.29 13.2  0    0        0;
%        13.2   13.2  22.68 0    0        0;
%         0      0     0    2*3.46 0        0;
%         0      0     0    0    2*3.46     0;
%         0      0     0    0    0    2*40000 ]/1000; 
        
% Cooper in GPa (Values are from [18])
% C     = [ 168   121.4 121.4    0      0      0;
%           121.4 168   121.4    0      0      0;
%           121.4 121.4 168      0      0      0;
%             0     0     0   2*75.4    0      0;
%             0     0     0      0   2*75.4    0;
%             0     0     0      0      0   2*75.4 ];

% Gold in GPa (Values are from [18])
% C     = [ 185 158 158  0    0    0;
%           158 185 158  0    0    0;
%           158 158 185  0    0    0;
%             0   0   0 2*39.7  0    0;
%             0   0   0  0   2*39.7  0;
%             0   0   0  0    0   2*39.7 ];

% Magnesium in GPa (Values are from [18])
% C     = [ 56.49 23.16 18.10    0      0        0;
%           23.16 56.49 18.10    0      0        0;
%           18.10 18.10 58.73    0      0        0;
%             0     0     0   2*16.81   0        0;
%             0     0     0      0   2*16.81     0;
%             0     0     0      0      0   (56.49-23.16) ];

% Feldspar in GPa (Values are from [18])
% C     = [ 61.9  43.4  36.8 -sqrt(2)*10       0        0;
%           43.4 158.3  21.8 -sqrt(2)* 1.8     0        0;
%           36.8  21.8 100.2 -sqrt(2)*12.1     0        0;
%            0     0   -sqrt(2)*12.1  2*14.1   0        0;
%            0     0   -sqrt(2)* 1.8   0    2*20.3  -2* 2.3;
%            0     0   -sqrt(2)*10     0   -2* 2.3   2*36 ];

% Pyrite in GPa (Values are from [18])
% C     = [ 361.88 -47.96 -47.96     0        0        0;
%           -47.96 361.88 -47.96     0        0        0;
%           -47.96 -47.96 361.88     0        0        0;
%             0   0  0           2*105.49     0        0;
%             0   0  0               0    2*105.49     0;
%             0   0  0               0        0    2*105.49 ];

% Weissblech in GPa (Values are from [18])
% C     = [ 76.19 71.10  67.68     0     0     0;
%           71.10 76.19  67.68     0     0     0;
%           67.68 67.68 116.24     0     0     0;
%            0     0      0     2*17.04  0     0;
%            0     0      0      0    2*17.04  0;
%            0     0      0      0     0    2*19.8 ];

% Isotropic testmaterial. Deriving the compliance matrix in line 172 must
% be deactivated by using this material.
% Material parameters are defined randomly. Just for debugging. Solutions
% are spheres with a radius of Y or nue.
% Y = 100000;
% nue = 0.25;
% S = [   1/Y, -nue/Y, -nue/Y,    0,           0,           0;
%      -nue/Y,    1/Y, -nue/Y,    0,           0,           0;
%      -nue/Y, -nue/Y,    1/Y,    0,           0,           0;
%         0,      0,      0,   2*(1+nue)/Y/2,    0,           0;
%         0,      0,      0,      0,        2*(1+nue)/Y/2,    0;
%         0,      0,      0,      0,           0,        2*(1+nue)/Y/2];

%% Definition of the Memory
% Angles
phi       = linspace(0,  pi,m);
theta     = linspace(0,2*pi,n);
psi       = linspace(0,2*pi,o);
% Young's modulus
E         = zeros(m,n);
% Bulk modulus
K         = zeros(m,n);
% Normal vector
nn        = zeros(3,o);
% Poisson ratio
nu        = zeros(1,o);
nu_max    = zeros(m,n);
nu_min    = zeros(m,n);
nu_avg    = zeros(m,n);
% Shear modulus
G         = zeros(1,o);
G_max     = zeros(m,n);
G_min     = zeros(m,n);
G_avg     = zeros(m,n);
% Vectors for the cartesian coordinates 
x_E       = zeros(m,n);
y_E       = zeros(m,n);
z_E       = zeros(m,n);
x_K       = zeros(m,n);
y_K       = zeros(m,n);
z_K       = zeros(m,n);
x_nu_max  = zeros(m,n);
y_nu_max  = zeros(m,n);
z_nu_max  = zeros(m,n);
x_nu_min  = zeros(m,n);
y_nu_min  = zeros(m,n);
z_nu_min  = zeros(m,n);
x_nu_avg  = zeros(m,n);
y_nu_avg  = zeros(m,n);
z_nu_avg  = zeros(m,n);
x_G_max   = zeros(m,n);
y_G_max   = zeros(m,n);
z_G_max   = zeros(m,n);
x_G_min   = zeros(m,n);
y_G_min   = zeros(m,n);
z_G_min   = zeros(m,n);
x_G_avg   = zeros(m,n);
y_G_avg   = zeros(m,n);
z_G_avg   = zeros(m,n);
% Calculation of the compliance matrix
S     = C^-1; 
%% Calculation of Young's and bulk modulus
% Definition of the unit tensor in Voigt-notation
U_V   = [ 1;
          1;
          1;
          0;
          0;
          0 ];
for i = 1:1:n
    for j = 1:1:m
    % Direction vector with radius r=1    
    d      = [ sin(phi(1,j))*cos(theta(1,i));
               sin(phi(1,j))*sin(theta(1,i));
               cos(phi(1,j)) ];      
    % Calculation of the tensor D, d dyadic d
    D      = d*d';
    % Transformation into Voigt-scheme
    d_V    = [ D(1,1);
               D(2,2);
               D(3,3);
               sqrt(2)*D(2,3);
               sqrt(2)*D(1,3);
               sqrt(2)*D(1,2) ];
    % Calculation of the radius 
    % Young's modulus
    E(j,i) = 1/(d_V'*(S*d_V));
    % Bulk modulus
    K(j,i) = 1/(3*U_V'*(S*d_V));
    end 
end
%% Calculation of Poisson ratio and shear modulus
% Activate parallel processing
poolobj = parpool;
for i = 1:1:n
    for j = 1:1:m
    % Direction vector d with radius r=1    
    d      = [ sin(phi(1,j))*cos(theta(1,i));
               sin(phi(1,j))*sin(theta(1,i));
               cos(phi(1,j)) ];
    % Transformation into tensor D, d dyadic d
    D      = d*d';
    % Transformation into Voigt-scheme
    d_V    = [ D(1,1);
               D(2,2);
               D(3,3);
               sqrt(2)*D(2,3);
               sqrt(2)*D(1,3);
               sqrt(2)*D(1,2) ];
    %       
               for k = 1:1:o
               % Normal vector with radius r=1    
               nn(1,k) = -cos(phi(1,j))*cos(theta(1,i))*cos(psi(1,k))+...
                          sin(theta(1,i))*sin(psi(1,k));
               nn(2,k) = -cos(phi(1,j))*sin(theta(1,i))*cos(psi(1,k))+...
                         -cos(theta(1,i))*sin(psi(1,k));
               nn(3,k) =  sin(phi(1,j))*cos(psi(1,k));  
                       % Very small vector elements are set to exact zero
                       for q = 1:1:3
                       if abs(nn(q,k)) <= 1e-6
                       nn(q,k) = 0;
                       else
                       end
                       end
               % Transformation of the normal vector into tensor N,
               % nn dyadic nn
               N   = nn(:,k)*nn(:,k)';
               % Transformation into Voigt-scheme
               n_V = [ N(1,1);
                       N(2,2);
                       N(3,3);
                       sqrt(2)*N(2,3);
                       sqrt(2)*N(1,3);
                       sqrt(2)*N(1,2) ];  
               % Calculation of Tensor M
               M   = sqrt(2)/2*(d*nn(:,k)'+nn(:,k)*d');
               % Transformation into Voigt-scheme
               m_V = [ M(1,1);
                       M(2,2);
                       M(3,3);
                       sqrt(2)*M(2,3);
                       sqrt(2)*M(1,3);
                       sqrt(2)*M(1,2) ]; 
               % Calculation of Poisson ratio for every psi     
               nu(1,k) = -E(j,i)*d_V'*(S*n_V);
               % Calculation of shear modulus for every psi
               G(1,k)  = 1/(2*m_V'*(S*m_V));
               end
    % Using of the criteria for Poisson ratio             
    nu_max(j,i) = max(nu);
    nu_min(j,i) = min(nu);
    nu_avg(j,i) = sum(nu)/o;
    % Using of the criteria for shear modulus              
    G_max(j,i)  = max(G);
    G_min(j,i)  = min(G);
    G_avg(j,i)  = sum(G)/o;
    end 
end
% Deactivate parallel processing
delete(poolobj);
%% Transformation of spherical coordiantes into cartesian coordinates 
for i = 1:1:n
    for j = 1:1:m
    % Young's modulus
    x_E(j,i)      = E(j,i)*sin(phi(1,j))*cos(theta(1,i));
    y_E(j,i)      = E(j,i)*sin(phi(1,j))*sin(theta(1,i));
    z_E(j,i)      = E(j,i)*cos(phi(1,j));
    % Bulk modulus
    x_K(j,i)      = K(j,i)*sin(phi(1,j))*cos(theta(1,i));
    y_K(j,i)      = K(j,i)*sin(phi(1,j))*sin(theta(1,i));
    z_K(j,i)      = K(j,i)*cos(phi(1,j));    
    % Poisson ratio 
    % 1. Maximum / 2. Minimum / 3. Average
    x_nu_max(j,i) = nu_max(j,i)*sin(phi(1,j))*cos(theta(1,i));
    y_nu_max(j,i) = nu_max(j,i)*sin(phi(1,j))*sin(theta(1,i));
    z_nu_max(j,i) = nu_max(j,i)*cos(phi(1,j));
    x_nu_min(j,i) = nu_min(j,i)*sin(phi(1,j))*cos(theta(1,i));
    y_nu_min(j,i) = nu_min(j,i)*sin(phi(1,j))*sin(theta(1,i));
    z_nu_min(j,i) = nu_min(j,i)*cos(phi(1,j));
    x_nu_avg(j,i) = nu_avg(j,i)*sin(phi(1,j))*cos(theta(1,i));
    y_nu_avg(j,i) = nu_avg(j,i)*sin(phi(1,j))*sin(theta(1,i));
    z_nu_avg(j,i) = nu_avg(j,i)*cos(phi(1,j));
    % Shear modulus
    % 1. Maximum / 2. Minimum / 3. Average
    x_G_max(j,i)  = G_max(j,i)*sin(phi(1,j))*cos(theta(1,i));
    y_G_max(j,i)  = G_max(j,i)*sin(phi(1,j))*sin(theta(1,i));
    z_G_max(j,i)  = G_max(j,i)*cos(phi(1,j));
    x_G_min(j,i)  = G_min(j,i)*sin(phi(1,j))*cos(theta(1,i));
    y_G_min(j,i)  = G_min(j,i)*sin(phi(1,j))*sin(theta(1,i));
    z_G_min(j,i)  = G_min(j,i)*cos(phi(1,j));
    x_G_avg(j,i)  = G_avg(j,i)*sin(phi(1,j))*cos(theta(1,i));
    y_G_avg(j,i)  = G_avg(j,i)*sin(phi(1,j))*sin(theta(1,i));
    z_G_avg(j,i)  = G_avg(j,i)*cos(phi(1,j));
    end
end
%% Plot of Young's modulus and Bulk modulus
% Colormap
map1 = [0.3, 0, 0
        0.4, 0, 0
        0.5, 0, 0
        0.6, 0, 0
        0.8, 0, 0
        1.0, 0, 0];
% Young's modulus    
a = figure(1);
set(a, 'Units', 'centimeters', 'Position', [15, 10, 12, 10]);
% Plot of the surface in cartesian coordiantes. The fourth argument 
% "E" respectively "K" is for the correct color.
hold on
surf(x_E, y_E, z_E, E,'FaceColor','interp','FaceAlpha',1.0,'LineWidth',...
     0.01,'EdgeAlpha',0.2,'EdgeColor','w');
hold off
grid off
box off
set(gca,...
        'XTick',[],...
        'YTick',[],...
        'ZTick',[],...
        'XColor','w',...
        'YColor','w',...
        'ZColor','w'); 
colormap(map1);
colormap(flipud(colormap));
c = colorbar('eastoutside','Position',[0.87 0.27 .02 .5],'FontSize',10); 
c.Label.String = '{\it{Y}} [GPa]';
c.Label.FontName = 'Times New Roman';
% Bulk modulus
b = figure(2);
set(b, 'Units', 'centimeters', 'Position', [15, 10, 12, 10]);
hold on
surf(x_K, y_K, z_K, K,'FaceAlpha',1.0,'EdgeAlpha',0.2,'LineWidth',....
     0.01,'EdgeColor','w');
hold off
grid off
box off
set(gca,...
        'XTick',[],...
        'YTick',[],...
        'ZTick',[],...
        'XColor','w',...
        'YColor','w',...
        'ZColor','w'); 
colormap(map1);
h = colorbar('eastoutside','Position',[0.8 0.3 .02 .4],'FontSize',10); 
set(h,'YTickLabel',cellstr(num2str(reshape(get(h,'YTick'),[],1),'%0.2f')));
h.Label.String = '{\it{K}} [GPa]';
h.Label.FontName = 'Times New Roman';
%% Plot of Poisson ratio
f = figure(3);
set(f, 'Units', 'centimeters', 'Position', [15, 10, 10, 8]);
hold on
h1 = surf(x_nu_max, y_nu_max, z_nu_max, nu_max,'FaceColor','interp',...
          'FaceAlpha',1.0,'LineWidth',0.01,'EdgeAlpha',0.2,'EdgeColor','w');
% h2 = surf(x_nu_min, y_nu_min, z_nu_min, nu_min,'FaceColor','interp',...
%           'FaceAlpha',1.0,'LineWidth',0.01,'EdgeAlpha',0.2,'EdgeColor','w');
% h3 = surf(x_nu_avg, y_nu_avg, z_nu_avg, nu_avg,'FaceColor','interp',...
%           'FaceAlpha',1.0,'LineWidth',0.01,'EdgeAlpha',0.2,'EdgeColor','w');
hold off
grid off
box off
set(gca,...
        'XTick',[],...
        'YTick',[],...
        'ZTick',[],...
        'XColor','w',...
        'YColor','w',...
        'ZColor','w');
colormap(map1);
g = colorbar('southoutside','Position',[0.32 0.27 .4 .02],'FontSize',10);
g.Label.String = '{\it{\nu}}_{max} [-]';
g.Label.FontName = 'Times New Roman';
%% Plot of shear modulus   
t = figure(4);
set(t, 'Units', 'centimeters', 'Position', [15, 10, 10, 8]);
hold on

h4 = surf(x_G_max, y_G_max, z_G_max, G_max,'FaceColor','interp',...
          'FaceAlpha',1.0,'LineWidth',0.01,'EdgeAlpha',0.2,'EdgeColor','w');
% h5 = surf(x_G_min, y_G_min, z_G_min, G_min,'FaceColor','interp',...
%           'FaceAlpha',1.0,'LineWidth',0.01,'EdgeAlpha',0.2,'EdgeColor','w');
% h6 = surf(x_G_avg, y_G_avg, z_G_avg, G_avg,'FaceColor','interp',...
%           'FaceAlpha',1.0,'LineWidth',0.01,'EdgeAlpha',0.2,'EdgeColor','w');
hold off
grid off
box off
set(gca,...
        'XTick',[],...
        'YTick',[],...
        'ZTick',[],...
        'XColor','w',...
        'YColor','w',...
        'ZColor','w');
colormap(map1);
d = colorbar('southoutside','Position',[0.33 0.20 .4 .02],'FontSize',10);
d.Label.String = '{\it{G}}_{max} [GPa]';
d.Label.FontName = 'Times New Roman';                  