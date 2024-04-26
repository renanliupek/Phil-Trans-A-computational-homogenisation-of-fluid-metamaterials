%% Renan Liupekevicius Carnielli TU/e
% 01-06-2022
% computational homogenization of fluid medium with internal dynamics
clear; close all;


% TOPOLOGY: Fluid at top-bottom-right-left 


% code lines in the sandwiched by comments of 
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% something you may change.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------



%% IMPORT MESH

% select unit cell design from  'meshes' folder
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
design = 'Liang4d1_willis';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

%% CONVERT TO QUADRATIC MESH

% read mesh text file comsol 5.4
%--------------------------------------------------------------------------
disp('codeblock: READ & CONVERT LINEAR TO QUADRATIC ELEMENT MESH 5.4')
% read text file
l   = read_mphtxt_54(design);

% convert to quadratic element
mesh=mesh2d_lin2qua_uc(l);
%--------------------------------------------------------------------------


% copy to local variables (please double click on 'mesh' in workspace to
% see the meaning of each cell)
%--------------------------------------------------------------------------
  x    = mesh.x{2,2};    % tensor form
  mx   = mesh.x{2,1};    % matrix form
  conn = mesh.elem{2,1}; % conn matrix of element
%--------------------------------------------------------------------------


% mesh parameters
%--------------------------------------------------------------------------
  m  = size(  conn, 1); % number of elements in the mesh
  n  = size(  x, 1);    % number of nodes in the mesh 
  fprintf('NUMBER OF MESH ELEMENTS: %d\n', m)
  fprintf('NUMBER OF MESH NODES: %d\n', n)
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% get unit cell dimensions ( UC net not be squared)
% unit cell size
  ax = 2*max(  mx(:,1));
  ay = 2*max(  mx(:,2));
  

if(max(mx(:,1))+min(mx(:,1))>1e-8); warning('RVE is not x centralized');end
if(max(mx(:,1))+min(mx(:,1))>1e-8); warning('RVE is not y centralized');end  

if(ax==ay); a=ax;end
%--------------------------------------------------------------------------


%% DESIGNS: LIST OF TAGS (MANUALLY ADD FOR EACH DESIGN)
disp('codeblock: LOAD DESIGN TAGS')
% MANUALLY import tags from comsol: check e.g. boundary
% list under 'boundary selection'. Same procedure for selecting tag-list of
% the fluid phase(s).

% For this article (P.Transactions A) only one fluid phase is considered.
multiple_fluid_phase = false;

% switch sscanf(design, '%c')
% 
%     case 'Kurzeja_squaregrain' % 
%     %----------------------------------------------------------------------
%     % fluid elements tag
%     list_fluid_tag = [1:180];
%     % fluid1 elements tag
%     list_fluid1_tag = [63 64 69 70 75 76 83 :86 97 :102 111 112 117 118];
%     % fluid elements tag
%     list_fluid2_tag = [1 :62 65 :68 71 :74 77 :82 87 :96 103 :110 113 :116 119 :180]; 
%     % multiples fluid phases
%     multiple_fluid_phase = true;
%     %----------------------------------------------------------------------
% 
%     case 'Kurzeja' % 
%     %----------------------------------------------------------------------
%     % fluid elements tag
%     list_fluid_tag = [1:16];
%     % fluid1 elements tag
%     list_fluid1_tag = [10 11];
%     % fluid elements tag
%     list_fluid2_tag = [1:9 12:16]; 
%     % multiples fluid phases
%     multiple_fluid_phase = true;
%     %----------------------------------------------------------------------
% 
%     case 'designXu_effective' % 
%     %----------------------------------------------------------------------
%     % fluid elements tag
%     list_fluid_tag = [1:12];
%     % fluid1 elements tag
%     list_fluid1_tag = [1 2 5:7 9 10 12];
%     % fluid elements tag
%     list_fluid2_tag = [3 4 8 11]; 
%     % multiples fluid phases
%     multiple_fluid_phase = true;
%     %----------------------------------------------------------------------
% 
%     case 'designXu_effective_modified' % 
%     %----------------------------------------------------------------------
%     % fluid elements tag
%     list_fluid_tag = [1:12];
%     % fluid1 elements tag
%     list_fluid1_tag = [1 2 5:7 9 10 12];
%     % fluid elements tag
%     list_fluid2_tag = [3 4 8 11]; 
%     % multiples fluid phases
%     multiple_fluid_phase = true;
%     %----------------------------------------------------------------------
% 
%     case 'designLiang4d1_effective' % 
%     %----------------------------------------------------------------------
%     % fluid elements tag
%     list_fluid_tag = [1:9];
%     % fluid1 elements tag
%     list_fluid1_tag = [1:4 6:9];
%     % fluid elements tag
%     list_fluid2_tag = [5]; 
%     % multiples fluid phases
%     multiple_fluid_phase = true;
%     %----------------------------------------------------------------------
% 
% 
%     case 'designLiu_50xcoiledeffective' % 
%     %----------------------------------------------------------------------
%     % fluid elements tag
%     list_fluid_tag = [1:12];
%     % fluid1 elements tag
%     list_fluid1_tag = [1 2 5:7 9 10 12];
%     % fluid elements tag
%     list_fluid2_tag = [3 4 8 11]; 
%     % multiples fluid phases
%     multiple_fluid_phase = true;
%     %----------------------------------------------------------------------
% 
% end

%% INITIAL DEFINITIONS: MATERIALS, TENSOR BASIS AND CANONICAL TENSORS
disp('codeblock: LOAD MATERIAL PROPERTIES')
%--------------------------------------------------------------------------
dim   = 2;  % problem dimension
thickout = 1;  % thickness (out-of-plane direction), in [m]
%--------------------------------------------------------------------------

if not(multiple_fluid_phase)
% fluid properties
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhof     = 1.225 ; % [kg/m3]
c        = 343 ;  % [m/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

elseif multiple_fluid_phase
    
    n_fluid_phases = 2;

    % define fluid phase parameters
    %----------------------------------------------------------------------
   
    % % fluid properties
    % %--------------------------------------------------------------------------
    % mat{1}.rhof     = 1.225  % [kg/m3]
    % mat{1}.c        = 343   % [m/s]
    % %-------------------------------------------------------------------------
    % 
    % 
    % % fluid properties (effective medium)
    % %--------------------------------------------------------------------------
    % mat{2}.rhof     = 3402.8   % [kg/m3]
    % mat{2}.c        = 6.86  % [m/s]
    % %------------------------------------------------------------------------

    % fluid properties (air)
    %--------------------------------------------------------------------------
    mat{1}.rhof     = 1.02  % [kg/m3]
    mat{1}.c        = 370.48   % [m/s]
    %-------------------------------------------------------------------------
    
    
    % fluid properties (water)
    %--------------------------------------------------------------------------
    mat{2}.rhof     = 1000   % [kg/m3]
    mat{2}.c        = 1414.2 % [m/s]
    %------------------------------------------------------------------------
end


% define tensor basis
%--------------------------------------------------------------------------
b  = {'e1'; 'e2'};
ee = cartesianbasis2d(  b{1},   b{2});
e1 = ee(1);
e2 = ee(2);
%--------------------------------------------------------------------------


% calculate fundamental tensors
%--------------------------------------------------------------------------
  I   = identity(2, ee);
  I4  = identity(4, ee);
  I4S = 1/2 * (  I4 + rtranspose( I4));
%--------------------------------------------------------------------------


%% SORTING FLUID ELEMENTS

% FLUID
%--------------------------------------------------------------------------
    % fluid elements
      fluid_elems = 1:m;
    % conn fluid
      connf       =   conn(  fluid_elems,:);
    % fluid nodes
      nodes_f     = unique(reshape(  connf,1,[]));
    % number of fluid elements in the mesh
      mf          = size(  connf,1);
    % number of fluid nodes
      nf          = length(  nodes_f);
%--------------------------------------------------------------------------

if multiple_fluid_phase

    
    % FLUID 1
    %---------------------------------------------------------------------- 
        % solid elements
          fluid1_elems =[];
          for i=list_fluid1_tag
          fluid1_elems = [fluid1_elems find( mesh.elem{2,2}==i).'];
          end
        % connectivity of the solid nodes
          connff{1}       =   conn(  fluid1_elems,:);
        % solid nodes
          nodes_ff{1}     = unique(reshape(  connff{1},1,[]));
        % number of solid elements in the mesh
          mff{1}          = size  (  connff{1}, 1);
        % number of solid nodes
          nff{1}          = length(  nodes_ff{1});
    %----------------------------------------------------------------------
    
    % FLUID 2
    %----------------------------------------------------------------------
        % solid elements
          fluid2_elems =[];
          for i=list_fluid2_tag
          fluid2_elems = [fluid2_elems find( mesh.elem{2,2}==i).'];
          end
        % pconnectivity of the solid nodes
          connff{2}       =   conn(  fluid2_elems,:);
        % solid nodes
          nodes_ff{2}     = unique(reshape(  connff{2},1,[]));
        % number of solid elements in the mesh
          mff{2}          = size  (  connff{2}, 1);
        % number of solid nodes
          nff{2}          = length(  nodes_ff{2});
    %----------------------------------------------------------------------

end


%% PLOT MESH
 
if multiple_fluid_phase
    
    %--------------------------------------------------------------------------
    % plot each fluid phase
    figure(2)
    clf
    daspect([1 1 1]);
    hold on;
    plotMesh(  mx,  connff{1},  [0 0.5 1] );  %blue
    plotMesh(  mx,  connff{2},  [.7 .7 .7] ); %gray
    hold off;
%     plot_edges(conni(:,[1 3]),x,ee)
    % exportgraphics(gca,'plot.png','BackgroundColor','none')
    xlabel('m');
    ylabel('m');
    %--------------------------------------------------------------------------

else
    % good looking mesh plot
    %--------------------------------------------------------------------------
    % plot each solid phase and fluid phase
    figure(2)
    clf
    daspect([1 1 1]);
    hold on;
    % plotMesh(  mx,  connf, [0 0.5 1] );   %blue
    plotMesh(  mx,  connf, [1 1 1] );   %blue
    hold off;
    xlabel('m');
    ylabel('m');
    %--------------------------------------------------------------------------
end

% % mesh where you can see the node numbering
% %--------------------------------------------------------------------------
% %plot each phase
% figure(2)
% clf
% daspect([1 1 1]);
% hold on;
% % femplot(  x,  conns,'Color' , 'black', 'Nodes', 'off', ...
% % 'NodeNumbers' , 'off', 'ElementNumbers', 'off');
% femplot(  x,  connf,'Color' , 'blue', 'Nodes', 'off', ...
% 'NodeNumbers' , 'off', 'ElementNumbers', 'off');
% hold off;
% % exportgraphics(gca,'plot.png','BackgroundColor','none')
% %--------------------------------------------------------------------------




%% BOUNDARY NODES SORTING
disp('codeblock: NODE SORTING')

precision = 1e-8; % mesh 

% nodes of the fluid boundaries
%--------------------------------------------------------------------------
% nodes on the left boundary of the mesh
[  l,~]        = find(abs(  mx(:,1) +   ax/2 )<precision );
  left_f       = intersect(  l,  nodes_f)';  
%   left_s       = intersect(  l,  nodes_s)';

% nodes on the right boundary of the mesh
[  r,~]        = find(abs(  mx(:,1) -   ax/2 )<precision );
  right_f      = intersect(  r,  nodes_f)';  
%   right_s      = intersect(  r,  nodes_s)';

% nodes on the bottom boundary of the mesh
[  bottom_f,~] = find(abs(  mx(:,2) +   ay/2 )<precision );
   bottom_f    = bottom_f.';

% nodes on the bottom boundary of the mesh
[  top_f,~]    = find(abs(  mx(:,2) -   ay/2 )<precision );
   top_f       = top_f.';

% corner nodes
  corners_f = [   intersect(  l,  bottom_f) ...
                  intersect(  r,  bottom_f) ...
                  intersect(  r,  top_f   ) ...
                  intersect(  l,  top_f   )];

% eliminate corners from the edges  
left_f       = setdiff(left_f,corners_f);
right_f      = setdiff(right_f,corners_f);
top_f      = setdiff(top_f,corners_f);
bottom_f      = setdiff(bottom_f,corners_f);


clear   l; clear   r;
%--------------------------------------------------------------------------



%% VOLUME FRACTIONS AND POROSITY COMPUTATION

if not(multiple_fluid_phase)

% VOLUME OF FLUID ELEMENTS 
%--------------------------------------------------------------------------

% initialize total volume of solid phase
Vf = 0;

    for e = 1:m % loop over all fluid elements
     
    % extract nodes
    iie =   connf(e, :); % nodes of the current element e
    xe  =   x(iie);            % coordinates of these nodes

    % get element area
    Vf=Vf+get_element_area(xe,ee);

    end % end of element loop


% compute porosity
phi= Vf/(ax*ay);
%--------------------------------------------------------------------------

end




%%  4 INTEGRATION POINTS WITHIN MASTER ELEMENT
%--------------------------------------------------------------------------
% Notation from TensorLab
% quadratic quadrilateral element:
%         4----7----3
%         |         |
%         8         6
%         |         |
%         1----5----2
%      then CONN = [1 2 3 4 5 6 7 8].
%
% The coordinates of integration points x (in the coordinate system of the
% master element). Exact integration up to polynomial of order 5.
%          --------- 
%         | x  x  x |
%         | x  x  x |
%         | x  x  x |
%          ---------
% and the weight factors
%--------------------------------------------------------------------------
xi = [ -sqrt(3/5) *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 + sqrt(3/5) *e2
       -sqrt(3/5) *e1 + sqrt(3/5) *e2
                0 *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 +         0 *e2       
                0 *e1 + sqrt(3/5) *e2                
       -sqrt(3/5) *e1 +         0 *e2
                0 *e1 +         0 *e2   ];
w  = [ 25/81
       25/81
       25/81
       25/81
       40/81
       40/81
       40/81
       40/81
       64/81     ];
%--------------------------------------------------------------------------


%% ASSEMBLE STIFFNESS AND MASS MATRICES IN A FLUID ELEMENT LOOP 
disp('codeblock: ASSEMBLING FLUID MASS AND STIFFNESS MATRICES')
%--------------------------------------------------------------------------
% following the notation in my notes
% governing wave equation's unit is [1/s^2], like in COMSOL
% A - [1/Pa]      equivalent "mass matrix", term that multiplies d^2/dt^2 p 
% B - [1/(kg/m3)] equivalent "stiffness matrix", term that multiplies  p 
% d - [m/s2]      RHS term "external forces"
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% A <--> Q (paper notation)
% B <--> H (paper notation)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
Ann = sparse(  nf,   nf); % zero matrix of nxn scalars 
Bnn = sparse(  nf,   nf); % zero matrix of nxn scalars
%--------------------------------------------------------------------------



if not(multiple_fluid_phase)
%--------------------------------------------------------------------------

   for e = 1:size(  connf,1)% loop over all fluid elements
    
    % display computation percentage
    if mod(floor(e/mf*100),10)==0
    fprintf('fluid: assembled %.2f ', e/mf*100); disp('%');
    end

    iie =   connf(e, :); % nodes of the current element e
    xe  =   x(iie);     % coordinates of these nodes
        
    % compute element matrices in integration point loop       
    Ae = zeros(8, 8);        % zero matrix of 8x8 scalars     
    Be = zeros(8, 8);        % zero matrix of 8x8 scalars

    for k = 1:length(w) % loop over 9 integration points 
        
        xi1 = dot(xi(k), e1); % x1 coordinate of the current integr. point
        xi2 = dot(xi(k), e2); % x2 coordinate of the current integr. point

       % column of the shape functions
        Ne = [ -1/4*(1-xi1  )*(1-xi2  )*(1+xi1+xi2)
               -1/4*(1+xi1  )*(1-xi2  )*(1-xi1+xi2)
               -1/4*(1+xi1  )*(1+xi2  )*(1-xi1-xi2)
               -1/4*(1-xi1  )*(1+xi2  )*(1+xi1-xi2)
                1/2*(1-xi1^2)*(1-xi2  )
                1/2*(1+xi1  )*(1-xi2^2)
                1/2*(1-xi1^2)*(1+xi2  )
                1/2*(1-xi1  )*(1-xi2^2)           ];
        
        % column of the gradient of the shape functions
        % with respect to the local coordinates of the mater element
        gradxiNe = [ - 1/4 * (-1 + xi2   )* (2*xi1 +   xi2) *e1 ...
                     - 1/4 * (-1 + xi1   )* (xi1   + 2*xi2) *e2
                     - 1/4 * (2*xi1 - xi2)* (-1    + xi2  ) *e1 ...
                     - 1/4 * (1 + xi1    )* (xi1   - 2*xi2) *e2
                       1/4 * (1 + xi2    )* (2*xi1 + xi2  ) *e1 ...
                     + 1/4 * (1 + xi1    )* (xi1   + 2*xi2) *e2
                       1/4 * (2*xi1 - xi2)* (1     + xi2  ) *e1 ...
                     + 1/4 * (-1 + xi1   )* (xi1 - 2*xi2  ) *e2
                                      xi1 * (-1 + xi2     ) *e1 ...
                                    + 1/2 * (-1 + xi1^2   ) *e2
                                      1/2 * ( 1 - xi2^2   ) *e1 ...
                                    - xi2 * ( 1 + xi1     ) *e2
                                    - xi1 * ( 1 + xi2     ) *e1 ...
                                    + 1/2 * ( 1 - xi1^2   ) *e2
                                      1/2 * (-1 + xi2^2   ) *e1 ...
                                    + xi2 * (-1 + xi1     ) *e2 ];
        % Jacobian
        J = gradxiNe' * xe;
       
        % column of the gradient of the shape functions
        % with respect to the global coordinates of the mesh
        gradNe = dot(inv(J), gradxiNe);
        % element matrix and right hand side
        Ae = Ae + w(k)*(1/  rhof)*Ne*(1/  c^2)*Ne'*  thickout * det(J);
        Be = Be + w(k)*(1/  rhof)*dot(gradNe,gradNe')*  thickout*det(J);     
   
    end % end of interation point loop
   
    % assembly
    Ann(iie, iie) = Ann(iie, iie) + Ae;
    Bnn(iie, iie) = Bnn(iie, iie) + Be;
   end % end of fluid element loop

% define actual stiffness and mass matrices.
A=Ann(  nodes_f,  nodes_f);
B=Bnn(  nodes_f,  nodes_f);

% force symmetrization (correct numerical errors)
A = (A+A')/2;
B = (B+B')/2;

end



if multiple_fluid_phase
    %--------------------------------------------------------------------------
    for i =1:n_fluid_phases % loop on fluid phases
        
           % material parameters of phase i
            rhof = mat{i}.rhof
            c    = mat{i}.c   
            

            for e = 1:mff{i}% loop over all fluid elements
            
            % display computation percentage
            if mod(floor(e/mff{1}*100),10)==0
            fprintf('fluid %d: assembled %.2f ', i, e/mff{i}*100); disp('%');
            end
        
            iie =   connff{i}(e, :); % nodes of the current element e
            xe  =   x(iie);     % coordinates of these nodes
                
            % compute element matrices in integration point loop       
            Ae = zeros(8, 8);        % zero matrix of 8x8 scalars     
            Be = zeros(8, 8);        % zero matrix of 8x8 scalars
        
            for k = 1:length(w) % loop over 9 integration points 
                
                xi1 = dot(xi(k), e1); % x1 coordinate of the current integr. point
                xi2 = dot(xi(k), e2); % x2 coordinate of the current integr. point
        
               % column of the shape functions
                Ne = [ -1/4*(1-xi1  )*(1-xi2  )*(1+xi1+xi2)
                       -1/4*(1+xi1  )*(1-xi2  )*(1-xi1+xi2)
                       -1/4*(1+xi1  )*(1+xi2  )*(1-xi1-xi2)
                       -1/4*(1-xi1  )*(1+xi2  )*(1+xi1-xi2)
                        1/2*(1-xi1^2)*(1-xi2  )
                        1/2*(1+xi1  )*(1-xi2^2)
                        1/2*(1-xi1^2)*(1+xi2  )
                        1/2*(1-xi1  )*(1-xi2^2)           ];
                
                % column of the gradient of the shape functions
                % with respect to the local coordinates of the mater element
                gradxiNe = [ - 1/4 * (-1 + xi2   )* (2*xi1 +   xi2) *e1 ...
                             - 1/4 * (-1 + xi1   )* (xi1   + 2*xi2) *e2
                             - 1/4 * (2*xi1 - xi2)* (-1    + xi2  ) *e1 ...
                             - 1/4 * (1 + xi1    )* (xi1   - 2*xi2) *e2
                               1/4 * (1 + xi2    )* (2*xi1 + xi2  ) *e1 ...
                             + 1/4 * (1 + xi1    )* (xi1   + 2*xi2) *e2
                               1/4 * (2*xi1 - xi2)* (1     + xi2  ) *e1 ...
                             + 1/4 * (-1 + xi1   )* (xi1 - 2*xi2  ) *e2
                                              xi1 * (-1 + xi2     ) *e1 ...
                                            + 1/2 * (-1 + xi1^2   ) *e2
                                              1/2 * ( 1 - xi2^2   ) *e1 ...
                                            - xi2 * ( 1 + xi1     ) *e2
                                            - xi1 * ( 1 + xi2     ) *e1 ...
                                            + 1/2 * ( 1 - xi1^2   ) *e2
                                              1/2 * (-1 + xi2^2   ) *e1 ...
                                            + xi2 * (-1 + xi1     ) *e2 ];
                % Jacobian
                J = gradxiNe' * xe;
               
                % column of the gradient of the shape functions
                % with respect to the global coordinates of the mesh
                gradNe = dot(inv(J), gradxiNe);
                % element matrix and right hand side
                Ae = Ae + w(k)*(1/  rhof)*Ne*(1/  c^2)*Ne'*  thickout * det(J);
                Be = Be + w(k)*(1/  rhof)*dot(gradNe,gradNe')*  thickout*det(J);     
           
            end % end of interation point loop
           
            % assembly
            Ann(iie, iie) = Ann(iie, iie) + Ae;
            Bnn(iie, iie) = Bnn(iie, iie) + Be;
           end % end of fluid element loop
        
        % define actual stiffness and mass matrices.
        A=Ann(  nodes_f,  nodes_f);
        B=Bnn(  nodes_f,  nodes_f);
        
        % force symmetrization (correct numerical errors)
        A = (A+A')/2;
        B = (B+B')/2;
    %----------------------------------------------------------------------
    end % of fluid-phase loop
end % if 


% clear big matrices
%--------------------------------------------------------------------------
clear Ae; clear Be;
clear Bnn; clear Ann;
%--------------------------------------------------------------------------


%% INDEX OF DISPLACEMENT/PRESSURE DOF UNDER NODES_S/NODES_F ORDERING
%--------------------------------------------------------------------------
% useless codeblock inherited from FSI problem:
% e.g vector 'doff_l' is the same as 'left_f' because there is no solid.  
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% fluid
    % left nodes
    doff_l = get_ind(  nodes_f,   left_f   ); 
    % right nodes
    doff_r = get_ind(  nodes_f,   right_f  ); 

    % bottom nodes
    doff_b = get_ind(  nodes_f,   bottom_f ); 
    % top nodes
    doff_t = get_ind(  nodes_f,   top_f    ); 
   
    % corner nodes
    doff_c = get_ind(  nodes_f,   corners_f);
%--------------------------------------------------------------------------





%% PRESCRIBED AND FREE NODES SPLIT
%--------------------------------------------------------------------------
% INDEX OF DISPLACEMENT DOF IN NODES_S (ELIMINATED DEPENDENT DOFS)
% same index change procedure must be executed but only in the solid phase.

%|__NOTES_NOTATION_|____ DESCRIPTION_______________________________________
%|                 |
%|        p'       | prescribed nodes in the fluid phase
%|_________________|_______________________________________________________
%|                 |
%|        f'       | free nodes in the fluid phase
%|_________________|_______________________________________________________
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% solid
p=[];
f=[];

% fluid
pp = [doff_c doff_l doff_r doff_b doff_t];
fp = setdiff(1:  nf,pp);
%--------------------------------------------------------------------------


%% PARTITIONING
%--------------------------------------------------------------------------
% A <--> Q (paper notation)
% B <--> H (paper notation)
%--------------------------------------------------------------------------

% fluid
    % 'mass' scalar matrix 
    % B =  Bp'p' Bp'f'   
    %      Bf'p' Bf'f' 
    A_pp_pp = A (pp, pp); A_pp_fp = A (pp, fp);
    A_fp_pp = A (fp, pp); A_fp_fp = A (fp, fp);

    % 'stiffness' scalar matrix 
    % B =  Bp'p' Bp'f'   
    %      Bf'p' Bf'f' 
    B_pp_pp = B (pp, pp); B_pp_fp = B (pp, fp);
    B_fp_pp = B (fp, pp); B_fp_fp = B (fp, fp);

%--------------------------------------------------------------------------




%% ASSEMBLE NON-SYMMETRIC SYSTEM OF EQUATIONS 
disp('codeblock: SOLVE')

% renaming (same notation as IJNME paper Transient computational 
% homogenization of heterogeneous poroelastic media with local resonances")


% lamda submatrices
%--------------------------------------------------------------------------
lam_P_P = [  A_pp_pp ];

lam_P_F = [ A_pp_fp ];

lam_F_P = [A_fp_pp ];

lam_F_F = [A_fp_fp ];
%--------------------------------------------------------------------------


% mu submatrices
%--------------------------------------------------------------------------
mu_P_P  = [  B_pp_pp    ];  

mu_P_F  = [B_pp_fp    ];  

mu_F_P  = [ B_fp_pp    ];  

mu_F_F  = [ B_fp_fp    ];  
%--------------------------------------------------------------------------


%% EIGEN: EIGEN STUDY, PRINT EIGENFREQUENCIES AND EIGENMODES


% for only-fluid problem the matrices are symmetric with a certain
% numerical precision.
%--------------------------------------------------------------------------
    disp('ABSOLUTE SYMMETRIZATION FORCED')
    mu_F_F   = (mu_F_F + mu_F_F.')/2 ;
    lam_F_F  = (lam_F_F+ lam_F_F.')/2;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% define desired number of modes to compute and display
n_modes = 5;
disp('ITERATIVE SOLVER')
    % Subset of eigenvalues -> CHEAPER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % eigs is recommended by matlab when sparse nonsymmetric matrices.
    % % compute eigen modes V, and eigen values D with smallest absolute
    % value right problem.
    [phi_F_Q, Dr] = eigs(mu_F_F ,lam_F_F ,n_modes,...
                   'smallestabs',...
                   'Display',0,...
                   'FailureTreatment','replacenan');
    % % left problem (psi appears always transposed along the expressions)
    % [psi_F_Q, Dl] = eigs(mu_F_F',lam_F_F',n_modes,...
    %                'smallestabs',...
    %                'Display',0,...
    %                'FailureTreatment','replacenan');
    
    % CHECK whether the frequencies in Dr and Dl are enoughly equal. For
    % the current example they are equal up to the 5th digit after comma.
    

    % print frequencies f
    freqs = diag(sqrt(Dr)/(2*pi))
    fprintf('Eigen frequencies:\n');
    fprintf('%f Hz\n', freqs);
%--------------------------------------------------------------------------


%direct solver
%--------------------------------------------------------------------------
% disp('DIRECT SOLVER')
% [phi_F_Q, Dr] = eig(full(mu_F_F) ,full(lam_F_F));
% % print frequencies f
% freqs = sort(diag(sqrt(Dr)/(2*pi)));
% fprintf('Eigen frequencies:\n');
% fprintf('%f Hz\n', freqs(1:n_modes));
%--------------------------------------------------------------------------
    

% EIGEN:  PHASE CORRECTION in case left eigenvector is diff than right one
%--------------------------------------------------------------------------
% IMPORTANT NOTE: PHASE CORRECTION BEFORE NORMALIZATION 
% Compute round(psi'*lamFF*phi). It should be positive diag matrix - note 
% that this product should be approx positive diag because the the
% frequencies must be positive-. Since the right and left eigen modes were
% computed with different calls from the function eigs, they can
% potentially converge to a shape mode with phase difference of pi. In 
% this case, psi'*lamFF*phi ( or round(psi'*lamFF*phi) for getting rid 
% of almost zero off diag elements) may have some of the diagonal elements
% as -negative instead of +positive. To correct this effect it is necessary
% to multiply either phi or psi by the 'sign' matrix 'I'= psi'*lamFF*phi.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% get diagonal
diag_psi_lam_phi   = diag(psi_F_Q'*lam_F_F*phi_F_Q);
% make it a vector with only the sign of each element
diag_psi_lam_phi   = diag_psi_lam_phi./abs(diag_psi_lam_phi);
% make 'identity' matrix with sign correction
I_phase_correction = diag(diag_psi_lam_phi);

% correct the sign of left eigen modes with respect to the right modes
psi_F_Q = psi_F_Q * I_phase_correction;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% save version of phi for plotting
phi_plot = phi_F_Q;
%--------------------------------------------------------------------------

% EIGEN:  MODE NORMALIZATION WITH RESPESCT TO THE LAMDA MATRIX
%--------------------------------------------------------------------------
% ready to mass normalization, now, the product psi'*lamFF*phi has only
% postive diagonal compared to previous case
% get norm of each mode
vec_norms = sqrt(diag(psi_F_Q'*lam_F_F*phi_F_Q));
phi_F_Q=phi_F_Q./vec_norms'; %divide each column of phi by element of vec_norms 
psi_F_Q=psi_F_Q./vec_norms'; %divide each column of phi by element of vec_norms 
% CONSISTENCY CHECK: 
% Compute again psi'*lamFF*phi. It should be approx. the identity matrix.
% when eigs is used the phase problem doen't exist. It can be seen by
% computing psi(:,1:n_modes)'*lamFF*phi(:,1:n_modes). The diagonal is made
% out of positive integers.
%--------------------------------------------------------------------------





% EIGEN: ASSIGN DISP. SHAPE MODES TO TENSOR VECTOR & PLOT
%--------------------------------------------------------------------------
% displacement 'ured' is followed is 'u' displacement reduced 'red' thanks
% to the elimitation of the dependend points due to periodic bc
% 'ured' stands for 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% prescribed (reduced) displacements are zero
    % declare & assign
ured_p = zeros(1, ee, length(p), size(phi_F_Q,2) ); 
        
% free (reduced) displacements are the eigenvectors    
ured_f = phi_plot(1:2:2*length(f)-1,:)*e1+ phi_plot(2:2:2*length(f),:)*e2;

% total reduced displacements (indices ordered as [dofs_pc dofs_un dofs_in])
    % declare ured(unconstrained+corners+independent nodes, modes)
    ured = zeros(1, ee, length(p)+length(f), size(phi_F_Q,2));
    % assign
    % assign prescribed nodes
    ured(p,:) = ured_p;
    % assign free nodes
    ured(f,:) = ured_f;
    % free unused variables
    clear ured_f; clear uref_p;
%--------------------------------------------------------------------------



% ASSIGN PRESSURE SHAPE MODES TO SCALAR VECTOR vp
%--------------------------------------------------------------------------
% extract pressure values from bc and eigen vectors
    % prescribed pressures (same for every mode)    
    p_pp     = zeros(length(pp), size(phi_F_Q,2)); 
    
    % free pressures (different for each mode)
    p_fp = phi_plot(end-length(fp)+1:end,:);

% declare pressure scalar vector
vp = zeros(  nf,size(phi_F_Q,2));

% assign to scalar vector solution of size
vp(pp,:)=p_pp;
vp(fp,:)=p_fp;

% prepare for plot
    % declare pressure scalar vector for plot
    vp_plot = zeros(  n,size(phi_F_Q,2));
    
    % assign vp to vector for plot
    vp_plot(  nodes_f,:) = vp;
%--------------------------------------------------------------------------




%% EIGEN: PLOT PRESSURE MODE SHAPE (RESCALED DISP.)
%--------------------------------------------------------------------------
disp('codeblock: PLOT')



sfreqs=freqs; %sorted frequencies vector
% mode=3;

for mode=1:n_modes

% plot the pressure
%--------------------------------------------------------------------------
figure
clf
sign=-1; % control modes shape phase ( 0 or 180 degrees)
scatter(dot( x(nodes_f),e1),dot(x(nodes_f),e2),50,sign*vp(:,mode),'filled')
colorbar;
pmax = max(abs(vp_plot(:,mode)));
caxis([-pmax pmax]);
colormap('jet');
axis equal
title(['Mode ' num2str(mode) ', frequency = ' num2str(sfreqs(mode)) ' Hz']);
hold off
%--------------------------------------------------------------------------
end

%% STEADY(LONG WAVELENGTH): PRESSURE AND ITS GRADIENTS

% reference position vector
%--------------------------------------------------------------------------
% xR = (    x(  corners_f(1)) +   x(  corners_f(3))   )/2;
xR = 0*e1;

% ones
ones_p = ones(length(p) ,1);
ones_pp = ones(length(pp),1);
%--------------------------------------------------------------------------




% precribed fluid points are [left_f right_f]
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL INPUT (JUST FOR VISUALIZATION OF STEADY MODE)
    pM     = 1;
    % gradpM =  1e5*(e1+2*e2);
    gradpM =  0*1e5*e2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %delta xf_pp
    dxf_pp = x([corners_f left_f right_f bottom_f top_f])-xR*ones_pp;
    % first-order comp. homogenization
    p_pp   =  pM *ones(length(pp),1) + ...
              dot(gradpM,  dxf_pp );
if length(pp)~=length([corners_f left_f right_f bottom_f top_f])
    error('something wrong');end

% change name pressure node vector
wP = p_pp;
%--------------------------------------------------------------------------



% %% STATIONARY: SOLVE RIGHT/LEFT LINEAR SYSTEM OF EQUATIONS

%--------------------------------------------------------------------------
% right constrain modes
S_F_P = - inv(mu_F_F  ) * mu_F_P  ;
% right stationary hybrid state vector
wr_f  =  S_F_P * wP;

% left constrain modes (Y appears always transposed along the expressions)
Y_F_P = - inv(mu_F_F.') * mu_P_F.';
% left stationary hybrid state vector
wl_f  =  Y_F_P * wP;
%--------------------------------------------------------------------------



% %% STATIONARY: ASSIGN PRESSURE SHAPE MODES TO SCALAR VECTOR vp
%--------------------------------------------------------------------------
% extract pressure values from state vector wr_f
    % free pressures 
    p_fp = wr_f(end-length(fp)+1:end);

% declare pressure scalar vector
vp_steady = zeros(  nf,1);

% assign to scalar vector solution of size
vp_steady(pp,:)=p_pp;
vp_steady(fp,:)=p_fp;

% prepare for plot length(vp_plot)>length(vp)
    % declare pressure scalar vector for plot - whole mesh vector
    vp_plot = zeros(  n,1);
    
    % assign vp to vector for plot
    vp_plot(  nodes_f) = vp_steady;
%--------------------------------------------------------------------------

disp('codeblock: PLOT')

% %% STATIONARY: PLOT PRESSURE/DISP SHAPE (RESCALED DISP.)




% plot the pressure solution 
%--------------------------------------------------------------------------
    figure(5)
    clf
    % scatter(dot(  x,e1),dot(  x,e2),30,vp_plot,'filled')
    scatter(dot( x(nodes_f),e1),dot(x(nodes_f),e2),200,vp_steady,'filled')
    colorbar;
    pmax = max(abs(vp_plot));
    caxis([-pmax pmax]);
    colormap('jet');
    axis equal
    title("Steady" );
    hold off
%--------------------------------------------------------------------------


%% REDUCED COUPLED DYNAMIC MODEL
% Comments
%--------------------------------------------------------------------------
% 1. check about transposition .' or hermitian transposition '
% regular transposition is being used.
% 2. notation is, for instance, tlam_P_P means matrix with second order
% tensor components of size = (P,P).
% t - matrix with second order tensor components
% v - matrix with first order (vector) tensor components 
% s - matrix with scalar components
%--------------------------------------------------------------------------


% compute reduced matrices on P nodes with Q eigenmodes
%--------------------------------------------------------------------------
% steady part
lam_qs  =  full( ...
                 lam_P_P   +  lam_P_F * S_F_P   +   Y_F_P.' * lam_F_P  ...
                + Y_F_P.'  *  lam_F_F * S_F_P ...
                                                                         );

mu_qs   =  full( ...
                 mu_P_P   +   mu_P_F * S_F_P   +   Y_F_P.' *  mu_F_P  ...
                + Y_F_P.' *   mu_F_F * S_F_P  ...
                                                                         );

lam_Q_P =   psi_F_Q.' * lam_F_P  +   psi_F_Q.' * lam_F_F * S_F_P  ;

lam_P_Q =   lam_P_F   * phi_F_Q  +     Y_F_P.' * lam_F_F * phi_F_Q;



% dynamics part
I_Q_Q    =    eye(n_modes);
LAM_Q_Q  =    psi_F_Q(:,1:n_modes)'*mu_F_F*phi_F_Q(:,1:n_modes);
O_Q_1    =    zeros(n_modes,1);
O_Q_2    =    zeros(n_modes,2);
%--------------------------------------------------------------------------

% convert equivalent scalar matrix back to tensor matrix form
%--------------------------------------------------------------------------
    % lam_qs hybrid matrix partitioning
    % lam_qs = tlam_p_p   vlam_p_p'
    %          vlam_p'_p  slam_p'_p' 
    tlam_p_p   = sm2tm(lam_qs(1:dim*length(p),1:dim*length(p)),ee,2);
    vlam_p_pp  = sm2tm(lam_qs(1:dim*length(p),dim*length(p)+1:end),ee,1);
    vlam_pp_p  = sm2tm(lam_qs(dim*length(p)+1:end,1:dim*length(p)).',ee,1).';
    slam_pp_pp =       lam_qs(dim*length(p)+1:end,dim*length(p)+1:end);

    % mu_qs hybrid matrix partitioning
    % mu_qs =  tmu_p_p   vmu_p_p'
    %          vmu_p'_p  mu_p'_p'
    % note vmu_p'_p is identically zero from the theory
    tmu_p_p   = sm2tm(mu_qs(1:dim*length(p),1:dim*length(p)),ee,2);
    vmu_p_pp  = sm2tm(mu_qs(1:dim*length(p),dim*length(p)+1:end),ee,1);
    vmu_pp_p  = sm2tm(mu_qs(dim*length(p)+1:end,1:dim*length(p)).',ee,1).';
    smu_pp_pp =       mu_qs(dim*length(p)+1:end,dim*length(p)+1:end);

    % lam_P_Q hybrid matrix partitioning
    % lam_P_Q = vm_p_Q   
    %           sm_p'_Q 
    vm_p_Q   = sm2tm(lam_P_Q(1:dim*length(p)    ,:),ee,1);
    sm_pp_Q  =       lam_P_Q(dim*length(p)+1:end,:);
    
    % lam_Q_P hybrid matrix partitioning
    % lam_Q_P = [vm_Q_p  sm_Q_p'] 
    % transpose because contraction happens on the column
    vm_Q_p   = sm2tm(lam_Q_P(: , 1:dim*length(p)).',ee,1).';
    sm_Q_pp  =       lam_Q_P(: , dim*length(p)+1:end  );
%--------------------------------------------------------------------------



% compute homogenized macroscopic fluid acceleration
%--------------------------------------------------------------------------
% note: h stands for 'homogenized' and m for 'momentum'
hvA = -1/(ax*ay)*                dxf_pp.' * vlam_pp_p  * ones_p  ;
% hvB = -1/(ax*ay)*                dxf_pp.' * vlam_pp_p  * dxs_p   ;
% hvB = -1/(ax*ay)*                0;
hvC = -1/(ax*ay)*                dxf_pp.' * slam_pp_pp * ones_pp ;
hvD = -1/(ax*ay)*                dxf_pp.' * slam_pp_pp * dxf_pp  ;
hvE = -1/(ax*ay)*                dxf_pp.' * vmu_pp_p   * ones_p  ;
% hvF = -1/(ax*ay)*                dxf_pp.' * vmu_pp_p   * dxs_p   ;
hvG = -1/(ax*ay)*                dxf_pp.' * smu_pp_pp  * ones_pp ;
hvH = -1/(ax*ay)*                dxf_pp.' * smu_pp_pp  * dxf_pp  ;
hvL = -1/(ax*ay)*               (sm_pp_Q.'* dxf_pp)             ; %corrected 20/04/2023-
% for the evolution equation
hvLs= -1/(ax*ay)*               (sm_Q_pp  * dxf_pp)             ; %corrected 04/05/2023
hvl =                          -(sm_Q_pp  * dxf_pp)             ; %corrected 04/05/2023

% full contraction of each tensor for an estimation of order
abshvA_BIOT = sqrt(    ddot(hvA,hvA) );
% abshvB = sqrt(   dddot(haB,haB) );
abshvC = sqrt(     dot(hvC,hvC) );
abshvD = sqrt(    ddot(hvD,hvD) );
abshvE = sqrt(    ddot(hvE,hvE) );
% abshvF = sqrt(   dddot(haF,haF) );
abshvG = sqrt(     dot(hvG,hvG) );
abshvH_BIOT = sqrt(    ddot(hvH,hvH) );

abshvL = [];
for i=1:size(hvL,1)
abshvL =  [abshvL  sqrt( dot(hvL(i),hvL(i)))];
end
%--------------------------------------------------------------------------


% compute homogenized macroscopic fluid volumetric rate rate
%--------------------------------------------------------------------------
% note: h stands for 'homogenized' and m for 'momentum'
heA = -1/(ax*ay)*                ones_pp.' * vlam_pp_p  * ones_p  ;
% heB = -1/(ax*ay)*                ones_pp.' * vlam_pp_p  * dxs_p   ;
heC = -1/(ax*ay)*                ones_pp.' * slam_pp_pp * ones_pp ;
heD = -1/(ax*ay)*                ones_pp.' * slam_pp_pp * dxf_pp  ;
heE = -1/(ax*ay)*                ones_pp.' * vmu_pp_p   * ones_p  ;
% heF = -1/(ax*ay)*                ones_pp.' * vmu_pp_p   * dxs_p   ;
heG = -1/(ax*ay)*                ones_pp.' * smu_pp_pp  * ones_pp ;
heH = -1/(ax*ay)*                ones_pp.' * smu_pp_pp  * dxf_pp  ;
heL = -1/(ax*ay)*               (ones_pp.' * sm_pp_Q).'           ;%corrected 20/04/2023-
% for the evolution equation
heLs = -1/(ax*ay)*              (sm_Q_pp   * ones_pp)             ;%corrected 04/05/2023
hel  =                         -(sm_Q_pp   * ones_pp)             ;%corrected 04/05/2023


% full contraction of each tensor for an estimation of order
absheA = sqrt(     dot(heA,heA) );
% absheB_BIOT = sqrt(    ddot(heB,heB) );
absheC_BIOT = sqrt(         heC*heC  );
absheD = sqrt(     dot(heD,heD) );
absheE = sqrt(     dot(heE,heE) );
% absheF = sqrt(    ddot(heF,heF) );
absheG = sqrt(         heG*heG  );
absheH = sqrt(     dot(heH,heH) );

absheL = [];
for i=1:size(heL,1)
absheL =  [absheL   sqrt( heL(i)*heL(i) )];
end
%--------------------------------------------------------------------------
%% PLOT COEFFICIENTS' MAGNITUDE VALUE
clear etaticklabels

% build the x tick labels
% mode labels
etaticklabels = cell(1,n_modes);
for i=1:n_modes
    etaticklabels{i} =  ['(' num2str(round(freqs(i),0)) ...
                         ' Hz)$ \eta_{' num2str(i) '} $'];
end     

% make artificial x axis for plot
xx=1:8+n_modes;

figure(6)
clf
set(groot,'defaultAxesTickLabelInterpreter','latex');  


% FLUID VELOCITY
%--------------------------------------------------------------------------
subplot(1,2,1)
Yv = [0 0 abshvC abshvD 0 0 abshvG ...
      abshvH_BIOT abshvL];
% Yv = [0 0 dimless_abshvC ...
%       dimless_abshvD 0 0 dimless_abshvG ...
%       dimless_abshvH_BIOT dimless_abshvL];
Xv = {'(BIOT) $ \ddot{\vec{u}}_M $', ...
      ' $\vec{\nabla} \ddot{\vec{u}}_{_{M}} $', ...
      ' $ \ddot{\mathrm{p}}_{_{M}} $', ...
      ' $\vec{\nabla} \ddot{\mathrm{p}}_{_{M}} $', ...
      ' $ \vec{u}_{_{M}} $', ...
      ' $\vec{\nabla} \vec{u}_{_{M}} $', ...
      ' $\mathrm{p}_{_{M}} $', ...
      ' (BIOT) $\vec{\nabla} \mathrm{p}_{_{M}} $'};
Xv = [Xv etaticklabels];
plot(xx,log10(Yv), 'kx');
hold on;
mincoeff = max(Yv(1),Yv(8)); % A and H are Biot coeffs
plot(xx([1 8]),log10(Yv([1 8])), 'cx');
plot(xx(9:8+n_modes),log10(Yv(9:8+n_modes)), 'gx');
plot(log10(mincoeff*ones(1,8)), "Color", 'c');
title( '$\dot{\vec{v}}_{_{M}} $ Fluid Velocity coefficients', 'interpreter', 'latex')
hold off;
grid on;
set(gca, 'XTick', 1:length( Xv), 'XTickLabels',  Xv);
ylabel('log10( absolute value ) ')
axi = gca;
axi.FontSize = 16; 
%--------------------------------------------------------------------------
clear Xv  mincoeff;



% FLUID VOLUMETRIC STRAIN
%--------------------------------------------------------------------------
subplot(1,2,2)

Ye = [0 0 absheC_BIOT absheD 0 0 absheG ...
      absheH absheL];
% Ye = [0 0 dimless_absheC_BIOT ...
%       dimless_absheD 0 0 ...
%       dimless_absheG dimless_absheH dimless_absheL];
Xe = {' $ \ddot{\vec{u}}_M $', ...
      ' (BIOT) $\vec{\nabla} \ddot{\vec{u}}_{_{M}} $', ...
      ' (BIOT) $ \ddot{\mathrm{p}}_{_{M}} $', ...
      ' $\vec{\nabla} \ddot{\mathrm{p}}_{_{M}} $', ...
      ' $ \vec{u}_{_{M}} $', ...
      ' $\vec{\nabla} \vec{u}_{_{M}} $', ...
      ' $\mathrm{p}_{_{M}} $', ...
      ' $\vec{\nabla} \mathrm{p}_{_{M}} $'};
Xe = [Xe etaticklabels];
plot(xx,log10(Ye), 'kx');
hold on;
mincoeff = max(Ye(2),Ye(3)); % B and C are Biot coeffs
plot(xx([2 3]),log10(Ye([2 3])), 'cx');
plot(xx(9:8+n_modes),log10(Ye(9:8+n_modes)), 'gx');
plot(log10(mincoeff*ones(1,8)), "Color", 'c');
% legend('coeff. abs', 'smallest biot');
title( '$\ddot{e}_{_{M}} $ Fluid Volumetric Strain coefficients', 'interpreter', 'latex')
hold off;
grid on;
set(gca, 'XTick', 1:length( Xe), 'XTickLabels',  Xe);
ylabel('log10( absolute value ) ')
axi = gca;
axi.FontSize = 16; 
%--------------------------------------------------------------------------
clear Xe  mincoeff;

clear etaticklabels axi;

%% DIMENSIONLESS CONSTITUTIVE EQUATIONS

% % characteristic values
% %--------------------------------------------------------------------------
% f0 = 600;
% l0 = sqrt(dot(e1,hvH,e1)/heC)/f0;% waveleng determined by homogenized par.
% % l0 = max(ax,ay);
% warning('l0 is the RVE size')
% p0 = 0.01; %check on wikipedia "sound pressure level"
% %--------------------------------------------------------------------------
% 
% 
% 
% % macroscopic fluid velocity absolute value
% %--------------------------------------------------------------------------
% % dimless_abshaA_BIOT = u0*f0^2     *abshvA_BIOT  ;
% % dimless_abshaB      = u0*f0^2/l0  *abshaB       ;
% dimless_abshvC      = p0*f0^2     *abshvC       ;
% dimless_abshvD      = p0*f0^2/l0  *abshvD       ;
% % dimless_abshvE      = u0          *abshvE       ;
% % dimless_abshvF      = u0/l0       *abshaF       ;
% dimless_abshvG      = p0          *abshvG       ;
% dimless_abshvH_BIOT = p0/l0       *abshvH_BIOT  ;
% dimless_abshvL      = f0^2        *abshvL       ;
% %--------------------------------------------------------------------------
% 
% % macroscopic fluid volumetric rate rate absolute value
% %--------------------------------------------------------------------------
% % dimless_absheA      = u0*f0^2     *absheA       ;
% % dimless_absheB_BIOT = u0*f0^2/l0  *absheB_BIOT  ;
% dimless_absheC_BIOT      = p0*f0^2     *absheC_BIOT       ;
% dimless_absheD      = p0*f0^2/l0  *absheD       ;
% % dimless_absheE      = u0          *absheE       ;
% % dimless_absheF      = u0/l0       *absheF       ;
% dimless_absheG      = p0          *absheG  ;
% dimless_absheH      = p0/l0       *absheH       ;
% dimless_absheL      = f0^2        *absheL       ;
% %--------------------------------------------------------------------------
% 
% 
% % PLOT DIMENSIONLESS COEFFICIENTS' MAGNITUDE VALUE
%  
% % mode labels
% etaticklabels = cell(1,n_modes);
% for i=1:n_modes
%     etaticklabels{i} =  ['(' num2str(round(freqs(i),0)) ...
%                          ' Hz)$ \eta_{' num2str(i) '}^* $'];
% end  
% 
% % make artificial x axis for plot
% xx=1:8+n_modes;
% 
% % without dimensions
% 
% figure(7)
% clf
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% sgtitle([ 'DIMENSIONLESS :'...
%          '$ f_0 $ = ' num2str(f0) ...
%          '[Hz], $l_0$ =' num2str(l0*100,'%.1f') ...
%          '[cm], $ p_0 $ =' num2str(p0) '[Pa]'], 'Interpreter','latex');  
% 
% % FLUID VELOCITY
% %--------------------------------------------------------------------------
% subplot(1,2,1)
% Yv = [0 0 dimless_abshvC ...
%       dimless_abshvD 0 0 dimless_abshvG ...
%       dimless_abshvH_BIOT dimless_abshvL];
% Xv = {'(BIOT) $ \ddot{\vec{u}}_M^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \ddot{\vec{u}}_{_{M}}^{*}  $', ...
%       ' $ \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
%       ' $ \vec{u}_{_{M}} ^{*} $', ...
%       ' $\vec{\nabla}^{*}  \vec{u}_{_{M}}^{*} $', ...
%       ' $\mathrm{p}_{_{M}}^{*} $', ...
%       ' (BIOT) $\vec{\nabla}^{*}  \mathrm{p}_{_{M}}^{*}  $'};
% Xv = [Xv etaticklabels];
% plot(xx,log10(Yv), 'kx');
% hold on;
% mincoeff = max(Yv(1),Yv(8)); % A and H are Biot coeffs
% plot(xx([1 8]),log10(Yv([1 8])), 'cx');
% plot(xx(9:8+n_modes),log10(Yv(9:8+n_modes)), 'gx');
% plot(log10(mincoeff*ones(1,8)), 'c');
% % legend('coeff. abs', 'smallest biot');
% title( '$\dot{\vec{v}}_{_{M}} $ Fluid Velocity coefficients', 'interpreter', 'latex')
% hold off;
% grid on;
% set(gca, 'XTick', 1:length( Xv), 'XTickLabels',  Xv);
% ylabel('log10( absolute value ) ')
% axi = gca;
% axi.FontSize = 16; 
% %--------------------------------------------------------------------------
% clear Xv Yv mincoeff;
% 
% 
% 
% % FLUID VOLUMETRIC STRAIN
% %--------------------------------------------------------------------------
% subplot(1,2,2)
% Ye = [0 0 dimless_absheC_BIOT ...
%       dimless_absheD 0 0 ...
%       dimless_absheG dimless_absheH dimless_absheL];
% Xe = {' $ \ddot{\vec{u}}_M^{*}  $', ...
%       ' (BIOT) $\vec{\nabla}^{*}  \ddot{\vec{u}}_{_{M}}^{*}  $', ...
%       ' (BIOT) $ \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
%       ' $ \vec{u}_{_{M}}^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \vec{u}_{_{M}}^{*}  $', ...
%       ' $\mathrm{p}_{_{M}}^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \mathrm{p}_{_{M}}^{*}  $'};
% Xe = [Xe etaticklabels];
% plot(xx,log10(Ye), 'kx');
% hold on;
% mincoeff = max(Ye(2),Ye(3)); % B and C are Biot coeffs
% plot(xx([2 3]),log10(Ye([2 3])), 'cx');
% plot(xx(9:8+n_modes),log10(Ye(9:8+n_modes)), 'gx');
% plot(log10(mincoeff*ones(1,8)), 'c');
% % legend('coeff. abs', 'smallest biot');
% title( '$\ddot{e}_{_{M}} $ Fluid Volumetric Strain coefficients', 'interpreter', 'latex')
% hold off;
% grid on;
% set(gca, 'XTick', 1:length( Xe), 'XTickLabels',  Xe);
% ylabel('log10( absolute value ) ')
% axi = gca;
% axi.FontSize = 16; 
% %--------------------------------------------------------------------------
% clear Xe Ye mincoeff;
% 
% clear etaticklabels axi;

%% HOMOGENIZED (CLASSIC) LONG WAVELENGTH MATERIAL PARAMETERS

% 'steady' homogenized bulk modulus (per volume of RVE)
    hKf     = - 1/ heC
   
% 'steady' homogenized density: mass per volume of RVE
    hrhof   = -inv(hvH)



% some important homogenized material parameters speed of sound
% longwavelength regime
hc_xx    = sqrt(dot(e1,hvH/heC,e1))
hc_yy    = sqrt(dot(e2,hvH/heC,e2))


%% SELECTION OF MODES

% if no mode selection
%--------------------------------------------------------------------------
n_modes_dispersion = n_modes;
%--------------------------------------------------------------------------
% 
% %     % manually select % designXu
% %     %--------------------------------------------------------------------------
% %     localized_modes = [1 4 5 10]; 
% %     %--------------------------------------------------------------------------
% %     
% %     
% %     % redefine number of modes
% %     %--------------------------------------------------------------------------
% %     n_modes_dispersion         = length(localized_modes);
% %     %--------------------------------------------------------------------------
% %     
% %     % refine some material material parameters
% %     %--------------------------------------------------------------------------
% %     hel     = hel(localized_modes);
% %     heL     = heL(localized_modes);
% %     LAM_Q_Q = LAM_Q_Q(localized_modes,localized_modes);
% %     I_Q_Q   = I_Q_Q  (localized_modes,localized_modes);
% %     O_Q_1   =  O_Q_1(localized_modes); 
% %     O_Q_2   =  O_Q_2(localized_modes);
% %     %--------------------------------------------------------------------------
% 
%     % manually select 
%     %--------------------------------------------------------------------------
%     localized_modes = [1 5]; 
%     %--------------------------------------------------------------------------
% 
% 
%     % redefine number of modes
%     %--------------------------------------------------------------------------
%     n_modes_dispersion         = length(localized_modes);
%     %--------------------------------------------------------------------------
% 
%     % refine some material material parameters
%     %--------------------------------------------------------------------------
%     hel     = hel(localized_modes);
%     heL     = heL(localized_modes);
%     LAM_Q_Q = LAM_Q_Q(localized_modes,localized_modes);
%     I_Q_Q   = I_Q_Q  (localized_modes,localized_modes);
%     O_Q_1   =  O_Q_1(localized_modes); 
%     O_Q_2   =  O_Q_2(localized_modes);
%     %--------------------------------------------------------------------------

%% DISPERSION BLOCH ANALYSIS

% figure(23)
hold on
grid on
% scatter(real(dispersion_bloch(:,1)),real(dispersion_bloch(:,2)),20);
% ylim([0 1000])
scatter(real(dispersion_bloch(:,1)),(real(dispersion_bloch(:,2))), 'Marker', '+', 'MarkerEdgeColor', 'blue', 'LineWidth', 2);
ylim([0 730])
box on

%% DISPERSION CURVE CLASSIC LONG WAVELENGTH (STEADY) TERMS


% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 1; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.001;
range_of_dimless_k = -1:kstep:1;
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda      =   zeros(number_of_waves,length(range_of_dimless_k)); % 
Gamma_index =   round(length(range_of_dimless_k)/3)+1;
X_index     = 2*round(length(range_of_dimless_k)/3)+1;
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;% index of current iteration
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else 
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end
    
    % Density matrix (k)
    
    % second line of M
    Mpp=        heC;
   
    
  

    % Elasticity matrix (k)    
    % second line of K
    Kpp=        dot(k,hvH,k)     ;
   

    % density matrix 
    M=Mpp ;
       
    % elasticity matrix 
    K=Kpp ;
    
    % compute eigevalues
%     [V,D]=eig(K,M)
    lambda(:,i)=sort(eig(K,M));       
    
i=i+1;
    end
%--------------------------------------------------------------------------


% plot
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
figure(23)
hold on
grid on
for i=1:number_of_waves
plot(range_of_dimless_k,real(frequencies(i,:)), 'k--', 'LineWidth', 2)
end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 3000])
% hold off
%--------------------------------------------------------------------------
%% DISPERSION CURVE CLASSIC LONG WAVELENGTH + elastic inertia


% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 1; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.01;
range_of_dimless_k = -1:kstep:2;
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda      =   zeros(number_of_waves,length(range_of_dimless_k)); % 
Gamma_index =   round(length(range_of_dimless_k)/3)+1;
X_index     = 2*round(length(range_of_dimless_k)/3)+1;
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;% index of current iteration
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else 
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end
    
    % Density matrix (k)
    
    % second line of M
    Mpp=        heC+dot(k,hvD,k)          ;
   
    
  

    % Elasticity matrix (k)    
    % second line of K
    Kpp=        dot(k,hvH,k)          ;
   

    % density matrix 
    M=Mpp ;
       
    % elasticity matrix 
    K=Kpp ;
    
    % compute eigevalues
%     [V,D]=eig(K,M)
    lambda(:,i)=sort(eig(K,M));       
    
    i=i+1;
    end
%--------------------------------------------------------------------------


% plot
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
figure(23)
hold on
grid on
for i=1:number_of_waves
scatter(range_of_dimless_k,real(frequencies(i,:)),'c.')
end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 2000])
% hold off
%--------------------------------------------------------------------------


%% DISPERSION CLASSIC LONG WAVELENGTH + e-MODES + v-MODES


% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 1+n_modes; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.01;
range_of_dimless_k = -1:kstep:2;
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda      =   zeros(number_of_waves,length(range_of_dimless_k)); % 
Gamma_index =   round(length(range_of_dimless_k)/3)+1;
X_index     = 2*round(length(range_of_dimless_k)/3)+1;
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;% index of current iteration
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else 
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end
    
    % Density matrix (k)
    
    % second line of M
    Mpp=       heC;
    Mpe=      (heL.'                 +1i*        dot(k,hvL.')         )       ;
    
    % third line of M
    Mep=       ( - hel                +1i*        dot( + hvl  ,k)     )     ;
    Mee= I_Q_Q                                                            ;

    % Elasticity matrix (k)    
    % second line of K
    Kpp=       dot(k,hvH,k) ;
    Kpe= O_Q_1.' ;
    
    % third line of K
    Kep= O_Q_1 ; 
    Kee= LAM_Q_Q ;                                                               

    % density matrix 
    M=[  Mpp Mpe;
         Mep Mee];
    % elasticity matrix 
    K=[Kpp Kpe;
       Kep Kee];
    
    % compute eigevalues
%     [V,D]=eig(K,M)
    lambda(:,i)=sort(eig(K,M));       
    
    i=i+1;
    end
%--------------------------------------------------------------------------


% plot
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
figure(23)
hold on
grid on
for i=1:number_of_waves
% scatter(range_of_dimless_k,log10(real(frequencies(i,:))), "MarkerEdgeColor","#77AC30", Marker=".")
% plot(range_of_dimless_k,(real(frequencies(i,:))), "MarkerEdgeColor","#77AC30", Marker=".")
plot(range_of_dimless_k,(real(frequencies(i,:))), 'red', 'LineWidth', 2)
end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
% ylim([0 2000])
% hold off
%--------------------------------------------------------------------------

%% DISPERSION CURVE CLASSIC LONG WAVELENGTH + WILLIS STEADY TERMS + e-MODES


% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 1+n_modes; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.01;
range_of_dimless_k = -1:kstep:2;
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda      =   zeros(number_of_waves,length(range_of_dimless_k)); % 
Gamma_index =   round(length(range_of_dimless_k)/3)+1;
X_index     = 2*round(length(range_of_dimless_k)/3)+1;
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;% index of current iteration
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else 
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end
    
    % Density matrix (k)
    
    % second line of M
    Mpp=        heC+dot(k,hvD,k)      +1i*       (dot(-heD,k)+dot(k,hvC) );
    Mpe=        heL.'                         ;
    
    % third line of M
    Mep=      - hel                           ;
    Mee= I_Q_Q                                                            ;

    % Elasticity matrix (k)    
    % second line of K
    Kpp=        heG+dot(k,hvH,k)      +1i*       (dot(-heH,k)+dot(k,hvG)     );
    Kpe= O_Q_1.' ;
    
    % third line of K
    Kep= O_Q_1 ; 
    Kee= LAM_Q_Q ;                                                                    

    % density matrix 
    M=[Mpp Mpe;
       Mep Mee];
    % elasticity matrix 
    K=[Kpp Kpe;
       Kep Kee];
    
    % compute eigevalues
%     [V,D]=eig(K,M)
    lambda(:,i)=sort(eig(K,M));       

    i=i+1;
    end
%--------------------------------------------------------------------------


% plot
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
figure(23)
hold on
grid on
for i=1:number_of_waves
scatter(range_of_dimless_k,real(frequencies(i,:)),'m.')
end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 2000])
% hold off
%--------------------------------------------------------------------------
%% DISPERSION CURVE ALL  TERMS


% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 1+n_modes; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.01;
range_of_dimless_k = -1:kstep:2;
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda      =   zeros(number_of_waves,length(range_of_dimless_k)); % 
Gamma_index =   round(length(range_of_dimless_k)/3)+1;
X_index     = 2*round(length(range_of_dimless_k)/3)+1;
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;% index of current iteration
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else 
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end
    
    % Density matrix (k)
    
    % second line of M
    Mpp=        heC+dot(k,hvD,k)      +1i*       (dot(-heD,k)+dot(k,hvC) );
    Mpe=        heL.'                 +1i*        dot(k,hvL.')            ;
    
    % third line of M
    Mep=      - hel                   +1i*        dot( + hvl  ,k)         ;
    Mee= I_Q_Q                                                            ;

    % Elasticity matrix (k)    
    % second line of K
    Kpp=        heG+dot(k,hvH,k)      +1i*       (dot(-heH,k)+dot(k,hvG)     );
    Kpe= O_Q_1.' ;
    
    % third line of K
    Kep= O_Q_1 ; 
    Kee= LAM_Q_Q ;                                                                    

    % density matrix 
    M=[Mpp Mpe;
       Mep Mee];
    % elasticity matrix 
    K=[Kpp Kpe;
       Kep Kee];
    
    % compute eigevalues
%     [V,D]=eig(K,M)
    lambda(:,i)=sort(eig(K,M));       

    i=i+1;
    end
%--------------------------------------------------------------------------


% plot
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(23)
hold on
grid on
for i=1:number_of_waves
% scatter(range_of_dimless_k,real(frequencies(i,:)),'r')
plot(range_of_dimless_k,real(frequencies(i,:)),'red', 'LineWidth', 2)
end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 2000])
% hold off
%--------------------------------------------------------------------------

