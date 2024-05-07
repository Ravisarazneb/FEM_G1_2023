clear all; close all;clear clc;

% Input modulus of elasticity and posson ratio
E0=22.36*10^9; %modulus of elasticity (Pa)
v=0.2;     %posson ratio 
g=-9.81;  %gravitational acceleration
SigmaMax= 20*10^6 %(Pa)
p = 2400; %kg/m3 density of concrete

% Input The cartesian coordinate for Tetrahedron with four nodes
%Uses pre-created mesh from gmsh
%Run mesh file name without .m:
tetra_size_one   % Name of the file without .m exported from Gmsh
coor = msh.POS;;
xp = coor(:,1)'; % [Node1 Node2 Node3 Node4 ]
yp = coor(:,2)'; % [Node1 Node2 Node3 Node4 ]
zp = coor(:,3)'; % [Node1 Node2 Node3 Node4 ]

ELEMCon=msh.TETS(:,1:4); %element connectivity
NE=size(ELEMCon,1);                %Number of elements
nNode=size(coor,1);                  %Number of nodes
nDof=nNode*3;                         %Total number degree of freedom

% Plot the tetrahedron using tetramesh
tetramesh(ELEMCon, coor);

% Customize the plot (optional)
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Tetrahedron');
%-------------------NODE container and ELEMENT container--------------------
for p=1:nNode
    NODE(p).X=xp(1,p);NODE(p).Y=yp(1,p);NODE(p).Z=zp(1,p);
end

for eNo=1:NE
    ELEMENT(eNo).con=ELEMCon(eNo,:);
end
%-------------------Input support matrix and Load matrix--------------------
%[ Node ux uy uz]          ux=1 Fixed ;ux=0 Free
support=[2 1 1 1; 
         4 1 1 1;
         6 1 1 1;
         8 1 1 1;
         10 1 1 1;
         12 1 1 1;
         13 1 1 1;
         16 1 1 1];
nSupport=size(support,1);

%IN put Load conditions [Node Fx Fy Fz]
load=[14   -1400000/2   0  0;
      15    -1400000/2  0  0]
nLoad=size(load,1);

%INput direction of force 2=X-direction ;3=Y-direction ;4=Z-direction 
F_Location=2

%INput Degree of freedome direction of interest to plot curve Force-Displacement
dof_direction=[3*load(1,1)-2;
               3*load(2,1)-2] 

%Shape Function
syms zeta xi eta
N1=zeta;
N2=eta;
N3=xi;
N4=1-xi-eta-zeta;
N=[N1 N2 N3 N4];
nN=length(N);

%To check the volume of strucutre
V_Ele=zeros(1,NE);

%Preparing to store each stiffness element in 3D matrix
kElem_temp=zeros(12,12,NE);
for eNo=1:NE
    ELEMENT(eNo).Ex=E0;
    ELEMENT(eNo).Ey=E0;
    ELEMENT(eNo).Ez=E0;
    ELEMENT(eNo).Exy=E0;
    ELEMENT(eNo).Eyz=E0;
    ELEMENT(eNo).Exz=E0;
end
uDisp=zeros(nDof,1);
NINC=30;                      %the number of calculation steps applied
uDel_inc=zeros(NINC,(nLoad));
uLoad_inc=zeros(NINC,nLoad);

for nin=1:NINC
for i = 1:NE
    %Mapping 
    x=(N1*coor(ELEMCon(i,1),1))+(N2*coor(ELEMCon(i,2),1))+(N3*coor(ELEMCon(i,3),1))+(N4*coor(ELEMCon(i,4),1));
    y=(N1*coor(ELEMCon(i,1),2))+(N2*coor(ELEMCon(i,2),2))+(N3*coor(ELEMCon(i,3),2))+(N4*coor(ELEMCon(i,4),2));
    z=(N1*coor(ELEMCon(i,1),3))+(N2*coor(ELEMCon(i,2),3))+(N3*coor(ELEMCon(i,3),3))+(N4*coor(ELEMCon(i,4),3));
    
    %Jacobian
    J=[diff(x,xi) diff(y,xi) diff(z,xi) ;
    diff(x,eta) diff(y,eta) diff(z,eta) ;
    diff(x,zeta) diff(y,zeta) diff(z,zeta)];
    % disp('The Jacobian is:');
    % J
    det_J=det(J);
    
    %Strain Displacement Matrix with dimension of 6x12
    %(Derivative of N with respect to x,y and z)
    dNs =zeros(3,4);
    for a =1:4
        dN =inv(J)*[diff(N(a),xi);diff(N(a),eta);diff(N(a),zeta)];
        dNs(:,a) = dN;
    end
    
    %Develop B matrix
    B=zeros(6,12);
    for b = 1:4
        B(:,3*b-2:3*b) = [dNs(1,b) 0 0;
                          0 dNs(2,b) 0;
                          0 0 dNs(3,b);
                          dNs(2,b) dNs(1,b) 0;
                          0 dNs(3,b) dNs(2,b);
                          dNs(3,b) 0 dNs(1,b)];
    end
    ELEMENT(i).B=B;
    % Elasticity matrix D
    d =(1+v)*(1-2*v);
    D = 1/d*[ELEMENT(i).Ex*(1-v) ELEMENT(i).Ey*v ELEMENT(i).Ez*v 0 0 0;
              ELEMENT(i).Ex*v ELEMENT(i).Ey*(1-v) ELEMENT(i).Ez*v 0 0 0;
              ELEMENT(i).Ex*v ELEMENT(i).Ey*v ELEMENT(i).Ez*(1-v) 0 0 0;
              0 0 0 ELEMENT(i).Exy*(1-2*v)/2 0 0;
              0 0 0 0 ELEMENT(i).Eyz*(1-2*v)/2 0;
              0 0 0 0 0 ELEMENT(i).Exz*(1-2*v)/2];
    
    % Volume of Tetrahedron
    V=det_J/6;
    V_Ele(1,i) = V;
    
    %local stifness matrix
    kElem=B'*D*B*V;
    % sprs_kElem=sparse(kElem)
    kElem_temp(:,:,i) = kElem;
    ELEMENT(i).stiffness=kElem;
end
kElem_temp;
%To counter check the negative stiffness effect by measuring the total
%volume of structure 
V_Ele;      %negative value indicates flipping of element, require reordering of nodes
Total_volume = sum(V_Ele);
%-------------------------------------------------------------------------%
%                  Global Stiffness Matrix Calculation                    %
%-------------------------------------------------------------------------%
% This part calculates the global stiffness matrix. Basically; for each
% element it takes 4x4 part from the element stiffness matrix and puts to
% the correct spot on the global stiffness matrix. This process loops until
% all elements all parts placed in to the global stiffness matrix.
KG = zeros(nDof,nDof);

for eNo=1:NE
    for j=1:nN
        for i=1:nN
            n = ELEMCon(eNo,i);
            m = ELEMCon(eNo,j);
            KG(3*n-2,3*m-2) = KG(3*n-2,3*m-2)+kElem_temp(3*i-2,3*j-2,eNo);              %K11
            KG(3*n-2,3*m-1) = KG(3*n-2,3*m-1)+kElem_temp(3*i-2,3*j-1,eNo);              %K12
            KG(3*n-2,3*m) = KG(3*n-2,3*m)+kElem_temp(3*i-2,3*j,eNo);                    %K13
            KG(3*n-1,3*m-2) = KG(3*n-1,3*m-2)+kElem_temp(3*i-1,3*j-2,eNo);              %K21
            KG(3*n-1,3*m-1) = KG(3*n-1,3*m-1)+kElem_temp(3*i-1,3*j-1,eNo);              %K22
            KG(3*n-1,3*m) = KG(3*n-1,3*m)+kElem_temp(3*i-1,3*j,eNo);                    %K23
            KG(3*n,3*m-2) = KG(3*n,3*m-2)+kElem_temp(3*i,3*j-2,eNo);                    %K13
            KG(3*n,3*m-1) = KG(3*n,3*m-1)+kElem_temp(3*i,3*j-1,eNo);                    %K23
            KG(3*n,3*m) = KG(3*n,3*m)+kElem_temp(3*i,3*j,eNo);                          %K33
        end
    end
end
KG=sparse(KG);
%% 
% 



%-------------------------------------------------------------------------%
%           Apply Support Conditions to the Global Stiffnes Matrix        %
%-------------------------------------------------------------------------%
% This part makes zeros all columns and rows where the supports are except
% the diagonal element of the matrix. Diagonal element set to 1. I choose
% this method because with this way sort of displacement evaulated are not
% changes. And later it is easier to use evaluated values. Only negative
% side of this approach is we have to be careful not to put force where the
% support is fixed.
for i=1:nSupport
    n = support(i,1);
    if (support(i,2) == 1) %%%ux =0
        KG(3*n-2,:) = 0;
        KG(:,3*n-2) = 0;
        KG(3*n-2,3*n-2) = 1;
    end
    if (support(i,3) == 1) %%%uy=0
        KG(3*n-1,:) = 0;
        KG(:,3*n-1) = 0;
        KG(3*n-1,3*n-1) = 1;
    end
    if   (support(i,4) == 1) %%%uz=0
        KG(3*n,:) = 0;
        KG(:,3*n) = 0;
        KG(3*n,3*n) = 1;
    end
end

KG;
%-------------------------------------------------------------------------%
%                       Load Vector Computation                           %
%-------------------------------------------------------------------------%
% In this part load traction vector created. If there is a load vector get the value
% from load matrix. If not not load value set to zero.
f_trac = zeros(nDof,1);
for i=1:nLoad
    n = load(i,1);
    f_trac(3*n-2) = load(i,2)/NINC;
    f_trac(3*n-1) = load(i,3)/NINC;
    f_trac(3*n) = load(i,4)/NINC;
end

%Body force applied at each node in z-axis
f_body = zeros(nDof,1);
for eNo=1:NE
        for i=1:nN
            n = ELEMCon(eNo,i);
            f_body(3*n,1) = (f_body(3*n,1)+g*p*V_Ele(1,eNo)/4)/NINC;
        end
end
fvec = f_trac;


% uDel_inc=KG\f_trac
uDisp= uDisp + KG\f_trac;

%uDel_inc =zeros(ninc,(nLoad))
    for q=1:(nLoad)          % Display all displacements on the free end
    %    fprintf('%10.7f \n',uDisp(load(i,1)*2));
        uDel_inc(nin,q)=uDisp(load(q,1)*3-2,1);
        if(nin==1)
            uLoad_inc(nin,q)=fvec(load(q,1)*3-2);
        else
            uLoad_inc(nin,q)=uLoad_inc(nin-1,q)+fvec(load(q,1)*3-2);
        end
    end
 %Calculation of stresses
    
    %Nodal displacements
    for i=1:nNode
        NODE(i).u=[uDisp(3*i-2); uDisp(3*i-1);uDisp(3*i)];
    end
    %Element displacements
    for i=1:NE
        temp=[];
        for j=1:nN
            temp=[temp;
                  NODE(ELEMCon(i,j)).u];
        end
                ELEMENT(i).u=temp;
    end
    %Element Stresses
    for i=1:NE
       uelem = ELEMENT(i).u;                                                % Retrieve the displacement vector for element 'i'
       Belem=ELEMENT(i).B;                                                  % Retrieve the B matrix for element 'i'
       sigma = D*Belem*uelem;                                               % Calculate stress using Hooke's law: stress = D * B * u            
       strain=Belem*uelem;                                                  % Calculate strain: strain = B * u
       for j=1:4
           ELEMENT(i).sigma(:,j) =sigma;                                    % Assigning 'sigma' to the j-th column of the 'sigma' field of ELEMENT(i)
           ELEMENT(i).strain(:,j) =strain;                                  % Assigning 'strain' to the j-th column of the 'strain' field of ELEMENT(i),This represents constant strain for the j-th column
       end
       Strain_x=mean(ELEMENT(i).strain(1,:));                               % Calculate the mean strain in the x-direction for the ith element
       Stress_x=mean(ELEMENT(i).sigma(1,:));                                % Calculate the mean stress in the x-direction for the ith element
       [Stress, Stiffness]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_x)); % Call the function GetHDStressStiffness to obtain the stress and stiffness for the material based on certain parameters (E0, SigmaMax, and the absolute value of Strain_x)
       ELEMENT(i).Actual_stress_x=Stress*sign(Strain_x);                    % Calculate the actual stress in the x-direction for the ith element,considering the sign of Strain_x
       ELEMENT(i).Ex=Stiffness;                                             % Assign the stiffness obtained from the function to the property Ex of the ith element

       Strain_y=mean(ELEMENT(i).strain(2,:));                               % Calculate the mean strain in the y-direction for the ith element
       Stress_y=mean(ELEMENT(i).sigma(2,:));                                % Calculate the mean stress in the y-direction for the ith element
       [Stress, Stiffness]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_y)); % Call the function GetHDStressStiffness to obtain the stress and stiffness for the material based on certain parameters (E0, SigmaMax, and the absolute value of Strain_y)
       ELEMENT(i).Actual_stress_y=Stress*sign(Strain_y);                    % Calculate the actual stress in the y-direction for the ith element,considering the sign of Strain_y
       ELEMENT(i).Ey=Stiffness;                                             % Assign the stiffness obtained from the function to the property Ey of the ith element
       
       Strain_z=mean(ELEMENT(i).strain(3,:));                               % Calculate the mean strain in the z-direction for the ith element
       Stress_z=mean(ELEMENT(i).sigma(3,:));                                % Calculate the mean stress in the z-direction for the ith element
       [Stress, Stiffness]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_z)); % Call the function GetHDStressStiffness to obtain the stress and stiffness for the material based on certain parameters (E0, SigmaMax, and the absolute value of Strain_z)
       ELEMENT(i).Actual_stress_z=Stress*sign(Strain_z);                    % Calculate the actual stress in the z-direction for the ith element,considering the sign of Strain_z
       ELEMENT(i).Ez=Stiffness;                                             % Assign the stiffness obtained from the function to the property Ez of the ith element

       Strain_xy=mean(ELEMENT(i).strain(4,:));                              % Calculate the mean strain in the xy-plane for the ith element
       Stress_xy=mean(ELEMENT(i).sigma(4,:));                               % Calculate the mean stress in the xy-plane for the ith element
       [Stress, Stiffness]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_xy));% Call the function GetHDStressStiffness to obtain the stress and stiffness for the material based on certain parameters (E0, SigmaMax, and the absolute value of Strain_xy)
       ELEMENT(i).Actual_stress_xy=Stress*sign(Strain_xy);                  % Calculate the actual stress in the xy-plane for the ith element,considering the sign of Strain_xy
       ELEMENT(i).Exy=Stiffness;                                            % Assign the stiffness obtained from the function to the property Ez of the ith element
       
       Strain_yz=mean(ELEMENT(i).strain(5,:));                              % Calculate the mean strain in the yz-plane for the ith element
       Stress_yz=mean(ELEMENT(i).sigma(5,:));                               % Calculate the mean stress in the yz-plane for the ith element
       [Stress, Stiffness]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_yz));% Call the function GetHDStressStiffness to obtain the stress and stiffness for the material based on certain parameters (E0, SigmaMax, and the absolute value of Strain_yz)
       ELEMENT(i).Actual_stress_yz=Stress*sign(Strain_yz);                  % Calculate the actual stress in the yz-plane for the ith element,considering the sign of Strain_yz
       ELEMENT(i).Eyz=Stiffness;                                            % Assign the stiffness obtained from the function to the property Eyz of the ith element
       
       Strain_xz=mean(ELEMENT(i).strain(6,:));                              % Calculate the mean strain in the xz-plane for the ith element
       Stress_xz=mean(ELEMENT(i).sigma(6,:));                               % Calculate the mean stress in the xz-plane for the ith element
       [Stress, Stiffness]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_xz));% Call the function GetHDStressStiffness to obtain the stress and stiffness for the material based on certain parameters (E0, SigmaMax, and the absolute value of Strain_xz)
       ELEMENT(i).Actual_stress_xz=Stress*sign(Strain_xz);                  % Calculate the actual stress in the xz-plane for the ith element,considering the sign of Strain_xz
       ELEMENT(i).Exz=Stiffness;                                            % Assign the stiffness obtained from the function to the property Exz of the ith element
       
    end
    if(nin==26) %If the condition 'nin==26' is met the requirement,Copy the actual values of stress and displacement for each element into their respective required fields.
       uglobal = uDisp
        for i = 1:NE
            ELEMENT(i).u_req = ELEMENT(i).u;                               % Displacement in x, y, and z directions
            ELEMENT(i).stress_x_req = ELEMENT(i).Actual_stress_x;          % Stress in the x-direction
            ELEMENT(i).stress_y_req = ELEMENT(i).Actual_stress_y;          % Stress in the y-direction
            ELEMENT(i).stress_z_req = ELEMENT(i).Actual_stress_z;          % Stress in the z-direction
            ELEMENT(i).stress_xy_req = ELEMENT(i).Actual_stress_xy;        % Shear stress in the xy-plane
            ELEMENT(i).stress_yz_req = ELEMENT(i).Actual_stress_yz;        % Shear stress in the yz-plane
            ELEMENT(i).stress_xz_req = ELEMENT(i).Actual_stress_xz;        % Shear stress in the xz-plane
        end
    end
    LoadStep=nin
end
%create incremental Load vector
%F=(abs(load(1,F_Location)/NINC)):(abs(load(1,F_Location)/NINC)):abs((load(1,F_Location)))
%Force Displacement Graph
%Force Displacement Graph
figure(2)
plot(abs(uDel_inc),abs(uLoad_inc));
xlabel('Displacement(m.)')
ylabel('Force(N)')
%%%%%%%%%%%%%%%%%%% POST PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to plot deform shape of struture
for i = 1:nNode
    udisp_coor(i,:) = [uglobal(3*i-2,1) uglobal(3*i-1,1) uglobal(3*i,1)];
end
deform_coor=coor+udisp_coor;
% Plot the tetrahedron using tetramesh
figure(3)
subplot(211)
tetramesh(ELEMCon, coor);
% Customize the plot (optional)
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Tetrahedron_before');
subplot(212)
tetramesh(ELEMCon, deform_coor);
% Customize the plot (optional)
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Tetrahedron_after');