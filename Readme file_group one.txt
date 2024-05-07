# Project: Tetrahedron Finite Element Analysis (FEA) Solver

## Description
This MATLAB script performs a Finite Element Analysis (FEA) on a tetrahedral structure subjected to external static loads and support conditions. It calculates nodal displacements to get the displacement that is equal to 4% of height of the structure. Besides, this file is customized for simple and regular 3D structure with isotropic material and also capable to capture non-linearity elasticity behavior. External meshing software is needed to perform structure meshing.

## Usage
1. Open MATLAB.
2. use the script into a MATLAB editor.
3. Modify the input parameters according to your specific problem.
4. Run the script.

## Input Parameters
- Modulus of Elasticity (E0): The material's modulus of elasticity in Pascal (Pa).
- Poisson's Ratio (v): The material's Poisson's ratio.
- Maximum Allowable Stress (SigmaMax): The maximum allowable stress in Pascal (Pa).
- Density of Concrete (p): The density of the concrete material in kilograms per cubic meter (kg/m^3).
- Cartesian Coordinates of Tetrahedron Nodes (coor): The coordinates of the four nodes defining each tetrahedral element.
- Element Connectivity (ELEMCon): The connectivity matrix defining the nodes of each tetrahedral element.
- Support Conditions (support): Defines which degrees of freedom are fixed for each node.
- External Loads (static load): Specifies the applied loads at certain nodes.

#Process
-step 1 Use external meshing software to generate a mesh of muto frame using tetrahedral elements base on Delaunay 3D algorithm
-step 2 Write MATLAB function to read GMSH file format and convert mesh data [Nodal coordinates and element connectivity matrix] and make the node container and element container
-step 3 input E0, v, sigma max and apply displacement boundary condition and force boundary condition
-step 4 Derive stiffness matrix and force vector for each of element incorporating material and geometric properties using iso-parametric formulation.
-step 5 Solve the equilibrium equations to obtain nodal displacements
-step 6 Solve strain for each of element and use HD-model for solve stress and new E for each direction for next step because the strain from each direction is not same.
-step 7 Modify D matrix from new E
-step 8 Loop from step 5 that calculate with new D matrix.
-step 9 After finish all loop, look at uDel_inc and uLoad_inc in workspace to find the displacement,load that we can get the displacement and which step that we will get that value. After 	that call all the desired displacement and stress in that step to see the result.

Example : in matrix uDel_inc =0.058 at step 26 is value that we want so use code to call the data 

if(nin==26)
        for i = 1:NE
            ELEMENT(i).u_req = ELEMENT(i).u;
            ELEMENT(i).stress_x_req = ELEMENT(i).Actual_stress_x;
            ELEMENT(i).stress_y_req = ELEMENT(i).Actual_stress_y;
            ELEMENT(i).stress_z_req = ELEMENT(i).Actual_stress_z;
            ELEMENT(i).stress_xy_req = ELEMENT(i).Actual_stress_xy;
            ELEMENT(i).stress_yz_req = ELEMENT(i).Actual_stress_yz;
            ELEMENT(i).stress_xz_req = ELEMENT(i).Actual_stress_xz;
        end
    end

From that particular desired loop step, the global displacement is stored and then applied on original nodes to obtain deformed nodes location.

## Output
The script provides the following outputs:
- Nodal displacements.
- Elemental stresses and strains.
- Plot of the tetrahedral structure before and after deforming shapes.

## Notes
- Ensure correct input of node coordinates and element connectivity for accurate results.
- Modify the input parameters and boundary conditions according to your specific problem.

## ASIAN INSTITUTE OF TECHNOLOGY
## FIRST YEAR MASTER STUDENT, STRUCTURAL DEPARTMENT
## THIS FILE IS SOLELY CONTRIBUTED BY GROUP 1 BATCH 2023 AUG INTAKE
