clear
clc
N=5647;
path='Structural_data.txt';
[entry_id, cell, symmetry_group, Components, Atoms]=load_structure_data(path);
[av, bv, cv]=calc_basevectors(cell(1), cell(2), cell(3)/7, cell(4), cell(5), cell(6));
bonds=load('Bonds.txt');
hydrogen=load('Wodory_Hyp_7.txt');
mol_list=load('Lista_molekul_atom_index.txt');
coefficients=load('coefficients_Hyp_7.txt');
[pozycje, global_ADP, extinction, dw_phason]=tables_coefficients(coefficients(:, 1), 14, N);
pozycje(:, 3:5)=pozycje(:, 3:5)*[av; bv; cv];

fp=fopen('Lista_atomow_id.txt');
licz_atom=0;
while ~feof(fp)
    licz_atom=licz_atom+1;
   line=fgetl(fp);
   words=strsplit(line);
   atomy_id(licz_atom, :)={str2num(words{1}), words{2}, words{3}, words{4}, str2num(words{5})};
end
fclose(fp);

choice=1:600;

[selected, bond, inter_mol]=select_mol(mol_list, atomy_id, bonds, choice);
x=pozycje(selected,3);
y=pozycje(selected,4);
z=pozycje(selected,5);

wodory_selected=select_hydrogen(hydrogen, choice);

figure
hold on

%%%%%%%%%%%Spheres
color_N='b';
color_C='k';
color_O='r';
color_S='y';
SPHERE_RES = 8; % resolution for your spheres
CYLINDER_RES=8;
SPHERE_RAD = 0.5; % radius of spheres

X_sphere=zeros(length(selected)*(SPHERE_RES+1)+length(selected), SPHERE_RES+1);
Y_sphere=zeros(length(selected)*(SPHERE_RES+1)+length(selected), SPHERE_RES+1);
Z_sphere=zeros(length(selected)*(SPHERE_RES+1)+length(selected), SPHERE_RES+1);
C_sphere=zeros(length(selected)*(SPHERE_RES+1)+length(selected), SPHERE_RES+1, 3);

X_cyl=zeros((length(bond)+length(inter_mol))*(4)+2*(length(bond)+length(inter_mol)), CYLINDER_RES);
Y_cyl=zeros((length(bond)+length(inter_mol))*(4)+2*(length(bond)+length(inter_mol)), CYLINDER_RES);
Z_cyl=zeros((length(bond)+length(inter_mol))*(4)+2*(length(bond)+length(inter_mol)), CYLINDER_RES);
C_cyl=zeros((length(bond)+length(inter_mol))*(4)+2*(length(bond)+length(inter_mol)), CYLINDER_RES, 3);

X_sphere_h=zeros(length(wodory_selected(:, 1))*(SPHERE_RES+1)+length(wodory_selected(:, 1)), SPHERE_RES+1);
Y_sphere_h=zeros(length(wodory_selected(:, 1))*(SPHERE_RES+1)+length(wodory_selected(:, 1)), SPHERE_RES+1);
Z_sphere_h=zeros(length(wodory_selected(:, 1))*(SPHERE_RES+1)+length(wodory_selected(:, 1)), SPHERE_RES+1);
C_sphere_h=zeros(length(wodory_selected(:, 1))*(SPHERE_RES+1)+length(wodory_selected(:, 1)), SPHERE_RES+1);

X_cyl_h=zeros((length(wodory_selected(:, 1)))*(4)+2*(length(wodory_selected(:, 1))), CYLINDER_RES);
Y_cyl_h=zeros((length(wodory_selected(:, 1)))*(4)+2*(length(wodory_selected(:, 1))), CYLINDER_RES);
Z_cyl_h=zeros((length(wodory_selected(:, 1)))*(4)+2*(length(wodory_selected(:, 1))), CYLINDER_RES);
C_cyl_h=zeros((length(wodory_selected(:, 1)))*(4)+2*(length(wodory_selected(:, 1))), CYLINDER_RES, 3);

[xb, yb, zb] = sphere(SPHERE_RES);
xb=single(xb); yb=single(yb); zb=single(zb);

list_index=1;
for i=1:length(x)
        X_sphere(list_index:(SPHERE_RES+list_index), :)=SPHERE_RAD*xb+x(i, 1);
        X_sphere(list_index+SPHERE_RES+1, :)=NaN;
        Y_sphere(list_index:(SPHERE_RES+list_index), :)=SPHERE_RAD*yb+y(i, 1);
        Y_sphere(list_index+SPHERE_RES+1, :)=NaN;
        Z_sphere(list_index:(SPHERE_RES+list_index), :)=SPHERE_RAD*zb+z(i, 1);
        Z_sphere(list_index+SPHERE_RES+1, :)=NaN;
    if pozycje(selected(i), 6)==1
        C_sphere(list_index:(SPHERE_RES+list_index), :, 1)=0;
        C_sphere(list_index:(SPHERE_RES+list_index), :, 2)=0;
        C_sphere(list_index:(SPHERE_RES+list_index), :, 3)=1;
        C_sphere(list_index+SPHERE_RES+1, :, 1:3)=NaN;
    end
    if pozycje(selected(i), 7)==1
        C_sphere(list_index:(SPHERE_RES+list_index), :, 1)=0;
        C_sphere(list_index:(SPHERE_RES+list_index), :, 2)=0;
        C_sphere(list_index:(SPHERE_RES+list_index), :, 3)=0;
        C_sphere(list_index+SPHERE_RES+1, :, 1:3)=NaN;
    end
    if pozycje(selected(i), 8)==1
        C_sphere(list_index:(SPHERE_RES+list_index), :, 1)=1;
        C_sphere(list_index:(SPHERE_RES+list_index), :, 2)=0;
        C_sphere(list_index:(SPHERE_RES+list_index), :, 3)=0;
        C_sphere(list_index+SPHERE_RES+1, :, 1:3)=NaN;
    end
    if pozycje(selected(i), 9)==1
        C_sphere(list_index:(SPHERE_RES+list_index), :, 1)=1;
        C_sphere(list_index:(SPHERE_RES+list_index), :, 2)=1;
        C_sphere(list_index:(SPHERE_RES+list_index), :, 3)=0;
        C_sphere(list_index+SPHERE_RES+1, :, 1:3)=NaN;
    end
    list_index=list_index+SPHERE_RES+1+1;
end
surf(X_sphere, Y_sphere, Z_sphere, C_sphere, 'edgealpha', 0);
list_index_cyl=1;
%%%%%%%%Cylinders as Bonds
for i=1:length(bond)
    [xc,yc,zc]=cylinder2P([SPHERE_RAD/2,SPHERE_RAD/2],CYLINDER_RES,pozycje(bonds(bond(i), 2), 3:5), 0.5*pozycje(bonds(bond(i), 3), 3:5)+0.5*pozycje(bonds(bond(i), 2), 3:5));
    X_cyl(list_index_cyl:list_index_cyl+1, :)=xc;
    Y_cyl(list_index_cyl:list_index_cyl+1, :)=yc;
    Z_cyl(list_index_cyl:list_index_cyl+1, :)=zc;
    X_cyl(list_index_cyl+2, :)=NaN;
    Y_cyl(list_index_cyl+2, :)=NaN;
    Z_cyl(list_index_cyl+2, :)=NaN;
    if pozycje(bonds(bond(i), 2), 6)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=1;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(bonds(bond(i), 2), 7)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(bonds(bond(i), 2), 8)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(bonds(bond(i), 2), 9)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    list_index_cyl=list_index_cyl+3;
    
    [xc,yc,zc]=cylinder2P([SPHERE_RAD/2,SPHERE_RAD/2],CYLINDER_RES,pozycje(bonds(bond(i), 3), 3:5), 0.5*pozycje(bonds(bond(i), 2), 3:5)+0.5*pozycje(bonds(bond(i), 3), 3:5));
    X_cyl(list_index_cyl:list_index_cyl+1, :)=xc;
    Y_cyl(list_index_cyl:list_index_cyl+1, :)=yc;
    Z_cyl(list_index_cyl:list_index_cyl+1, :)=zc;
    X_cyl(list_index_cyl+2, :)=NaN;
    Y_cyl(list_index_cyl+2, :)=NaN;
    Z_cyl(list_index_cyl+2, :)=NaN;
    if pozycje(bonds(bond(i), 3), 6)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=1;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(bonds(bond(i), 3), 7)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(bonds(bond(i), 3), 8)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(bonds(bond(i), 3), 9)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    list_index_cyl=list_index_cyl+3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if ~isempty(inter_mol)
for i=1:length(inter_mol(:, 1))
    [xc,yc,zc]=cylinder2P([SPHERE_RAD/2,SPHERE_RAD/2],CYLINDER_RES,pozycje(inter_mol(i, 3), 3:5), 0.5*pozycje(inter_mol(i, 3), 3:5)+0.5*pozycje(inter_mol(i, 4), 3:5));
    X_cyl(list_index_cyl:list_index_cyl+1, :)=xc;
    Y_cyl(list_index_cyl:list_index_cyl+1, :)=yc;
    Z_cyl(list_index_cyl:list_index_cyl+1, :)=zc;
    X_cyl(list_index_cyl+2, :)=NaN;
    Y_cyl(list_index_cyl+2, :)=NaN;
    Z_cyl(list_index_cyl+2, :)=NaN;
    if pozycje(inter_mol(i, 3), 6)==1
        color1=color_N;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=1;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(inter_mol(i, 3), 7)==1
        color1=color_C;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(inter_mol(i, 3), 8)==1
        color1=color_O;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(inter_mol(i, 3), 9)==1
        color1=color_S;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    list_index_cyl=list_index_cyl+3;
    
    [xc,yc,zc]=cylinder2P([SPHERE_RAD/2,SPHERE_RAD/2],CYLINDER_RES,pozycje(inter_mol(i, 4), 3:5), 0.5*pozycje(inter_mol(i, 4), 3:5)+0.5*pozycje(inter_mol(i, 3), 3:5));
    X_cyl(list_index_cyl:list_index_cyl+1, :)=xc;
    Y_cyl(list_index_cyl:list_index_cyl+1, :)=yc;
    Z_cyl(list_index_cyl:list_index_cyl+1, :)=zc;
    X_cyl(list_index_cyl+2, :)=NaN;
    Y_cyl(list_index_cyl+2, :)=NaN;
    Z_cyl(list_index_cyl+2, :)=NaN;
    if pozycje(inter_mol(i, 4), 6)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=1;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(inter_mol(i, 4), 7)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(inter_mol(i, 4), 8)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=0;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    if pozycje(inter_mol(i, 4), 9)==1
        C_cyl(list_index_cyl:list_index_cyl+1, :, 1)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 2)=1;
        C_cyl(list_index_cyl:list_index_cyl+1, :, 3)=0;
        C_cyl(list_index_cyl+2, :, 3)=NaN;
    end
    list_index_cyl=list_index_cyl+3;
end
end
surf(X_cyl, Y_cyl, Z_cyl, C_cyl, 'edgealpha', 0);

list_index_h=1;
list_index_cyl_h=1;
for i=1:length(wodory_selected(:, 1))
    [xc,yc,zc]=cylinder2P([SPHERE_RAD/2,SPHERE_RAD/2],CYLINDER_RES,pozycje(wodory_selected(i, 1), 3:5), 0.5*pozycje(wodory_selected(i, 1), 3:5)+0.5*wodory_selected(i, 2:4));
    X_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :)=xc;
    Y_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :)=yc;
    Z_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :)=zc;
    X_cyl_h(list_index_cyl_h+2, :)=NaN;
    Y_cyl_h(list_index_cyl_h+2, :)=NaN;
    Z_cyl_h(list_index_cyl_h+2, :)=NaN;
    if pozycje(wodory_selected(i, 1), 6)==1
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 1)=0;
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 2)=0;
        C_cyl_H(list_index_cyl_h:list_index_cyl_h+1, :, 3)=1;
        C_cyl_h(list_index_cyl_h+2, :, 3)=NaN;
    end
    if pozycje(wodory_selected(i, 1), 7)==1
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 1)=0;
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 2)=0;
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 3)=0;
        C_cyl_h(list_index_cyl_h+2, :, 3)=NaN;
    end
    if pozycje(wodory_selected(i, 1), 8)==1
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 1)=1;
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 2)=0;
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 3)=0;
        C_cyl_h(list_index_cyl_h+2, :, 3)=NaN;
    end
    if pozycje(wodory_selected(i, 1), 9)==1
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 1)=1;
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 2)=1;
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 3)=0;
        C_cyl_h(list_index_cyl_h+2, :, 3)=NaN;
    end
    list_index_cyl_h=list_index_cyl_h+3;

    [xc,yc,zc]=cylinder2P([SPHERE_RAD/2,SPHERE_RAD/2],CYLINDER_RES,wodory_selected(i, 2:4), 0.5*wodory_selected(i, 2:4)+0.5*pozycje(wodory_selected(i, 1), 3:5));
    X_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :)=xc;
    Y_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :)=yc;
    Z_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :)=zc;
    X_cyl_h(list_index_cyl_h+2, :)=NaN;
    Y_cyl_h(list_index_cyl_h+2, :)=NaN;
    Z_cyl_h(list_index_cyl_h+2, :)=NaN;
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 1)=0.5;
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 2)=0.5;
        C_cyl_h(list_index_cyl_h:list_index_cyl_h+1, :, 3)=0.5;
        C_cyl_h(list_index_cyl_h+2, :, 3)=NaN;
        list_index_cyl_h=list_index_cyl_h+3;
    
     X_sphere_h(list_index_h:(SPHERE_RES+list_index_h), :)=SPHERE_RAD*xb+wodory_selected(i, 2);
     X_sphere_h(list_index_h+SPHERE_RES+1, :)=NaN;
     Y_sphere_h(list_index_h:(SPHERE_RES+list_index_h), :)=SPHERE_RAD*yb+wodory_selected(i, 3);
     Y_sphere_h(list_index_h+SPHERE_RES+1, :)=NaN;
     Z_sphere_h(list_index_h:(SPHERE_RES+list_index_h), :)=SPHERE_RAD*zb+wodory_selected(i, 4);
     Z_sphere_h(list_index_h+SPHERE_RES+1, :)=NaN;
        C_sphere_h(list_index_h:(SPHERE_RES+list_index_h), :, 1)=0.5;
        C_sphere_h(list_index_h:(SPHERE_RES+list_index_h), :, 2)=0.5;
        C_sphere_h(list_index_h:(SPHERE_RES+list_index_h), :, 3)=0.5;
        C_sphere_h(list_index_h+SPHERE_RES+1, :, 1:3)=NaN;
        list_index_h=list_index_h+SPHERE_RES+1+1;

end
surf(X_cyl_h, Y_cyl_h, Z_cyl_h, C_cyl_h, 'edgealpha', 0);
surf(X_sphere_h, Y_sphere_h, Z_sphere_h, C_sphere_h, 'edgealpha', 0);



camlight; lighting gouraud;
axis tight
daspect([1 1 1]);
xlabel('$X [\AA]$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$Y [\AA]$', 'Interpreter', 'Latex', 'FontSize', 15);
zlabel('$Z [\AA]$', 'Interpreter', 'Latex', 'FontSize', 15);
set(gca, 'FontSize', 15);
