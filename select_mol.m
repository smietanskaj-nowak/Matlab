function [list, bond, inter_mol]=select_mol(list_mol, atomy_id, bonds, choice)
rozmiar=0;
info=fopen('Konwersja_data.txt', 'r');
while ~feof(info)
   line=fgetl(info);
   words=strsplit(line);
   if strcmp(words{1}, 'Liczba_molekul_aminokwasow')
       amin_lim=str2num(words{5});
   end
end
for i=1:length(choice)
    atomy=list_mol(choice(i), 2):list_mol(choice(i), 3);
    rozmiar_p=rozmiar;
    rozmiar=rozmiar+length(atomy);
list(rozmiar_p+1:rozmiar)=atomy;
end

licz_bond=0;
for i=1:length(bonds)
    for j=1:length(choice)
        if bonds(i, 1)==choice(j)
            licz_bond=licz_bond+1;
          bond(licz_bond)=i;
        end
    end
end
licz_inter=0;
inter_mol=[];
for i=1:length(choice)-1
    if ((choice(1, i+1)-1)==choice(1, i) && choice(1, i)<=amin_lim)
        flag_p=0;
        flag_k=0;
        for j=list_mol(choice(i), 2):list_mol(choice(i), 3)
            if strcmp(atomy_id{j, 3}, 'C')==1
                pocz=j;
                flag_p=1;
            end
        end
        for j=list_mol(choice(i+1), 2):list_mol(choice(i+1), 3)
            if strcmp(atomy_id{j, 3}, 'N')==1
                kon=j;
                flag_k=1;
            end
        end
        if flag_p==1 && flag_k==1
            if atomy_id{pocz, 5}==atomy_id{kon, 5}
        licz_inter=licz_inter+1;
        inter_mol(licz_inter, :)=[list_mol(choice(i), 1) list_mol(choice(i+1), 1) pocz kon];
            end
        end
    end
end
fclose(info);
end