function [entry_id, cell, symmetry_group, Components, Atoms]=load_structure_data(path)

fp=fopen(path, 'r');
while ~feof(fp)
newline=fgetl(fp);
    if strcmp(newline, 'entry_id')==1
        entry_id=fgetl(fp);
    end
    if strcmp(newline, 'cell_data')==1
        for i=1:6
            words=strsplit(fgetl(fp));
            cell(i)=str2num(words{2});
        end
        words=strsplit(fgetl(fp));
        symmetry_group=str2num(words{2});
    end
    if strcmp(newline, 'Components')==1
        Components=strsplit(fgetl(fp));
    end
    if strcmp(newline, 'Atoms')
        Atoms=strsplit(fgetl(fp));
    end
end

end
        
        