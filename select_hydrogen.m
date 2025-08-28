function wodory_selected=select_hydrogen(wodory, choice)
licz_wodory=0;
for i=1:length(wodory(:, 1))
    for j=1:length(choice)
        if wodory(i, 2)==choice(j)
            licz_wodory=licz_wodory+1;
            wodory_selected(licz_wodory, :)=[wodory(i, 7), wodory(i, 4:6)];
        end
    end
end

end