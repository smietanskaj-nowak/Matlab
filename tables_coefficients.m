function [pozycje, global_ADP, extinction, dw_phason]=tables_coefficients(coefficients, num_of_elem, N)
pozycje=zeros(N, num_of_elem);
for i=1:N
   for j=1:num_of_elem
       pozycje(i, j)=coefficients((i-1)*num_of_elem+j, 1);
   end
end

global_ADP=coefficients(end-2, 1);
extinction=coefficients(end-1, 1);
dw_phason=coefficients(end, 1);
end
