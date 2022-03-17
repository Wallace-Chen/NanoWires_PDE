% a script code to save all current variable to the file.

name_of_save_file = sprintf('%d_charge_neutrality', degree_of_polygon);
for i = 1 : number_of_interfaces
    name_of_save_file = [name_of_save_file, '_', num2str(vector_of_side_lengths(i))];
end
save([name_of_save_file, '.mat']);
fprintf("Variables saved to %s.\n", [name_of_save_file, '.mat']);