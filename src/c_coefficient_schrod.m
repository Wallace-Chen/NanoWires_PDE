function c_coeff_schrod = c_coefficient_schrod(sd)

    global vector_of_masses;
    %%
    
    c_coeff_schrod = 1 ./ vector_of_masses(sd);
    
    %% debug purpose
    global save_file;
    save(save_file, 'sd', 'c_coeff_schrod', '-append');
    %%
    
end