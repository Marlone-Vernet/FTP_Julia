""" CHECK avancement de la ftp - analyse """


function count_elements_in_folder(folder_path::String)::Int
    # Check if the provided path is a directory
    if !isdir(folder_path)
        throw(ArgumentError("The provided path is not a directory: $folder_path"))
    end
    
    # Get the list of elements in the directory
    elements = readdir(folder_path)
    
    # Return the number of elements
    return length(elements)
end

path = "E:/DATA_FTP/150425/h_map_20Hz"
n = 4000

for i=1:n

    m = count_elements_in_folder(path)
    percentage_m = round(m*100/n)
    print("\rAvancement : $(percentage_m) %")

    if m==n
        break
    end

    sleep(10)

end
