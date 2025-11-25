def getFileNames(extension):
    # List all files in the current working directory
    current_directory = os.getcwd()
    files_in_directory = os.listdir(current_directory)
    # Filter out directories from the list
    filenames = [file for file in files_in_directory
                 if os.path.isfile(os.path.join(current_directory, file))
                 and file.endswith(extension)]
    filenames.sort()
    return filenames
