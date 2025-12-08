import os
from pathlib import Path


def getFileNames(extension=None, path=None, include_hidden=False):
    """
    Returns a list containing the filenames in the working directory or somewhere else
    """
    if ( extension is not None or extension != '' ) and not extension.startswith('.'):
        extension = '.' + extension
    elif extension is None:
        extension = ''
    if path is None:  # Use current dir
        path = os.getcwd()
    path_obj = Path(path)
    assert path_obj.is_dir(), "The path argument is invalid."
    # Filter out directories from the list
    return sorted(
        file.name
        for file in path_obj.iterdir()
        if file.is_file()
        and (extension is None or file.suffix == extension)
        and (include_hidden or not file.name.startswith('.'))
    )


#def getFileNames(extension=None, path=None, include_hidden=False):
#    """
#    Returns a list containing the filenames in the working directory or somewhere else
#    """
#    if path is None:  # Use current dir
#        path = os.getcwd()
#    else:
#        assert os.path.exists(path)
#    files_in_directory = os.listdir(path)
#    # Filter out directories from the list
#    if extension is None:
#        filenames = [file for file in files_in_directory
#                     if os.path.isfile(os.path.join(path, file))]
#    else:
#        filenames = [file for file in files_in_directory
#                     if os.path.isfile(os.path.join(path, file))
#                     and file.endswith(extension)]
#    if not include_hidden:
#        for index in range(len(filenames))[::-1]:
#            if filenames[index].startswith('.'):
#                del filenames[index]
#    filenames.sort()
#    return filenames
#def getFileNames(extension):
#    # List all files in the current working directory
#    current_directory = os.getcwd()
#    files_in_directory = os.listdir(current_directory)
#    # Filter out directories from the list
#    filenames = [file for file in files_in_directory
#                 if os.path.isfile(os.path.join(current_directory, file))
#                 and file.endswith(extension)]
#    filenames.sort()
#    return filenames
