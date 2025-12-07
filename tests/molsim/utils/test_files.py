from molsim.utils.files import getFileNames
import tempfile
from unittest.mock import patch
from pathlib import Path

def test_getFileNames():
    # Create some test files
    txt_files = ["file1.txt", "file2.txt", "file4.txt"]
    log_files = ["file3.log"]
    filenames = txt_files + log_files
    with tempfile.TemporaryDirectory() as tempdir:
        path = Path(tempdir)
        for f in filenames:
            (path / f).touch()

        with patch("os.getcwd", return_value=tempdir):
            # Test getting all .txt files
            returned_txt_files = getFileNames(".txt")
            assert set(txt_files) == set(returned_txt_files)

            # Test getting all .log files
            returned_log_files = getFileNames(extension="log")
            assert isinstance(returned_log_files, list)
            assert log_files == returned_log_files

        # Test getting all files without filtering by extension
        #all_files = getFileNames(tmp_path)
        #assert set(all_files) == {str(tmp_path / f) for f in filenames}A

