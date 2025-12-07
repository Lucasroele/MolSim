from molsim.parsers.xvgParser import XvgParser
from molsim.utils.files import getFileNames


def test_XvgParser(data_dir):
    data_1 = XvgParser("xvgs/dist_1.xvg", path=data_dir)
    filepaths = []
    for f in getFileNames(".xvg", path=data_dir / "xvgs"):
        filepaths.append(data_dir / "xvgs" / f)
    assert data_1.filenames == ["xvgs/dist_1.xvg"]

    data_2 = XvgParser(filepaths)
    assert data_2.filenames[0] == "dist_1.xvg"
    assert len(data_2) == len(filepaths)
    
    first_times = []
    for dataset in data_2:
        first_times.append(dataset[0][0])
    assert first_times == [0.001, 0.002, 0.003, 0.004, 0.005, 0.006]
    


    
