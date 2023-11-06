import os

def pvtu_files(dir0):
    fname = os.listdir(dir0)
    filenames,foldernames = [],[]

    for name in fname:
        if os.path.isfile(dir0+name):
            filenames.append(name)
        if os.path.isdir(dir0+name):
            foldernames.append(name)

    pvtus = []

    for name in filenames:
        if ".pvtu"in name:
            pvtus.append(name)

    pvtus = sorted(pvtus, key=lambda x: int(x.split('.')[1]))
    return pvtus
