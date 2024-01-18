from main import main
from matplotlib import pyplot as plt

DirList=["Measurement_Data/Glock/Data/"]

MaxOmega=1e6
for directory in DirList:
    print(directory)
    main(directory,MaxOmega)
    plt.show()