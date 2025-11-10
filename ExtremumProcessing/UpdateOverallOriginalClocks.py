import sys
import pandas as pd

args = sys.argv
newfilename = args[1]
orifilename = args[2]

newfile = pd.read_csv(newfilename, header = None)
orifile = pd.read_csv(orifilename, header = None)
newfile = newfile.set_index(newfile.columns[0], drop = False)
number = len(orifile)

for i in range(number):
    site = orifile.iat[i, 0]
    orimax = orifile.iat[i, 1]
    orimin = orifile.iat[i, 2]

    if site in newfile.index.values:
        newmax = newfile.loc[site].values[1]
        newmin = newfile.loc[site].values[2]
        if newmax > orimax or newmin < orimin:
            if newmax > orimax:
                orifile.iat[i, 1] = newmax
            if newmin < orimin:
                orifile.iat[i, 2] = newmin
            orifile.to_csv(orifilename, header = False, index = False)
            print("Site", i + 1, "Updated")
        else:
            print("Site", i + 1, "Maintained")
    else:
        print("Site", i + 1, "Maintained")

print("All Finished")
