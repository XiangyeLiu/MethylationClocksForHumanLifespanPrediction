import sys
import pandas as pd

args = sys.argv
filename = args[1]
mlfilename = args[2]
if args[3] == "1":
    septype = " "
elif args[3] == "2":
    septype = "\t"
elif args[3] == "3":
    septype = ","

if ".txt" in filename:
    scfilename = str.replace(filename, ".txt", "_screened.csv")
elif ".csv" in filename:
    scfilename = str.replace(filename, ".csv", "_screened.csv")

mlfile = pd.read_csv(mlfilename, header = None)
number = len(mlfile)
file = pd.read_csv(filename, header = 0, index_col = 0, sep = septype)

for i in range(number):
    site = mlfile.iat[i, 0]
    orimax = mlfile.iat[i, 1]
    orimin = mlfile.iat[i, 2]

    if site in file.index.values:
        sclist = file.loc[site]
        if max(sclist) < 0 or min(sclist) > 1:
            print("Site", i + 1, "Extremum Unaccepted")
        else:
            scmax = max([scvalue for scvalue in sclist if scvalue <= 1])
            scmin = min([scvalue for scvalue in sclist if scvalue >= 0])
            scline = pd.DataFrame(columns = ["site", "max", "min"])
            scline.loc[0] = [site, scmax, scmin]
            scline.to_csv(scfilename, header = False, index = False, mode = "a")
            print("Site", i + 1, "Screened")
            if scmax > orimax or scmin < orimin:
                if scmax > orimax:
                    mlfile.iat[i, 1] = scmax
                if scmin < orimin:
                    mlfile.iat[i, 2] = scmin
                mlfile.to_csv(mlfilename, header = False, index = False)
                print("Site", i + 1, "Extremum Changed")
            else:
                print("Site", i + 1, "Extremum Unchanged")
    else:
        print("Site", i + 1, "Unmatched")

print("All Finished")
