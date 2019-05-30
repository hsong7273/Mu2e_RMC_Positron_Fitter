import csv #For Converting csv file to space-separated-value file
 # RooDataSet Read takes uses space as the delimiter


filename = "fits_kmax_90.9.csv"
rows = []
with open("/home/hsong/Dropbox/Mu2e/RMCFitting/CDF_Interpolation/fitdata/"+filename,"r") as f:
    reader = csv.reader(f,delimiter=',')
    next(reader)
    #row_count = sum(1 for row in reader)
    for row in reader:
        rows.append(row)

with open("fitdata/fits_kmax_90.9.ssv","w") as newfile:
    writer = csv.writer(newfile, delimiter=" ")
    for r in rows:
        writer.writerow([r[0],r[1],r[2],r[3],r[4]])
newfile.close()
