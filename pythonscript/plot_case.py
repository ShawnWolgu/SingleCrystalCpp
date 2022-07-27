import os, sys
from posix import getcwd
import numpy as np
import pandas as pd
from snplot.snplot import plotdata

def getplot_database(folders,filename,xcol,ycol):
    database = []
    for ifol in folders:
        folderp = os.path.join(os.getcwd(),ifol)
        filep = os.path.join(folderp, filename)
        df = pd.read_csv(filep)
        dfcsv = pd.concat([df.iloc[:,xcol],df.iloc[:,ycol]],axis=1)
        database.append({'linedata':dfcsv,'label':ifol})
    return database

folders = sys.argv[1:]
plotdata(getplot_database(folders,'stress_grain.csv',2,8),'stress')


