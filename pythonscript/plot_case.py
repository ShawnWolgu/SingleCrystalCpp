import os, sys
from posix import getcwd
import numpy as np
import pandas as pd
from snplot.snplot import plotdata

def getplot_database(folders,filename,xcol,ycol):
    database = []
    xlim = [0,0]
    ylim = [0,0]
    for ind, ifol in enumerate(folders):
        folderp = os.path.join(os.getcwd(),ifol)
        filep = os.path.join(folderp, filename)
        df = pd.read_csv(filep)
        interval = max(round(len(df.iloc[:,xcol])/200),1)
        dfcsv = pd.concat([df.iloc[::interval,xcol],df.iloc[::interval,ycol]],axis=1)
        xlim[0] = min(dfcsv.iloc[0,0],xlim[0])
        xlim[1] = max(dfcsv.iloc[-1,0],xlim[1])
        ylim[0] = min(dfcsv.iloc[0,1],ylim[0])
        ylim[1] = max(dfcsv.iloc[:,1].max(),ylim[1])
        if ind != 10000:
            database.append({'linedata':dfcsv,'label':ifol})
        else:
            database.append({'dashdata':dfcsv,'label':ifol})
    return database, xlim, ylim

name = sys.argv[1]
folders = sys.argv[2:]
datab,xlim,ylim = getplot_database(folders,'stress_grain.csv',3,9)
plotdata(datab,name,xlim=(xlim[0],xlim[1]),ylim=(ylim[0],1.2*ylim[1]))
datab,xlim,ylim = getplot_database(folders,'stress_step.csv',0,9)
plotdata(datab,name + '_step',xlim=(xlim[0],xlim[1]),ylim=(ylim[0],1.2*ylim[1]))


