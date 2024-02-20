
import pandas as pd

def Main():
    df = ReadData()
    df = CalculateRadius(df)
    print(df)
    PlotAgeRadius(df)

def ReadData():
    input = "SNR.csv"
    try:
        df = pd.read_csv(input,sep=',',skiprows=0,header=0)
    except IOError:
        print(" cannot open " + path )
        sys.exit()
    #print(df)
    return df


def CalculateRadius(df):
    arcmin = 0.000290888 # radian
    dtheta = (df["size x[arcmin]"]+df["size y[arcmin]"])/4.0*arcmin
    df["radius [pc]"] =  dtheta * df["Distance [kpc]"]*1000.0
    #print(df["radius [pc]"])
    return df


#######################################
# just a preparation for the matplot lib
#######################################
import matplotlib.pyplot as plt
from cycler import cycler
fnameforfig="sans-serif"
fsizeforfig=14
fsizeforlabel=16
plt.rcParams['font.family'] = fnameforfig
plt.rcParams['font.size'] = fsizeforfig
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
# red, blue, green, yellow, sky-blue azure, pink, orange, purple, brown
cmapudc =  cycler(color=["#ff2800","#0041ff" ,"#35a16B","#faf500","#66ccff", "#ff99a0","#ff9900" ,"#9a0079", "#663300"])
plt.rcParams['axes.prop_cycle'] = cmapudc
cmap = ["#ff2800","#0041ff" ,"#35a16B","#faf500","#66ccff", "#ff99a0","#ff9900" ,"#9a0079", "#663300"]

def PlotAgeRadius(df):
    from adjustText import adjust_text
    otputfile = "Age-Radius.png"
    x = df['age [kyr]']*1000
    y = df["radius [pc]"]
    text=df["Name"]

    fig=plt.figure(figsize=(6.4, 5.2),layout="tight")
    ax = fig.add_subplot(1,1,1)
    im =ax.scatter(x,y,marker="o", color=cmap[1])
    texts = [plt.text(x[i], y[i], text[i], ha='center', va='center') for i in range(len(x))]
    adjust_text(texts)
    ax.grid(color='lightgray')
    ax.set_xlabel(r'Age [year]', fontsize=fsizeforlabel)
    ax.set_ylabel(r'Radius [pc]' ,fontsize=fsizeforlabel)
    fig.savefig(otputfile)
Main()
