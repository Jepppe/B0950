import numpy as np
import matplotlib.pyplot as plt
import os, sys, time, math
from scipy.optimize import curve_fit
from lmfit import Model, Parameters, Minimizer, minimize, printfuncs
from lmfit.models import LinearModel
from matplotlib.font_manager import FontProperties
import matplotlib.font_manager

plt.rcParams['font.family']='Times New Roman'

binwid = eval(input("\nBin: "))
#beam = input("\nBeam: ")

name = 'Beam 2 Width ' + str(binwid)

with open(name, 'r') as f:
    lines = f.readlines()
    
with open('0Load', 'w') as f:
    for line in lines:
        f.write(line[1:])
        
with open('0Load', 'r') as f:
    x = np.loadtxt(f, comments = ')', dtype = {'names': ('Time', 'DM', 'SNR', 'Bin'), 'formats': ('f8', 'f8', 'f8', 'int')}, delimiter = ', ')

os.remove('0Load')

y0 = 10 + 0.75*np.log(binwid)/np.log(2)

plt.plot((x['Time'] - 57453.901365074)*24*60, x['SNR'], 'ro')
plt.ylabel('SNR', fontsize = 18)
plt.xlabel('Time (minutes)', fontsize = 18)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)


print(y0)

plt.axis([0, 350, y0, 60])
plt.show()
    
rownum = 0
rowdel = [0]

for row in x:
    if x['Time'][rownum] > 57453.936:
        rowdel = np.append(rowdel, rownum)
    rownum += 1
    
rowdel = np.delete(rowdel, 0)
x = np.delete(x, rowdel)

plt.plot((x['Time'] - x[0]['Time'])*24*60+5, x['SNR'], 'ro')
plt.ylabel('SNR', fontsize = 18)
plt.xlabel('Time (minutes)', fontsize = 18)
plt.axis([0, 55, y0, 40])
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.show()

rownum = 0
for row in x:
    rownum += 1
    
nrows = rownum

print(nrows)

numpul = 1
rownum = 0

while rownum < nrows - 1:
    if not x[rownum + 1]['Time'] - x[rownum]['Time'] < 0.25306/(2*60*60*24):
        numpul += 1
    rownum += 1
    
print(numpul)

x = np.sort(x, order = 'Time')

timeint = x[-1]['Time'] - x[0]['Time']

periods = timeint/(0.25306/(60*60*24))

print(numpul/periods)

#Count the number of pulses in individual pulse windows

numpul = 1
rownum = 0

while rownum < nrows - 1:
    if not x[rownum + 1]['Time'] - x[rownum]['Time'] < 0.25306/(2*60*60*24):
        numpul += 1
    rownum += 1
    
print(rownum+1)
print(numpul)

'''This section deals with the integrated pulse intensities'''

#Create an array to store the data for the integrated pulse intensities

intsnr = np.zeros(numpul, dtype = {'names': ('Pulse', 'Time', 'Int SNR'), 'formats': ('int', 'f8', 'f8')})

#Loop over the array to fill up the pulse numbers

pulse = 0

while pulse < numpul:
    intsnr[pulse]['Pulse'] = pulse + 1
    pulse += 1

#Work out the integrated SNRs of pulses
    
rownum = 0
pulse = 0

intsnr[0]['Time'] = x[0]['Time']
intsnr[0]['Int SNR'] = x[0]['SNR']

while rownum < nrows - 1:
    if not x[rownum + 1]['Time'] - x[rownum]['Time'] <= 0.253/(2*60*60*24):
        pulse += 1
        intsnr[pulse]['Time'] = x[rownum + 1]['Time']
        intsnr[pulse]['Int SNR'] = intsnr[pulse]['Int SNR'] + x[rownum + 1]['SNR']
    else:
        intsnr[pulse]['Int SNR'] = intsnr[pulse]['Int SNR'] + x[rownum + 1]['SNR']
    rownum += 1
    
print(np.sum(x['SNR']))
print(np.sum(intsnr['Int SNR']))

plt.plot(intsnr['Time'], intsnr['Int SNR'], 'ro')
plt.show()
     
#Choose the number of bins and size of bins and generate an array of appropriate size

numbin = eval(input('Choose the number of points (default 30): '))
binsize = eval(input('Choose the length between points (default 2): '))

#numbin = 30
#binsize = 2

binarray = np.zeros(numbin, dtype = np.int)

#Loop which counts number of pulses in each bin

rownum = 0
while rownum < numpul:
    bin = 0
    while bin < numbin:
        if 10 + 0.75*np.log(binwid)/np.log(2) + binsize*bin <= intsnr[rownum]['Int SNR']:
            binarray[bin] += 1
        bin += 1
    rownum += 1

binx = np.zeros(numbin)

bin = 0
while bin < numbin:
    binx[bin] = 10 + 0.75*np.log(binwid)/np.log(2) + binsize*bin
    bin += 1
    
#Plot the histogram

plt.bar(binx, binarray, width = binsize)
plt.ylabel('Number of pulses')
plt.xlabel('Scaled SNR')
plt.title('Scaled SNR of pulses')
plt.xticks(np.arange(10 + 0.75*np.log(binwid)/np.log(2), binsize*numbin, binsize))

plt.show()

#Remove infinite values

logbin = np.log10(binarray)
logx = np.log10(binx)
rowdel = [0]
r = 0

for row in logbin:
    if math.isfinite(logbin[r]) == False:
        rowdel = np.append(rowdel, r)
    r += 1
        
rowdel = np.delete(rowdel, 0)        
logbin = np.delete(logbin, rowdel)
logx = np.delete(logx, rowdel)
erx = np.delete(binx, rowdel)

#Plot the log-log plot

plt.plot(logx, logbin, 'ro')
plt.ylabel('log(Number of pulses)')
plt.xlabel('log(Scaled SNR)')
plt.title('Log-log plot of Scaled SNR of pulses')

plt.show()

#Plot the log-log plot with a line of best fit

decision = eval(input('\nEnter 1 for a single power law exponent or 2 for a piecewise power law: '))

width = eval(input('\nEnter the thickness of the lines on the plot: '))

if decision == 1:

    fitlogbin = np.copy(logbin)
    fitlogx = np.copy(logx)

    rowdel = [0]
    r = 0

    for row in logbin:
        if fitlogbin[r] <= 1: 
            rowdel = np.append(rowdel, r)
        r += 1
        
    rowdel = np.delete(rowdel, 0)
    fitlogx = np.delete(fitlogx, rowdel)
    fitlogbin = np.delete(fitlogbin, rowdel)
    fiterx1 = np.delete(erx, rowdel)
     
    fitval = np.polyfit(fitlogx, fitlogbin, 1)
    a = np.linspace(fitlogx[0] - 0.02, fitlogx[-1] + 0.02, 100)

    #Use lmfit instead, for the plot have run the program twice and manually input the values to be used from what was found using the fitting function

    mod = LinearModel()
    result = mod.fit(fitlogbin, x = fitlogx, slope = -3, intercept = 10)
    print(result.fit_report())

    '''plt.plot(logx, logbin, 'ro', b, b*-3.62212153 + 6.96592205, 'k-', a, a*-6.39557577 + 10.6634723, 'k-')
    plt.errorbar(logx, logbin, xerr = 1/erx, yerr = 0, fmt = 'ro')
    plt.ylabel('log(Cumulative Number of pulses > Scaled SNR)')
    plt.xlabel('log(Scaled SNR)')

    plt.show()'''

    print(fitval[0])
    print(fitval[1])

    #Try writing the function

    x = fitlogbin

    data = fitlogx

    error = 1/fiterx1

    def residual(pars, x, data=None):
            a=pars['intercept'].value
            b=pars['slope'].value
            model = -a/b + ((1/b)*x)
            resids = model - data
            weighted = np.sqrt(resids ** 2 / error ** 2)
            return weighted

    params=Parameters()
    params.add('intercept', value=10.0)
    params.add('slope', value=-3.0)

    mi=minimize(residual, params, args=(x, data))
        
    print(printfuncs.report_fit(mi.params, min_correl=0.5))

    slope = eval(input("\nSlope: "))
    sloperr = eval(input("Slope error: "))
    int = eval(input("Intercept: "))

    plt.plot(logx, logbin, 'ro')
    plt.plot(a, a*slope + int, 'k-', linewidth = width, label = r'$ \alpha = %.2f  \pm  %.2f$' %(slope, sloperr))
    plt.errorbar(logx, logbin, xerr = 1/erx, yerr = 0, fmt = 'ro')
    plt.ylabel('log(Cumulative number of pulses > SNR)', fontsize = 18)
    plt.xlabel('log(SNR)', fontsize = 18)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(frameon = False, fontsize = 18)

    plt.show()
elif decision == 2:

    fitlogbin = np.copy(logbin)
    fitlogx = np.copy(logx)
    fitlogx1 = np.copy(logx)
    fit2 = np.copy(fitlogbin)
    
    #Choose the cutoff point for the piecewise power law
    
    cutoff = eval(input('\nChoose a point on the y axis to use to differentiate between the regimes (default 2.3): '))

    rowdel = [0]
    r = 0

    for row in logbin:
        if fitlogbin[r] > cutoff or fitlogbin[r] <= 1: 
            rowdel = np.append(rowdel, r)
        r += 1
        
    rowdel = np.delete(rowdel, 0)
    fitlogx = np.delete(fitlogx, rowdel)
    fitlogbin = np.delete(fitlogbin, rowdel)
    fiterx1 = np.delete(erx, rowdel)
     
    fitval = np.polyfit(fitlogx, fitlogbin, 1)
    a = np.linspace(fitlogx[0] - 0.02, fitlogx[-1] + 0.02, 100)

    rowdel = [0]
    r = 0

    for row in logbin:
        if fit2[r] < cutoff: 
            rowdel = np.append(rowdel, r)
        r += 1
        
    rowdel = np.delete(rowdel, 0)
    fitlogx1 = np.delete(fitlogx1, rowdel)
    fit2 = np.delete(fit2, rowdel)
    fiterx2 = np.delete(erx, rowdel)
     
    fitval1 = np.polyfit(fitlogx1, fit2, 1)
    b = np.linspace(fitlogx1[0] - 0.02, fitlogx1[-1] + 0.02, 100)

    #Use lmfit instead, for the plot have run the program twice and manually input the values to be used from what was found using the fitting function

    mod = LinearModel()
    result = mod.fit(fitlogbin, x = fitlogx, slope = -3, intercept = 10)
    #print(result.fit_report())

    result2 = mod.fit(fit2, x = fitlogx1, slope = -3, intercept = 10)
    #print(result2.fit_report())

    '''plt.plot(logx, logbin, 'ro', b, b*-3.62212153 + 6.96592205, 'k-', a, a*-6.39557577 + 10.6634723, 'k-')
    plt.errorbar(logx, logbin, xerr = 1/erx, yerr = 0, fmt = 'ro')
    plt.ylabel('log(Cumulative Number of pulses > Scaled SNR)')
    plt.xlabel('log(Scaled SNR)')

    plt.show()

    print(fitval[0])
    print(fitval[1])'''

    #Try writing the function

    x = fitlogbin

    data = fitlogx

    error = 1/fiterx1

    def residual(pars, x, data=None):
            a=pars['intercept'].value
            b=pars['slope'].value
            model = -a/b + ((1/b)*x)
            resids = model - data
            weighted = np.sqrt(resids ** 2 / error ** 2)
            return weighted

    params=Parameters()
    params.add('intercept', value=10.0)
    params.add('slope', value=-3.0)

    mi=minimize(residual, params, args=(x, data))
        
    print(printfuncs.report_fit(mi.params, min_correl=0.5))

    x = fit2

    data = fitlogx1

    error = 1/fiterx2

    def residual(pars, x, data=None):
            a=pars['intercept'].value
            b=pars['slope'].value
            model = -a/b + ((1/b)*x)
            resids = model - data
            weighted = np.sqrt(resids ** 2 / error ** 2)
            return weighted

    params=Parameters()
    params.add('intercept', value=10.0)
    params.add('slope', value=-3.0)

    mi=minimize(residual, params, args=(x, data))
        
    print(printfuncs.report_fit(mi.params, min_correl=0.5))

    stslope = eval(input("\nSteep slope: "))
    stsloperr = eval(input("Steep slope error: "))
    stint = eval(input("Steep intercept: "))
    shslope = eval(input("\nShallow slope: "))
    shsloperr = eval(input("Shallow slope error: "))
    shint = eval(input("Shallow intercept: "))

    plt.plot(logx, logbin, 'ro')
    plt.plot(b, b*shslope + shint, 'k-', linewidth = width, label = r'$ \alpha = %.2f \pm %.2f$' %(shslope, shsloperr))
    plt.plot(a, a*stslope + stint, 'k--', linewidth = width, label = r'$ \alpha = %.2f \pm %.2f$' %(stslope, stsloperr))
    plt.errorbar(logx, logbin, xerr = 1/erx, yerr = 0, fmt = 'ro')
    plt.ylabel('log(Cumulative number of pulses > SNR)', fontsize = 18)
    plt.xlabel('log(SNR)', fontsize = 18)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(frameon = False, fontsize = 18)

    plt.show()
else:
    print('\nIncorrect entry.')
