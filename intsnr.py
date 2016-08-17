import numpy as np
import matplotlib.pyplot as plt
import os, sys, time, math
from scipy.optimize import curve_fit
from lmfit import Model, Parameters, Minimizer, minimize, printfuncs
from lmfit.models import LinearModel
from matplotlib.font_manager import FontProperties
import matplotlib.font_manager

plt.rcParams['font.family']='Times New Roman'

#Define the guassian in order to scale the data later
    
def gaussian(x, A, mu, sigma):
    return A*np.exp(-(x - mu)**2./(2.*sigma**2.))

#Decide which bin width to look at, open the right file and use loadtxt to create the array
    
bin = input("\nBin: ")

name = str(bin) + 'Combined'

with open(name, 'r') as f:
    lines = f.readlines()
    
with open('0Load', 'w') as f:
    for line in lines:
        f.write(line[1:])
        
with open('0Load', 'r') as f:
    x = np.loadtxt(f, comments = ')', dtype = {'names': ('Time', 'DM', 'SNR', 'Bin'), 'formats': ('f8', 'f8', 'f8', 'int')}, delimiter = ', ')

#Create two copies of x; y and z. y has it's time changed to LST and has the SNR scaled, z only has its SNR scaled    

y = np.copy(x)
z = np.copy(x)

x = np.sort(x, order = 'Time')
z = np.sort(z, order = 'Time')
y = np.sort(y, order = 'Time')

#Loop which looks at all the data and decides if it is relevant. A loop works out the MJD for a sidereal time of 9h53m
    
rownum = 0
a1 = [0]
failures = 0

for row in x:
    MJD = x[rownum]['Time']

    MJD0 = math.floor(MJD)
    
    if MJD0 <= 57030: #Excluding these because it isn't the right MJD (16767)
        a1 = np.append(a1, rownum)
        rownum += 1
        continue

    H = 24*(MJD - MJD0)

    JD0 = MJD0 + 2400000.5

    D0 = JD0 - 2451545.0

    GMST = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*H

    LST = GMST%24

    LSTHour = math.floor(LST)

    LSTMinutes = (LST - LSTHour)*60

    H = 27
    r = 0

    while (H/24) > 1.0 or (H/24) < 0.0:
        H = ((9.897859735-0.00005 + 24*r) - 6.697374558 - 0.06570982441908*D0)/1.00273790935   #was 9.895
        r += 1
    
    Decimal = H/24
    
    Centre = MJD0 + Decimal

    #Entries that do not correspond to the right time are removed
 
    if not -5/(60*24) < x[rownum]['Time'] - Centre < 5/(60*24):
        a1 = np.append(a1, rownum)
        
    if (0.9972695663 - 10/(60*24) < x[rownum]['Time'] - Centre < 0.9972695663 + 10/(60*24)) or (-0.9972695663 + 10/(60*24) < x[rownum]['Time'] - Centre < -0.9972695663 - 10/(60*24)):
        failures += 1
    
    rownum += 1
    
a1 = np.delete(a1, 0)

x = np.delete(x, a1)
y = np.delete(y, a1)
z = np.delete(z, a1)

#The times are converted from MJD to LST

rownum = 0

for row in y:
    MJD = y[rownum]['Time']

    MJD0 = math.floor(MJD)
    
    H = 24*(MJD - MJD0)

    JD0 = MJD0 + 2400000.5

    D0 = JD0 - 2451545.0

    GMST = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*H

    LST = GMST%24

    y[rownum]['Time'] = LST
    
    rownum += 1

#Plot measured SNR vs LST
    
plt.plot(y['Time'], z['SNR'], 'ro')
plt.ylabel('SNR')
plt.xlabel('LST')
plt.title('SNR vs LST')

plt.show()

#Perform the scaling

rownum = 0
for row in y:
    z[rownum]['SNR'] = y[rownum]['SNR']/gaussian(y[rownum]['Time'], 1, 9.897859735-0.00005, 0.09)
    y[rownum]['SNR'] = y[rownum]['SNR']/gaussian(y[rownum]['Time'], 1, 9.897859735-0.00005, 0.09)
    rownum += 1

print(10/gaussian(9.897859735-0.00005 + 5/60, 1, 9.897859735-0.00005, 0.09))

nrows = rownum
    
#Plot scaled SNR vs LST
    
'''plt.plot(y['Time'], z['SNR'], 'ro')
plt.ylabel('Normalized SNR')
plt.xlabel('LST')
plt.title('Normalized SNR vs LST')

plt.show()'''

#Plot scaled SNR vs MJD

plt.plot(z['Time'], z['SNR'], 'ro')

plt.show()

#Count the number of pulses in individual pulse windows

numpul = 1
rownum = 0

while rownum < nrows - 1:
    if not z[rownum + 1]['Time'] - z[rownum]['Time'] < 0.25306/(2*60*60*24):
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

intsnr[0]['Time'] = z[0]['Time']
intsnr[0]['Int SNR'] = z[0]['SNR']

while rownum < nrows - 1:
    if not z[rownum + 1]['Time'] - z[rownum]['Time'] <= 0.253/(2*60*60*24):
        pulse += 1
        intsnr[pulse]['Time'] = z[rownum + 1]['Time']
        intsnr[pulse]['Int SNR'] = intsnr[pulse]['Int SNR'] + z[rownum + 1]['SNR']
    else:
        intsnr[pulse]['Int SNR'] = intsnr[pulse]['Int SNR'] + z[rownum + 1]['SNR']
    rownum += 1
    
print(np.sum(z['SNR']))
print(np.sum(intsnr['Int SNR']))

'''plt.plot(intsnr['Time'], intsnr['Int SNR'], 'ro')
plt.ylabel('SNR')
plt.xlabel('Number of pulse')
plt.show()'''
     
#Choose the number of bins and size of bins and generate an array of appropriate size

numbin = eval(input('Choose the number of points (default 30): '))
binsize = eval(input('Choose the length between points (default 3): '))

#numbin = 30
#binsize = 3

binarray = np.zeros(numbin, dtype = np.int)

#Loop which counts number of pulses in each bin

rownum = 0
while rownum < numpul:
    bin = 0
    while bin < numbin:
        if 10 + binsize*bin <= intsnr[rownum]['Int SNR']:
            binarray[bin] += 1
        bin += 1
    rownum += 1

binx = np.zeros(numbin)

bin = 0
while bin < numbin:
    binx[bin] = 10 + binsize*bin
    bin += 1
    
#Plot the histogram

'''plt.bar(binx, binarray, width = binsize)
plt.ylabel('Number of pulses')
plt.xlabel('Scaled SNR')
plt.title('Scaled SNR of pulses')
plt.xticks(np.arange(10, binsize*numbin, binsize))

plt.show()'''

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
        if fitlogx[r] <= np.log10(10/gaussian(9.897859735-0.00005 + 5/60, 1, 9.897859735-0.00005, 0.09)) or fitlogbin[r] < 1: 
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
    plt.xticks(fontsize = 15)
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
        if fitlogx[r] <= np.log10(10/gaussian(9.897859735-0.00005 + 5/60, 1, 9.897859735-0.00005, 0.09)) or fitlogbin[r] > cutoff or fitlogbin[r] < 1: 
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
        if fitlogx1[r] <= np.log10(10/gaussian(9.897859735-0.00005 + 5/60, 1, 9.897859735-0.00005, 0.09)) or fit2[r] < cutoff: 
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