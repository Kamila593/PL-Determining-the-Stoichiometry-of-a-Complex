import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['text.usetex'] = True


def linear(x, a, b):
    return a*x + b


def linear_function(x, a):
    return a*x


def r_squared(y_true, y_pred):
    residuals = y_true - y_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_true - np.mean(y_true))**2)
    r2 = 1 - (ss_res / ss_tot)
    return r2


print('Metoda stosunków molowych')
C_0ofen = 0.000315
print(f'Początkowe stężenie o-fenantroliny: {C_0ofen} M')
C_0Fe = 0.000315
print(f'Początkowe stężenie soli żelaza(II): {C_0Fe} M')
print('Wielkości stałe: ')
print(f'   Objętość całkowita roztworu i: 4 ml')

df1 = pd.read_csv(f'R3_spectro_1.txt', delimiter=' ', header=None)
df1 = df1.drop(0)
W1 = df1[0].astype(float)
A1 = df1[1].astype(float)
A2 = df1[2].astype(float)
A3 = df1[3].astype(float)
A4 = df1[4].astype(float)
A5 = df1[5].astype(float)
A6 = df1[6].astype(float)

blue_shades = ['#008400', '#0be70a','#10d9b7', '#00e7ff', '#004aff', '#703ec6']


plt.figure(figsize=(12/2.54, 8/2.54))
plt.plot(W1, A1, color = blue_shades[0], label = r'$x_{\mathrm{Fe}}$ = 0')
plt.plot(W1, A2, color = blue_shades[1], label = r'$x_{\mathrm{Fe}}$ = 0,2')
plt.plot(W1, A3, color = blue_shades[2], label = r'$x_{\mathrm{Fe}}$ = 0,4')
plt.plot(W1, A4, color = blue_shades[3], label = r'$x_{\mathrm{Fe}}$ = 0,6')
plt.plot(W1, A5, color = blue_shades[4], label = r'$x_{\mathrm{Fe}}$ = 0,8')
plt.plot(W1, A6, color = blue_shades[5], label = r'$x_{\mathrm{Fe}}$ = 1,0')
plt.xlim(450, 700)
plt.xlabel(r'Długość fali, nm')
plt.ylabel(r'Absorbancja')
plt.legend()
plt.subplots_adjust(bottom=0.2, top=0.9, left = 0.2, right = 0.9)
plt.savefig('ChFizLab_R3_widma1', dpi = 300)
# plt.show()
plt.clf()

A1_510 = df1[(df1[0] == '510.0')][1].values[0]
A2_510 = df1[(df1[0] == '510.0')][2].values[0]
A3_510 = df1[(df1[0] == '510.0')][3].values[0]
A4_510 = df1[(df1[0] == '510.0')][4].values[0]
A5_510 = df1[(df1[0] == '510.0')][5].values[0]
A6_510 = df1[(df1[0] == '510.0')][6].values[0]

# We're creating the plot of absorbance vs. x_Fe
Abs1 = [A1_510, A2_510, A3_510, A4_510, A5_510, A6_510]
xFe = [0, 0.2, 0.4, 0.6, 0.8, 1.0]

xFe_l = [xFe[0], xFe[1]]
xFe_r = xFe[2:5]    # We don't take into consideration the 6th value
Abs1_l = [Abs1[0], Abs1[1]]
Abs1_r = Abs1[2:5]

fitl, covl = curve_fit(linear, xFe_l, Abs1_l)
fitr, covr = curve_fit(linear, xFe_r, Abs1_r)

xfit_l = np.linspace(0, 0.25, 50)
xfit_r = np.linspace(0.13, 1.0, 100)
yfit_l = fitl[0]*xfit_l + fitl[1]
yfit_r = fitr[0]*xfit_r + fitr[1]

print(f'Dopasowana lewa prosta: {fitl[0]:.2f}*x_A + {fitl[1]:.2f}')
print(f'Dopasowana prawa prosta: {fitr[0]:.2f}*x_A + {fitr[1]:.2f}')

Abs1_l_check = fitl[0]*np.array(xFe_l) + fitl[1]
Abs1_r_check = fitr[0]*np.array(xFe_r) + fitr[1]
R2_l = r_squared(Abs1_l, Abs1_l_check)
R2_r = r_squared(Abs1_r, Abs1_r_check)

print(f'Współczynnik dopasowania lewej prostej: {R2_l}')
print(f'Współczynnik dopasowania prawej prostej: {R2_r:.3f}')

x = (fitr[1] - fitl[1])/(fitl[0]-fitr[0])

print(f'Ułamek molowy jonów żelaza(II), dla których proste się przecinają: {x:.3f}')

plt.figure(figsize=(12/2.54, 8/2.54))
plt.scatter(xFe, Abs1, zorder=2, color='#00a5e7')
plt.plot(xfit_l, yfit_l, zorder=1, color='#94e0ff')
plt.plot(xfit_r, yfit_r, zorder=1, color='#94e0ff')
plt.ylim(-0.02, 0.1)
plt.xlabel(r'Ułamek molowy jonów żelaza(II) $x_{\mathrm{Fe}}$')
plt.ylabel(r'Absorbancja')
# plt.legend()
plt.subplots_adjust(bottom=0.2, top=0.9, left = 0.2, right = 0.9)
plt.savefig('ChFizLab_R3_A(xFe)', dpi = 300)
# plt.show()
plt.clf()
