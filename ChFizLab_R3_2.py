import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import ChFizLab_R3_1 as R1

plt.rcParams['text.usetex'] = True

print('Metoda stosunków molowych')

# Volumes and concentrations of initial solutions
C_0ofen = 0.0044
print(f'Początkowe stężenie o-fenantroliny: {C_0ofen} M')
C_0Fe = 0.000315
print(f'Początkowe stężenie soli żelaza(II): {C_0Fe} M')
print('Wielkości stałe: ')
V_H = 1500
print(f'Objętość roztworu kwasu: {V_H} mikro moli')
V_Fe = 1500
print(f'   Objętość roztworu soli żelaza(II): {V_Fe} mikro moli')
V_ofen = [0, 25, 50, 100, 150, 200, 250, 300, 400]
n_Fe = V_Fe*C_0Fe
print(f'   Liczba mikro moli Fe2+: {n_Fe:.4f}')


# Creating columns with data for wavelengths and absorbance
df2 = pd.read_csv(f'R3_spectro_2.txt', delimiter=' ', header=None)
df2 = df2.drop(0)
W2 = df2[0].astype(float)
A0 = df2[1].astype(float)
A25 = df2[2].astype(float)
A50 = df2[3].astype(float)
A100 = df2[4].astype(float)
A150 = df2[5].astype(float)
A200 = df2[6].astype(float)
A250 = df2[7].astype(float)
A300 = df2[8].astype(float)
A400 = df2[9].astype(float)

plot_colours = ['#002cff', '#a54af5', '#ed00ff', '#d609ae', '#ff0000', '#ff8e00', '#ffc600', '#ffef00', '#c0ff00']

# Plot with spectras for all the solutions
plt.figure(figsize=(12/2.54, 8/2.54))
plt.plot(W2, A0, color=plot_colours[0], label=r'0 $\mathrm{\mu l}$')
plt.plot(W2, A25, color=plot_colours[1], label=r'25 $\mathrm{\mu l}$')
plt.plot(W2, A50, color=plot_colours[2], label=r'50 $\mathrm{\mu l}$')
plt.plot(W2, A100, color=plot_colours[3], label=r'100 $\mathrm{\mu l}$')
plt.plot(W2, A150, color=plot_colours[4], label=r'150 $\mathrm{\mu l}$')
plt.plot(W2, A200, color=plot_colours[5], label=r'200 $\mathrm{\mu l}$')
plt.plot(W2, A250, color=plot_colours[6], label=r'250 $\mathrm{\mu l}$')
plt.plot(W2, A300, color=plot_colours[7], label=r'300 $\mathrm{\mu l}$')
plt.plot(W2, A400, color=plot_colours[8], label=r'400 $\mathrm{\mu l}$')
plt.xlim(370, 590)
plt.ylim(-0.05, 0.65)
plt.xlabel(r'Długość fali, nm')
plt.ylabel(r'Absorbancja')
# plt.legend(ncol=2, bbox_to_anchor=(1.5, 1))
plt.legend(bbox_to_anchor=(1.0, 0.93))
plt.subplots_adjust(bottom=0.2, top=0.9, left = 0.1, right = 0.8)
plt.savefig('ChFizLab_R3_widma2', dpi = 300)
#plt.show()
plt.clf()

# We're looking for the absorbance values for 510.0 nm (max value)
A0_510 = df2[(df2[0] == '510.0')][1].values[0]
A25_510 = df2[(df2[0] == '510.0')][2].values[0]
A50_510 = df2[(df2[0] == '510.0')][3].values[0]
A100_510 = df2[(df2[0] == '510.0')][4].values[0]
A150_510 = df2[(df2[0] == '510.0')][5].values[0]
A200_510 = df2[(df2[0] == '510.0')][6].values[0]
A250_510 = df2[(df2[0] == '510.0')][7].values[0]
A300_510 = df2[(df2[0] == '510.0')][8].values[0]
A400_510 = df2[(df2[0] == '510.0')][9].values[0]
# print(f'Absorbance values for wavelength 510.0 nm: {A0_510}, {A25_510}, {A150_510}')

Abs2 = [A0_510, A25_510, A50_510, A100_510, A150_510, A200_510, A250_510, A300_510, A400_510]
Abs2_reversed = Abs2[::-1]

df2r = pd.DataFrame({0: Abs2})
df2r[1] = V_ofen                    # volume of o-phenantroline solution
df2r[2] = df2r[1] + V_Fe + V_H      # total volume
df2r[3] = df2r[1]*C_0ofen           # n_ofen
df2r[4] = df2r[3]/(n_Fe)            # n_ofen / n_Fe
df2r[5] = df2r[4][::-1].values      # n_ofen / n_Fe reversed order
df2r[6] = n_Fe/df2r[3]              # n_Fe / n_ofen
df2r[7] = ((1/3)*df2r[3])/df2r[2]   # ferroin concentration

df2r.columns = ['A', 'V_ofen', 'V_c', 'n_ofen', 'n_ofen/n_Fe', 'n_ofen/n_Fe reversed', 'n_Fe/n_ofen', 'C_ferroin']

df2r['n_ofen/n_Fe'] = df2r['n_ofen/n_Fe'].round(2)
df2r['n_ofen/n_Fe reversed'] = df2r['n_ofen/n_Fe reversed'].round(2)

print(df2r)


# We're fitting a model to the data
nAnB_l = df2r['n_ofen/n_Fe'][0:8]
Abs2_l = Abs2[0:8]
fitl, covl = curve_fit(R1.linear_function, nAnB_l, Abs2_l)

xfit_l = np.linspace(0, 5.5, 50)
xfit_r = np.linspace(3.0, 5.5, 100)
yfit_l = fitl[0]*xfit_l
yfit_r = 0*xfit_r + Abs2[8]

x = Abs2[8] / fitl[0]

Abs2_l_pred = np.array(nAnB_l)*fitl[0]
R2l1 = R1.r_squared(Abs2_l, Abs2_l_pred)

print(f'Dopasowana lewa prosta: {fitl[0]:.3f}*x')       # 0.092*x
print(f'Proste przecinają się dla x = {x:.2f}')         # x = 5.33
print(f'Współczynnik dopasowania prostej: {R2l1:.3f}')  # R2 = 0.772

plt.figure(figsize=(12/2.54, 8/2.54))
plt.scatter(df2r['n_ofen/n_Fe'], df2r['A'], color='#cd09ca', zorder=2)
plt.plot(xfit_l, yfit_l, color='#ffa0fe', zorder=1)
plt.plot(xfit_r, yfit_r, color='#ffa0fe', zorder=1)
plt.xlabel(r'Stosunek molowy $n_{\mathrm{fen}} / n_{\mathrm{Fe}}$')
plt.ylabel(r'Absorbancja')
plt.xlim(-0.1, 5.5)
plt.ylim(-0.05, 0.55)
plt.subplots_adjust(bottom=0.2, top=0.9, left = 0.2, right = 0.9)
plt.savefig('ChFizLab_R3_A(nAnB)', dpi = 300)
# plt.show()
plt.clf()


# Second attempt on fitting the model
# nAnB_l2 = [df2r['n_ofen/n_Fe'][0], df2r['n_ofen/n_Fe'][7]]
# Abs2_l2 = [Abs2[0], Abs2[7]]
# fitl2, covl2 = curve_fit(R1.linear_function, nAnB_l2, Abs2_l2)

nAnB_l2 = [df2r['n_ofen/n_Fe'][i] for i in [5, 6, 7]]
Abs2_l2 = [Abs2[i] for i in [5, 6, 7]]
fitl2, covl2 = curve_fit(R1.linear, nAnB_l2, Abs2_l2)

xfit_l2 = np.linspace(0, 5.5, 50)
xfit_r2 = np.linspace(3.0, 5.5, 100)
yfit_l2 = fitl2[0]*xfit_l2 + fitl2[1]
yfit_r2 = 0*xfit_r2 + Abs2[8]

x2 = (Abs2[8] - fitl2[1]) / fitl2[0]

Abs2_l2_pred = np.array(nAnB_l2)*fitl2[0] + fitl2[1]
R2l2 = R1.r_squared(Abs2_l2, Abs2_l2_pred)

print(f'Dopasowana lewa prosta: {fitl2[0]:.3f}*x + {fitl2[1]:.3f}')      # 0.128*x
print(f'Proste przecinają się dla x = {x2:.2f}')        # x = 3.83
print(f'Współczynnik dopasowania prostej: {R2l2:.3f}')  # R2 = 0.963

plt.figure(figsize=(12/2.54, 8/2.54))
plt.scatter(nAnB_l2, Abs2_l2, color='#cd09ca', zorder=2)
plt.scatter(df2r['n_ofen/n_Fe'][8], Abs2[8], color='#cd09ca', zorder=2)
plt.plot(xfit_l2, yfit_l2, color='#ffa0fe', zorder=1)
plt.plot(xfit_r2, yfit_r2, color='#ffa0fe', zorder=1)
plt.xlabel(r'Stosunek molowy $n_{\mathrm{fen}} / n_{\mathrm{Fe}}$')
plt.ylabel(r'Absorbancja')
plt.xlim(-0.1, 5.5)
plt.ylim(-0.05, 0.55)
plt.subplots_adjust(bottom=0.2, top=0.9, left = 0.2, right = 0.9)
plt.savefig('ChFizLab_R3_A(nAnB)2', dpi = 300)
# plt.show()
