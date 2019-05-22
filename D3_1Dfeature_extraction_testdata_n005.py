####################################################################
# Author: Shan Qu, Delft University of Technology                ###
# First created: May 2019                                        ###
# This code is propriatary under the Delphi Research Consortium  ###
#                                                                ###
# product: ML workflow for event detection                       ###
# generate 1D features
####################################################################


import pandas as pd
import numpy as np
from sklearn import svm, datasets
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.utils import shuffle
import librosa
import scipy
import math
import cmath
from statsmodels import robust

def cal_energy(signal):

    energy = np.sum(signal ** 2) / np.float64(len(signal))

    return energy


def cal_energyentropy(signal, num_blocks=10):
    """Computes entropy of energy"""

    tot = np.sum(signal ** 2)

    length = len(signal)

    subwinlength = int(np.floor(length / num_blocks))

    if length != subwinlength * num_blocks:
            signal = signal[0:subwinlength * num_blocks]

    subwindows = signal.reshape(subwinlength, num_blocks, order='F').copy()

    # Compute normalized sub-frame energies:
    s = np.sum(subwindows ** 2, axis=0) / (tot + 0.00000001)

    # Compute entropy of the normalized sub-frame energies:
    entropy = -np.sum(s * np.log2(s + 0.00000001))

    return entropy


def cal_SpectralCentroidAndSpread(X, fs):
    """Computes spectral centroid of frame (given abs(FFT))"""

    ind = (np.arange(1, X.shape[0] + 1)) * (fs/(2.0 * X.shape[0]))

    Xt = X.copy()
    Xt = Xt / Xt.max()
    NUM = np.sum(ind * Xt)
    DEN = np.sum(Xt) + 0.00000001

    # Centroid:
    C = (NUM / DEN)

    # Spread:
    S = np.sqrt(np.sum(((ind - C) ** 2) * Xt) / DEN)

    # Normalize:
    C = C / (fs / 2.0)
    S = S / (fs / 2.0)

    return (C, S)


def cal_SpectralEntropy(X, numOfShortBlocks=10):
    """Computes the spectral entropy"""
    L = X.shape[0]                         # number of frame samples
    Eol = np.sum(X ** 2)            # total spectral energy

    subWinLength = int(np.floor(L / numOfShortBlocks))   # length of sub-frame
    if L != subWinLength * numOfShortBlocks:
        X = X[0:subWinLength * numOfShortBlocks]

    subWindows = X.reshape(subWinLength, numOfShortBlocks, order='F').copy()  # define sub-frames (using matrix reshape)
    s = np.sum(subWindows ** 2, axis=0) / (Eol + 0.00000001)                      # compute spectral sub-energies
    En = -np.sum(s*np.log2(s + 0.00000001))                                    # compute spectral entropy

    return En


def cal_SpectralRollOff(X, c, fs):
    """Computes spectral roll-off"""
    totalEnergy = np.sum(X ** 2)
    fftLength = X.shape[0]
    Thres = c*totalEnergy
    # Ffind the spectral rolloff as the frequency position where the respective spectral energy is equal to c*totalEnergy
    CumSum = np.cumsum(X ** 2) + 0.00000001
    [a, ] = np.nonzero(CumSum > Thres)
    if a.shape[0] > 0:
        mC = np.float64(a[0]) / (float(fftLength))
    else:
        mC = 0.0
    return (mC)


def cal_SpectralFlatness(powerSpectrum):
    geometricMean = np.float64(0)
    arithmeticMean = 0
    for i in range(len(powerSpectrum)):
        y = np.float64(powerSpectrum[i])
        geometricMean += np.float64(math.log(y))
        arithmeticMean += y
    geometricMean /= np.float64(len(powerSpectrum))
    geometricMean = math.exp(geometricMean)
    arithmeticMean /= float(len(powerSpectrum))
    spectralFlatness = geometricMean / arithmeticMean

    return spectralFlatness


# load and concatenate data and noise
mat_contents = scipy.io.loadmat('Data/data2_n005.mat')
full_data = mat_contents['data']
mat_contents = scipy.io.loadmat('Data/label2_n005.mat')
label_full_data = mat_contents['label']

Nt = full_data.shape[0]
N = full_data.shape[1]
dt = 0.0005
Fs = 1. / dt
df = 1. / (Nt * dt)

# fourier transform
FULL_DATA = np.zeros_like(full_data)
for i in range(N):
    FULL_DATA[:, i] = abs(scipy.fftpack.fft(full_data[:, i]))
    FULL_DATA[:, i] /= Nt

print "1"

# different features
N_feature = 63
feature = np.zeros((N, N_feature))

feature[:, 0] = np.mean(full_data, axis=0)

feature[:, 1] = np.median(full_data, axis=0)

feature[:, 2] = np.std(full_data, axis=0)

feature[:, 3] = robust.mad(full_data, axis=0)

feature[:, 4] = np.percentile(full_data, q=25, axis=0)

feature[:, 5] = np.percentile(full_data, q=75, axis=0)

feature[:, 6] = feature[:, 6] - feature[:, 5]

feature[:, 7] = scipy.stats.skew(full_data, axis=0)

feature[:, 8] = scipy.stats.kurtosis(full_data, axis=0)

for i in range(N):
    feature[i, 9] = librosa.feature.zero_crossing_rate(full_data[:, i])[0, 0]

print "2"

for i in range(N):
    feature[i, 10] = cal_energyentropy(full_data[:, i])
    feature[i, 11] = cal_energy(full_data[:, i])

print "3"

for i in range(N):
    feature[i, 12 : 25] = np.mean(librosa.feature.mfcc(full_data[:, i], Fs, n_mfcc=13).T, axis=0)

    feature[i, 25] = np.max(FULL_DATA[:, i])

    [feature[i, 26], feature[i, 27]] = cal_SpectralCentroidAndSpread(FULL_DATA[:, i], Fs)

    feature[i, 28] = cal_SpectralEntropy(FULL_DATA[:, i])

    feature[i, 29] = np.mean(librosa.feature.spectral_rolloff(full_data[:, i], Fs, roll_percent=0.90).T, axis=0)

    feature[i, 30] = np.mean(librosa.feature.rmse(y=full_data[:, i], frame_length=Nt).T, axis=0)

    feature[i, 31] = np.mean(librosa.feature.spectral_bandwidth(full_data[:, i], Fs).T, axis=0)

    feature[i, 32 : 36] = np.mean(librosa.feature.poly_features(full_data[:, i], Fs, order=3).T, axis=0)

    feature[i, 36 : 48] = np.mean(librosa.feature.chroma_stft(full_data[:, i], Fs).T, axis=0)

    feature[i, 48] = np.std(feature[i, 36 : 48])

    feature[i, 49 : 56] = np.mean(librosa.feature.spectral_contrast(full_data[:, i], Fs, fmin=5, n_bands=6).T, axis=0)

    feature[i, 56] = cal_SpectralFlatness(FULL_DATA[:, i])

    feature[i, 57 : 63] = np.mean(librosa.feature.tonnetz(full_data[:, i], Fs, chroma=librosa.feature.chroma_stft(full_data[:, i], Fs)).T, axis=0)


# plt.figure()
# plt.plot(feature[:, 24])
# plt.figure()
# plt.plot(feature[:, 25])
# plt.figure()
# plt.plot(feature[:, 54])
# plt.figure()
# plt.plot(feature[:, 55])
# plt.figure()
# plt.plot(feature[:, 56])
# plt.figure()
# plt.plot(feature[:, 57])
# plt.figure()
# plt.plot(feature[:, 58])
# plt.figure()
# plt.plot(feature[:, 59])
# plt.figure()
# plt.plot(feature[:, 51])
# plt.figure()
# plt.plot(feature[:, 60])
# plt.figure()
# plt.plot(feature[:, 61])
# plt.figure()
# plt.plot(feature[:, 62])
# plt.figure()
# plt.plot(feature[:, 63])
# plt.figure()
# plt.plot(feature[:, 64])
# plt.figure()
# plt.plot(feature[:, 65])
# plt.show()

# scipy.io.savemat('full_data.mat', {'full_data':full_data})
# scipy.io.savemat('label_full_data.mat', {'label_full_data':label_full_data})
scipy.io.savemat('Data/testdata_1dfeature_n005.mat', {'feature':feature})


