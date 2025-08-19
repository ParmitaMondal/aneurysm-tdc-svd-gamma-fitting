import os
import numpy as np
import pandas as pd
import pydicom
import matplotlib.pyplot as plt
from scipy.signal import medfilt
from scipy.signal import correlate
from PIL import Image
from scipy import interpolate
from scipy.linalg import toeplitz, svd, pinv
from openpyxl import load_workbook
from scipy.stats import gamma
from scipy import signal
from scipy.optimize import curve_fit 
import math


def API_parameters(x, y):
    j = 0
    BAT_index = 0
    while j < len(x):
        if (y[j] > np.max(y)*0.1):
            BAT_index = j-1
            break
        j += 1
  
    PH=np.max(y) 
    j = 0
    PH_index = 0
    while j < len(x):
        if (y[j] == PH):
            PH_index = j
            break
        j += 1  
    TTP=x[PH_index]-x[BAT_index]
    AUC=np.trapz(y, x)
    MTT=np.trapz(y*x, x)/AUC
    derivative=np.gradient(y[BAT_index:PH_index])
    max_Df=derivative.max()
    return x[BAT_index], MTT, TTP, PH, AUC, max_Df

def API_SVD(x_inl, y_inl, x_ANY, y_ANY):
    j = 0
    BAT_index = 0
    while j < len(x_inl):
        if (y_inl[j] > np.max(y_inl)*0.1):
            BAT_index = j-1
            break
        j += 1
    P_SVD=0.1
   
    inlet_zeros = np.zeros([len(y_inl)])
    inlet_toeplitz = toeplitz(y_inl, inlet_zeros) 
    V, S, U = svd(np.matrix.transpose(inlet_toeplitz))
    S = np.diag(S)
    truncate=P_SVD*np.max(S)
    lam=0.1
    S[S<truncate]=0
    S[S>truncate]= S[S>truncate]/(S[S>truncate]**2 + lam**2) #Tikhonov regression
    IRF = abs(np.matmul(np.matmul(V, S), np.matmul(U, y_ANY)))
    IRF1 = IRF/np.max(IRF) 
    RBF = np.max(IRF)
    RBV = np.trapz(IRF, x_inl)
    MTT_SVD = RBV / RBF 
    return RBF, RBV, MTT_SVD, IRF

mask_aneurysm_path = 'path...' 
mask_inlet_path='path...'
excel_file = 'path...'
sheet = 'Sheet1'
save_file = 'path...'
df=load_workbook(save_file)
ws=df.worksheets[0]
caseInfo = pd.read_excel(excel_file)
cases_mask = caseInfo.iloc[:,0]
cases_dicom=caseInfo.iloc[:,1]

for x in range(0,1): 
    a= 5
    y= 300
    z= 600 
    import re
    import glob  
    numbers = re.compile(r'(\d+)')
    def numericalSort(value):
      parts = numbers.split(value)
      parts[1::2] = map(int, parts[1::2])
      return parts

   filelist = sorted(glob.glob('path...'+"*.raw"), key=numericalSort)
    proj = np.empty((512,366,1000), dtype=np.float32) 
    for file, idx in enumerate(filelist):
      proj[:,:,file] = np.fromfile(open(filelist[file], 'rb'), dtype=np.float32).reshape((512,366))
    plt.imshow(proj[:,:,500], cmap='gray')
    plt.show()
     
    ANY_case=cases_dicom[x]
    ANY_case_aneurysm_mask=cases_mask[x]
    print(ANY_case)
    dicom_temp=proj
    dicom_image=dicom_temp
    dicom_image[dicom_image==-np.inf]=0
    dicom_mean=np.mean(dicom_image[:,:,400:np.shape(dicom_image)[2]-100], axis=2)
    dimension=dicom_temp.shape[0]


    TDC_average=np.zeros(shape=np.shape(dicom_image)[2], dtype=np.float32)
    TDC_inlet_average=np.zeros(shape=np.shape(dicom_image)[2], dtype=np.float32)
    time_vector=np.arange(0,1.0,0.001) 
    
    mask_aneurysm=np.zeros(shape=(dimension,dimension), dtype=np.float32)
    mask_inlet=np.zeros(shape=(dimension,dimension), dtype=np.float32)
    mask_aneurysm = Image.open(os.path.join(mask_aneurysm_path,  ANY_case_aneurysm_mask +'_0.tif'))  
    mask_inlet = Image.open(os.path.join(mask_inlet_path,  ANY_case_aneurysm_mask +'_0_inl.tif'))  
    im=np.asarray(mask_aneurysm, dtype=np.float32)
    im2=np.asarray(mask_inlet, dtype=np.float32)
    mask_display=np.asarray(mask_aneurysm)+np.asarray(mask_inlet)
      
    im=abs((im-255)/255)
    im2=abs((im2-255)/255)
    ind = np.transpose(np.nonzero(im))
    ind_inlet = np.transpose(np.nonzero(im2))
    print(im.sum(), im2.sum())

    for mm in range(0, len(ind)):
        TDC_average = TDC_average+dicom_image[ind[mm][0],ind[mm][1] ,: ]/im.sum()
    for nn in range(0, len(ind_inlet)):
        TDC_inlet_average = TDC_inlet_average+dicom_image[ind_inlet[nn][0], ind_inlet[nn][1], :]/im2.sum()
 
    TDC_average=medfilt(TDC_average, kernel_size=1)
    TDC_average[TDC_average<0]=0 
    
    plt.plot(time_vector, TDC_average, time_vector, TDC_inlet_average)
    plt.xlabel('Time')
    plt.ylabel('Contrast density')
    plt.gca().legend(('Original  TDC of the aneurysm dome ','Original  TDC of the inlet '))
    plt.show()
    
    new_time_vector=np.arange(0,3.0,0.001) 
    residual_func=API_SVD(time_vector, TDC_inlet_average, time_vector, TDC_average)
    residual_function= residual_func[3]
    residual_function1 =residual_function/np.max(residual_function) 

    from scipy.stats import gamma
    y_values = 4800*new_time_vector**2*np.exp(-new_time_vector/0.1)
    plt.plot(new_time_vector, y_values, label='Gamma Variate Function', color='blue')
    plt.title('perfect Inlet ')
    plt.show() 

    res = np.concatenate([residual_function, np.full(2000, residual_function[-1])])
    k=signal.fftconvolve(res, y_values) 
    k1=k[0:len(new_time_vector),] 
    plt.plot(new_time_vector, k1)
    plt.title('IRF convolved with ideal inlet') 
    plt.show()
    TDC_inlet= y_values
    TDC_aneu= k1
    
    
    ####################### adding patient motion #####################################                     
    y=300
    z=600
    time= new_time_vector   
    TDC_in= TDC_inlet
    TDC_an=TDC_aneu
    for i in range(len(TDC_inlet)):
      if i>y and i<z: #start from 300to 600
        TDC_in[i] = TDC_inlet[i] -2
        TDC_an[i] = TDC_aneu[i]-1
        for i in range(1000):
          if i>y and i<z: #start from 300to 600
            TDC_average[i] = TDC_average[i]-2
    plt.plot(time, TDC_an, time, TDC_in)
    plt.show()
             
    ws2=df.worksheets[1]
    column1=["C", "D","E","F","G","H"]
    column2=["K", "L","M","N","O","P"]
    popt, pcov = curve_fit(func, time, TDC_an)
    new_TDC_average=func(time,*popt) 
    plt.plot(time, new_TDC_average)
    ANY_results3=np.array(API_parameters(time, new_TDC_average))
    popt, pcov = curve_fit(func, time, TDC_in)
    new_TDC_inlet_average=func(time,*popt) #gamma variate inlet
    plt.plot(time,  new_TDC_inlet_average)
    I_inlet_results3=np.array(API_parameters(time, new_TDC_inlet_average))
    I_inlet_results3[0]=1
    cor_pre=correlate(new_TDC_inlet_average-np.mean(new_TDC_inlet_average), new_TDC_average-np.mean(new_TDC_average))/(len(new_TDC_average) * np.std(new_TDC_inlet_average) * np.std(new_TDC_average))
    ANY_results3[0]=cor_pre.max()
    plt.plot(time, new_TDC_average, time,  new_TDC_inlet_average)
    plt.xlabel('Time')
    plt.ylabel('Contrast density')
    plt.gca().legend(('gamma:TDC of the aneurysm','gamma:TDC of the inlet'), loc='upper left')
    plt.show()
  
    API_all_parameters = ANY_results3
    API_all_parameters = np.reshape(API_all_parameters, (1, -1)) 
    
    k=0
    for col in column1:
      ws2[col + str(a)] = API_all_parameters[0,k]  # Column D #starts with 5
       k= k+1
    df.save(save_file)
  
    API_parameters = I_inlet_results3
    API_parameters = np.reshape(API_parameters, (1, -1)) 
    k=0
    for col in column2:  
      ws2[col + str(a)] = API_parameters[0,k]  # Column C #starts with 5
      k= k+1 
    df.save(save_file)
    
    
    
    
    
    
    
    
    
    
