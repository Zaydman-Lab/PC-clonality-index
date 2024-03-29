"""
visualize.py

Created by: Mark A. Zaydman 
11/29/2022

Purpose: create visualizations for PC clonality package

"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import transforms
from matplotlib.lines import Line2D   
from typing import Callable
import pandas as pd 

def plot_sflc(X_nonmg: np.array, pc2_RI: list[float,float], z_transform: Callable[[np.array,bool],np.array], pc_transform: Callable[[np.array,bool],np.array], X_mg: np.array = [], X_cases: np.array = [], parameters: dict = False, performance: pd.DataFrame = pd.DataFrame())->None:
  """Saves .png image of plot of non-MG cohort with manufacturer's sFLC-ratio-based and PC2-based reference intervals superimposed
  If case_2 is invoked, saves a .png image of plot of MG cohort and non-MG cohort with reference intervals superimposed
  If case_3 is invoked, saves a .png image of cases with reference intervals superimposed"""

  # plot non-mg cohort as scatter plot
  fig,ax=plt.subplots()
  fig.patch.set_facecolor('xkcd:white')
  fig.set_size_inches(6,6)
  sns.scatterplot(x=X_nonmg[:,0],y=X_nonmg[:,1], color='grey')
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlabel(r'$[\kappa]\quad$(g/dL)')
  ax.set_ylabel(r'$[\lambda]\quad$(g/dL)')
  ax.spines.right.set_visible(False)
  ax.spines.top.set_visible(False)
  xlims=ax.get_xlim()
  ylims=ax.get_ylim()

  # Add manufacturer's sFLC-ratio-based reference interval
  Katz_lb=0.26
  Katz_ub=1.65
  Katz_LB=np.column_stack((np.linspace(min(X_nonmg[:,0]),max(X_nonmg[:,0]),num=50),1/Katz_lb*np.linspace(min(X_nonmg[:,0]),max(X_nonmg[:,0]),num=50)))
  Katz_UB=np.column_stack((np.linspace(min(X_nonmg[:,0]),max(X_nonmg[:,0]),num=50),1/Katz_ub*np.linspace(min(X_nonmg[:,0]),max(X_nonmg[:,0]),num=50)))    
  ax.plot(Katz_LB[:,0],Katz_LB[:,1],color='k', linestyle=':')
  ax.plot(Katz_UB[:,0],Katz_UB[:,1],color='k', linestyle=':')

  # Add pc2-based interval
  pc2_LB=np.column_stack((np.linspace(-3,7,num=50),np.ones(50)*pc2_RI[0]))
  pc2_UB=np.column_stack((np.linspace(-3,7,num=50),np.ones(50)*pc2_RI[1]))    
  LB=transforms.log_transform(z_transform(pc_transform(pc2_LB,inverse=True),inverse=True),inverse=True)
  UB=transforms.log_transform(z_transform(pc_transform(pc2_UB,inverse=True),inverse=True),inverse=True)
  if isinstance(LB,np.ndarray): #WHY IS THIS NEEDED?
    ax.plot(LB[:,0],LB[:,1],'--k')
  if isinstance(UB,np.ndarray):
    ax.plot(UB[:,0],UB[:,1],'--k')
  line_katz = Line2D([0,1],[0,1],linestyle=':', color='black')
  line_pc2 = Line2D([0,1],[0,1],linestyle='--', color='black')
  ax.legend([line_katz, line_pc2],['Katzmann: RI = 0.26≤sFLC-ratio≤1.65', 'PCCI: RI = %.2f≤PCCI≤%.2f' %(pc2_RI[0],pc2_RI[1])])
  if parameters:
    ax.set_title(r'PCCI=$%.2f*(\frac{log(kappa)-{%.2f}}{%.2f}) + %.2f*(\frac{log(lambda)-{%.2f}}{%.2f})$' %(parameters['A'],parameters['B'], parameters['C'], parameters['D'], parameters['E'], parameters['F']))
  if len(performance)>0:
    ax.legend([line_katz, line_pc2],['Katzmann: RI = 0.26≤sFLC-ratio≤1.65\n    Sensitivity=%s, Specificity=%s' %(performance.loc["Manufacturer's sFLC-ratio interval",'Sensitivity'], performance.loc["Manufacturer's sFLC-ratio interval",'Specificity']), 'PCCI: RI = %.2f≤PCCI≤%.2f\n    Sensitivity=%s, Specificity=%s' %(pc2_RI[0], pc2_RI[1], performance.loc['PC clonality index interval','Sensitivity'], performance.loc['PC clonality index interval','Specificity'])])
  plt.tight_layout()
  plt.savefig('./Output/case_1.png')

  if len(X_mg)>0:
    sns.scatterplot(x=X_mg[:,0],y=X_mg[:,1])
    plt.tight_layout()
    plt.savefig('./Output/case_2.png')

  if len(X_cases)>0:
    sns.scatterplot(x=X_cases[:,0],y=X_cases[:,1])
    plt.tight_layout()
    plt.savefig('./Output/case_3.png')
  

