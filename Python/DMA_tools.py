
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.signal import savgol_filter
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline

import math

import Aerosol_tools

class DMA():
    """
        This class is used to gather all the parameters
        of the studied DMA.
    """
    def __init__(self, DMA_props):
        self.Reset_DMA(DMA_props)
        return
    def Reset_DMA(self, DMA_props):
        self.model = DMA_props["model"]
        if (self.model == "long"):
            self.Dp_range = np.array([10, 487])       # nm
            self.n_range = np.array([1e+2, 1e+07])    # part./cm^3
        self.v_range = np.array([0.998e+0, 1.0001e+04])
        self.flow_aerosol = DMA_props["flow_aerosol"]   # L/min
        self.flow_sheath = DMA_props["flow_sheath"]     # L/min
        # self.flow_sampling = DMA_props["flow_sampling"] # L/min
        self.Load_parameters_model()
        self.L_eff = DMA_props["L_eff"]   # m
        # nano DMA: self.L = 4.987e-02    # m
        self.A = self.L/np.log(self.R2/self.R1)
        self.V = DMA_props["Voltage"]     # V
        self.Zp = (2*self.flow_sheath)*(1.66667e-5)/(4*np.pi*self.A*self.V)
        # Ambient properties
        self.Pressure = DMA_props["Pressure"]         # Pa
        self.Temperature = DMA_props["Temperature"]   # K
        self.Charge_limit = DMA_props["Charge_limit"]
        return
    def Load_parameters_model(self):
        if(self.model == "nano"):
            self.R1 = 0.937e-02
            self.R2 = 1.905e-02
            self.L = 4.987e-02
        elif(self.model == "long"):
            self.R1 = 9.37e-03
            self.R2 = 19.61e-03
            self.L = 44.37e-02
        return
    def Calculate_Zp(self, voltage):
        """
         Determine the electrical mobility (Zp) from
         the properties of the DMA and voltage.
         
         Be careful, flow_sheath is assumed in units: lpm
        """
        Zp = (2*self.flow_sheath)*(1.66667e-5)/(4*np.pi*self.A*voltage)
        return Zp
    def Update_voltage(self,new_V):
        if ((new_V >= self.v_range[0]) and (new_V <= self.v_range[1])):
            self.V = new_V
            self.Zp = self.Calculate_Zp(new_V)
        else:
            print("Error: voltage out of range for this DMA, V=",new_V)
            self.V = 0
            self.Zp = 0
        return
    def Electric_mobility(self,Dp,charge):
        f = Aerosol_tools.friction(Dp,self.Temperature)
        Z = charge * Aerosol_tools.q_e/f
        return Z
    def Diffusion_coefficient(self,Dp):
        """
        Determines the diffusion coefficient of a particle
        given its diameter Dp, temperature and pressure.
        """
        f = Aerosol_tools.friction(Dp,self.Temperature)
        Diff = Aerosol_tools.k_B * self.Temperature/f
        return Diff
    def Penetration_efficiency_Tl(self,Dp):
        """
        Penetration efficiency through the DMA is computed based on the
        parameterized results by Reineking & Porstendörfer (1986).
        $$T_l = 0.82\exp(-11.5u)+0.1\exp(-70.0u)+0.03\exp(-180.0u)+0.02\exp(-340.0u)$$
        where $u = \frac{D_{ab} l_{eff}}{q_{sa}}$, $l_{eff}$ is the parameterized
        effective diffusion length, and $q_{sa}$ is the aerosol flow rate through the DMA,
        and $D_{ab}$ is the particle diffusion coefficient.
        """
        Diff = self.Diffusion_coefficient(Dp)
        q_sa = self.flow_aerosol * (1.66667e-5)
        u = Diff * self.L_eff/q_sa
        t = np.zeros(4)
        t[0] = 0.82 * np.exp(-11.5 * u)
        t[1] = 0.10 * np.exp(-70.0 * u)
        t[2] = -0.03 * np.exp(-180.0 * u)
        t[3] = -0.02 * np.exp(-340.0 * u)
        return np.sum(t)
    def Charging_efficiency_Tc(self,Dp,charge):
        """
        Charging efficiency (charge equilibrium) obtained in the bipolar charger
        is computed based on the parameterized measurements by Wiedensohler et al. (1988)
        with coefficients taken from the TSI 3080 Manual (2009). 
        $$T_c(k) = 10^\left\{ \sum_{i=1}^6 a_i (k) \left[ \ln \left(\frac{D_p}{nm}\right) \right]^{i-1} \right \}$$
        where $k = -2,-1,1,2$ is the number and polarity of particle charge and $a_i$ are
        empirical coefficients based on Fuchs theory. This equation is valid for 
        $1\, nm<D_p<1000\, nm$ for $k=[-1, 1]$ and $20\, nm<D_p<1000\, nm$ for $k=[-2, 2]$.
        
        ![alt text](Figures/params_Wiedensohler.png)
        
        For $k \ge \pm 3$, the formula from the TSI manual is used:
        $$T_c(k) = \frac{e}{\sqrt{4\pi^2\epsilon D_pk_bTK_e}} \exp \left( \frac{-\frac{\left[|k| - 2\pi\epsilon D_pk_bT K_e\ln(0.875)\right]^2}{e^2}}{ \frac{4\pi\epsilon D_pk_bT}{e^2}} \right) $$
        where $e$ is the elementary charge and $\epsilon$ is the dielectric constant for air.
        """
        # Charges [-2, +2]: Wiedensohler et al. (1988) parameters
        if((charge >= -2) and (charge <= 2)):
            a = {
            "-2": [-26.3328,35.9044,-21.4608,
                  7.0867,-1.3088,0.1051],
            "-1": [-2.3197,0.6175,0.6201,
                  -0.1105,-0.1260,0.0297],
            "0": [-0.0003,-0.1014,0.3073,
                  -0.3372,0.1023,-0.0105],
            "1": [-2.3484,0.6044,0.4800,
                 0.0013,-0.1544,0.0320],
            "2": [-44.4756,79.3772,-62.8900,
                 26.4492,-5.7480,0.5059] }
            ch_exp = 0
            ai = a[str(charge)]
            for i in range(len(ai)):
                ch_exp += ai[i] * np.power(np.log10(Dp*1e+09),i)
            ch_frac = np.power(10,ch_exp)
        else:
            K_e = 9e+09 # N*m^2/C^2
            temp10 = np.pi * Dp *Aerosol_tools.k_B * self.Temperature/np.power(Aerosol_tools.q_e,2)/K_e
            temp11 = np.power(charge,2)
            temp12 = temp10 / np.pi
            ch_frac = (1/np.sqrt(temp10)) * np.exp(-temp11/temp12)
        return ch_frac
    def erf_func(self,x):
        er = x*math.erf(x)+np.exp(-x**2)/np.sqrt(np.pi)
        return er
    def Transfer_function_Tf(self,Dp,voltage,charge):
        """
        The DMA transfer function is the probability that a particle of a given size
        exits the classifier via the sample flow. The diffusive broadened DMA transfer
        function is computed assuming blanced sheath and excess flows using the expression 
        of Stolzenburg and McMurry (2008).
        $$ \Omega(\tilde{z},\beta,\sigma) = \frac{\sigma}{\sqrt{2}\beta}\left[\epsilon \left( \frac{\tilde{z}-(1+\beta)}{\sqrt{2}\sigma} \right) + \epsilon \left (\frac{\tilde{z}-(1-\beta)}{\sqrt{2}\sigma} \right) - 2\epsilon \left ( \frac{\tilde{z}-1}{\sqrt{2}\sigma}\right)  \right]$$
    
        where $\tilde{z} = \frac{z}{z^s}$ is the dimensionless mobility, $z$ is the particle mobility $z^s$ is the centroid mobility selected by the DMA, $\epsilon = x \mathrm{erf}(x) +\left(\exp(-x^2)/\sqrt{\pi}\right)$, $\mathrm{erf}$ is the error function, and $\beta = \frac{q_{sa}}{q_{sh}}$. The parameter $\sigma$ accounts for diffusional broading of the transfer function. Assuming plug flow, $\sigma$ can be computed using the following equations Hagwood (1999) <br> <br>
        $\gamma = \left(\frac{r_1}{r_2}\right)^2$ <br>
        $I = \frac{1}{2}(1+γ)$ <br>
        $\kappa = \frac{lr_2}{r_2^2-r_1^2}$ <br>
        $G = \frac{4(1+\beta)^2}{(1-γ)} \left[I+\{2(1+\beta)\kappa\}^{-2} \right ] $<br>
        $\sigma = \sqrt{\frac{2G\pi ld_{ab}}{q_{sh}}}$
        """
        Diff = self.Diffusion_coefficient(Dp)
        beta = self.flow_aerosol/self.flow_sheath
        gamma = np.power(self.R1/self.R2,2)
        I = 0.5 * (1 + gamma)
        kappa = self.L * self.R2/(np.power(self.R2,2)-np.power(self.R1,2))
        G = 4 *np.power(1+beta,2)/(1-gamma) * (I+np.power(2*(1+beta)*kappa,-2))
        sigma = np.sqrt(2*G*np.pi*self.L*Diff/(self.flow_sheath*1.66667e-5))
        Z = self.Electric_mobility(Dp,charge)
        Zp = self.Calculate_Zp(voltage)#/charge
        z_til = np.abs(Z/Zp)
        # TODO: Check why this is giving negative values in extreme cases
        tf = sigma/(np.sqrt(2)*beta)*(self.erf_func((z_til-(1+beta))/(np.sqrt(2)*sigma))+\
                                      self.erf_func((z_til-(1-beta))/(np.sqrt(2)*sigma))-\
                                      2*self.erf_func((z_til-1)/(np.sqrt(2)*sigma)))
        return np.max([tf,1e-60])
    def Convolution_matrix(self, Dp, Voltage):
        """
        This function determined the convolution matrix consisting
        ot the product of the penetration efficiency (Tl), the charging
        efficiency (Tc), and the DMA transfer function (Tf).
        """
        n = len(Dp)
        k_vec = np.array(range(self.Charge_limit)) + 1
        A = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                Tl = self.Penetration_efficiency_Tl(Dp[j])
                for charge in k_vec:
                    Tc = self.Charging_efficiency_Tc(Dp[j],charge)
                    Tf = self.Transfer_function_Tf(Dp[j],Voltage[i],charge)
                    A[i,j] += Tf * Tc * Tl
        self.A_mat = A
        return A
    def Convolution_matrix_simplified(self, Dp, Voltage):
        """
        This function determined the convolution matrix consisting
        ot the product of the penetration efficiency (Tl), the charging
        efficiency (Tc), and the DMA transfer function (Tf).
        """
        n = len(Dp)
        k_vec = np.array(range(self.Charge_limit)) + 1
        A_mat = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                Tl = 1 #self.Penetration_efficiency_Tl(Dp[j])
                for charge in k_vec:
                    Tc = self.Charging_efficiency_Tc(Dp[j],charge)
                    Tf = self.Transfer_function_Tf(Dp[j],Voltage[i],charge)
                    A_mat[i,j] += Tf * Tc * Tl
        self.A_mat = A_mat
        return A_mat
    def Total_PSD_probability(self, Dp_j, Dp_vec, f_vec):
        """
        This function is used to the determine the probability to
        have a given particle size.
        """
        y = np.interp(Dp_j, Dp_vec, f_vec)
        return y
    def Convolution_matrix_B(self, Dp, Voltage, Dp_vec, f_vec):
        """
        This function determined the B convolution matrix consisting
        ot the product of the penetration efficiency (Tl), the size
        probability (Ts), and the DMA transfer function (Tf).
        """
        m = len(Dp)
        k_vec = np.array(range(m)) + 1
        B = np.zeros((m,m))
        for i in range(m):
            for j in range(m):
                charge = k_vec[j]
                for k in range(m):
                    Ts = self.Total_PSD_probability(Dp[k], Dp_vec, f_vec)
                    Tl = self.Penetration_efficiency_Tl(Dp[k])
                    Tf = self.Transfer_function_Tf(Dp[k],Voltage[i],charge)
                    B[i,j] += Tf * Ts * Tl
        return B,k_vec
    def Convolution_matrix_C(self, Dp, Charges, voltage):
        """
        This function determined the convolution matrix consisting
        of the product of the penetration efficiency (Tl), and the
        DMA transfer function (Tf).
        """
        n = len(Dp)
        m = len(Charges)
        C = np.zeros((n,m))
        for i in range(n):
            Tl = self.Penetration_efficiency_Tl(Dp[i])
            for j in range(m):
                Tf = self.Transfer_function_Tf(Dp[i],voltage,Charges[j])
                C[i,j] = Tf * Tl
        return C
    def Convolution(self, A, N):
        """
        In order to check the inversion process, this function gives
        the convolution of the original signal (N) according to the
        determined convolution matrix (A):
        $$ R = A\cdot N$$
        Where R is the convoluted signal.
        """
        R_conv = np.zeros_like(N)
        for i in range(len(N)):
            R_conv[i] = np.sum(A[i,:] * N)
        return R_conv
    # http://g2s3.com/labs/notebooks/inverseProblemPrototype.html
    def Deconvolution(self, R, A, alpha):
        """
        The inversion is ill-posed and therefore a regularization
        approach is necessary. Here we apply Tikhonov with parameter
        alpha (assumed to be known).
        """
        H = np.dot( A.transpose(), A) + alpha*np.identity(A.shape[1])
        rhs = np.dot( A.transpose(), R)
        return np.linalg.solve(H, rhs)
    def Boltzmann_charge(self, charge, Dp):
        """
        Charge distribution selected according to Boltzmann equation.
        """
        K_e = 9e+09 # N*m^2/C^2
        temp10 = np.pi * Dp *Aerosol_tools.k_B * self.Temperature/np.power(Aerosol_tools.q_e,2)/K_e
        temp11 = np.power(charge,2)
        temp12 = temp10 / np.pi
        ch_frac = (1/np.sqrt(temp10)) * np.exp(-temp11/temp12)
        return ch_frac
    def Error_alpha_parameter(self, x):
        """
        Determines the error of the current alpha value (x)
        according to the sum of the squared differences of
        the observed and the convoluted values.
        """
        N_alpha = self.Deconvolution(self.R, self.A_mat, x)
        R_conv = self.Convolution(self.A_mat, N_alpha)
        error = np.power(R_conv/self.R - 1, 4)
        return np.sum(error)
    def Find_optimum_alpha(self, R, A_mat, x0=1e-01):
        """
        Conducts an optimization process where the alpha parameter
        that minimizes the sum of the squared errors between the
        convoluted and the observed values is searched.
        """
        self.R = R
        self.A_mat = A_mat
        # methods: BFGS, nelder-mead, Newton-CG (requires Jacobian), trust-ncg (requires Jacobian)
        res = minimize(self.Error_alpha_parameter, x0=x0, method='BFGS',
               options={'xatol': 1e-11, 'disp': True})
        return res.x
        
def Noise_reduction(y,window_size,poly_order,mod):
    yhat = savgol_filter(y,window_size,poly_order,mode=mod)
    return yhat
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def Deconvolute_singleDMA(dma,measurements_data,smooth_n=10,interp_points=100,alpha = 1e-03,check=False):
    Dp0 = measurements_data.Dp.values * 1e-09
    Voltage0 = measurements_data.V2.values
    R0 = measurements_data.R.values
    
    cs_v = CubicSpline(Dp0, Voltage0)
    cs_R = CubicSpline(Dp0, R0)
    
    #Dp = np.linspace(np.min(Dp0), np.max(Dp0), interp_points)
    Dp = np.logspace(np.log10(np.min(Dp0)), np.log10(np.max(Dp0)), interp_points)
    Voltage = cs_v(Dp)
    R = cs_R(Dp)
    
    A = dma.Convolution_matrix(Dp, Voltage)
    #A = dma0.Convolution_matrix_simplified(Dp, Voltage)
    R_nr = Noise_reduction(R,5,2,"interp")
    
    N_alpha = dma.Deconvolution(R_nr, A, alpha)
    R_conv = dma.Convolution(A, N_alpha)
    
    dNdLogD = np.zeros_like(N_alpha)
    dNdLogD[0] = N_alpha[0] / np.log10(Dp[1]/Dp[0])
    for i in range(1,len(N_alpha)):
        dNdLogD[i] = N_alpha[i] / np.log10(Dp[i]/Dp[i-1])
    
    dNdLogD_r = smooth(dNdLogD,smooth_n)
    
    if (check==True):
        fig, ax = plt.subplots(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
        plt.plot(Dp0 * 1e+09, R0,"o b", label="Read DMA (original)")
        plt.plot(Dp * 1e+09, R,". r", label="Read DMA (cubic spline)")
        plt.plot(Dp * 1e+09, R_conv,"-g", label="R_conv")
        plt.ylabel("Number concentration (cm⁻³)", fontsize=20)
        plt.xlabel("Diameter (nm)", fontsize=20)
        ax.tick_params(direction='in', length=6, width=1, colors='k',
               grid_color='k', grid_alpha=0.5)
        ax.tick_params(axis='x', which='minor', direction='in')#,bottom=False)
        ax.tick_params(axis='y', which='minor', direction='in')#,bottom=False)
        plt.rc('xtick', labelsize=16); plt.rc('ytick', labelsize=16)
        plt.legend(fontsize=20); plt.grid()
        plt.show()
        
        fig, ax = plt.subplots(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
        plt.plot(Dp * 1e+09, dNdLogD,"-g", label="N_alpha")
        plt.ylabel("dN/dLog(Dp) (cm$^{-3}$)", fontsize=20)
        plt.xlabel("Diameter (nm)", fontsize=20)
        ax.tick_params(direction='in', length=6, width=1, colors='k',
               grid_color='k', grid_alpha=0.5)
        ax.tick_params(axis='x', which='minor', direction='in')#,bottom=False)
        ax.tick_params(axis='y', which='minor', direction='in')#,bottom=False)
        plt.rc('xtick', labelsize=16); plt.rc('ytick', labelsize=16)
        plt.legend(fontsize=20); plt.grid()
        plt.show()
    
    dNdLogD = np.interp(Dp0,Dp,dNdLogD)
    R_conv = np.interp(Dp0,Dp,R_conv)
    dNdLogD_r = np.interp(Dp0,Dp,dNdLogD_r)
    return dNdLogD,dNdLogD_r,R_conv
