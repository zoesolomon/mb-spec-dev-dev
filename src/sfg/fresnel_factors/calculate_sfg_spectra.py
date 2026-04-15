import numpy as np
import math
import cmath
import sys
#import matplotlib.pyplot as plt

# Adapted by Kelly Hunter from matlab code 
# provided by Muhammet Mammetkuliyev (grad student of Alex Benderskii)

# J. Chem. Phys. 144, 244711 (2016)   shen
# (J. Phys. Chem. B, Vol. 107, No. 2, 2003) richmond
# considers the effect of removing Fresnel factors from SFG spectra
# J. Phys. Chem. Lett. 2011, 2, 1946-1949 calculates for air-water

def fresnel_factor(pol_type, n1, n2, ni, theta):
    gamma_here = cmath.asin(n1*cmath.sin(theta)/n2)
    if pol_type == 'xx':
        return 2*n1*cmath.cos(gamma_here) / (n1*cmath.cos(gamma_here) + n2*cmath.cos(theta))
    elif pol_type == 'yy':
        return 2*n1*cmath.cos(theta) / (n1*cmath.cos(theta)+n2*cmath.cos(gamma_here))
    else:
        return 2*n2*cmath.cos(theta) * (n1/ni)**2 / (n1*cmath.cos(gamma_here)+n2*cmath.cos(theta))

if (len(sys.argv) != 8):
  print("Usage: ./script.py xxz_file zzz_file xzx_file vis_angle ir_angle NR_backgrd redshift")
  sys.exit()

xxz_imag = []
xxz_real = []
freq_axis = []
zzz_imag = []
zzz_real = []
xzx_imag = []
xzx_real = []
zxx_imag = []
zxx_real = []
xW_n = []
yW_n = []
xW_k = []
yW_k = []

# Nonresonant background (in a.u.)
nonres_bg = float(sys.argv[6])

# Redshift applied to spectrum (if classical)
shift = float(sys.argv[7])

xxz = open(sys.argv[1]).readlines()
for i in range(0,len(xxz)):
    xxz_imag.append(float(xxz[i].split()[1]))
    xxz_real.append(float(xxz[i].split()[2]) + nonres_bg)
zzz = open(sys.argv[2]).readlines()
for i in range(0,len(zzz)):
    freq_axis.append(float(zzz[i].split()[0]) - shift)
    zzz_imag.append(float(zzz[i].split()[1]))
    zzz_real.append(float(zzz[i].split()[2]) + nonres_bg)
xzx = open(sys.argv[3]).readlines()
for i in range(0,len(xzx)):
    xzx_imag.append(float(xzx[i].split()[1]))
    xzx_real.append(float(xzx[i].split()[2]) + nonres_bg)
zxx = xzx
for i in range(0,len(zxx)):
    zxx_imag.append(float(zxx[i].split()[1]))
    zxx_real.append(float(zxx[i].split()[2]) + nonres_bg)

# Load in refractive index files
W_n = open("refractive_index_water_n.txt").readlines()
for i in range(0,len(W_n)):
    xW_n.append(float(W_n[i].split()[0]))
    yW_n.append(float(W_n[i].split()[1]))
W_k = open("refractive_index_water_k.txt").readlines()
for i in range(0,len(W_k)):
    xW_k.append(float(W_k[i].split()[0]))
    yW_k.append(float(W_k[i].split()[1]))

ssp_amplitude = []
ssp_intensity = []
pss_amplitude = []
pss_intensity = []
sps_amplitude = []
sps_intensity = []
ppp_amplitude = []
ppp_intensity = []

# Output files
ssp_stream = open("sfg_ssp_nr_"+str(nonres_bg)+"_rs_"+str(shift)+".dat", 'w')
pss_stream = open("sfg_pss_nr_"+str(nonres_bg)+"_rs_"+str(shift)+".dat", 'w')
sps_stream = open("sfg_sps_nr_"+str(nonres_bg)+"_rs_"+str(shift)+".dat", 'w')
ppp_stream = open("sfg_ppp_nr_"+str(nonres_bg)+"_rs_"+str(shift)+".dat", 'w')

for i in range(0,len(freq_axis)):
    vis_um = 0.8                 #vis wavelength of laser in micron (0.8 = 800 nm)
    IR_um = 10**4 / float(freq_axis[i])     
    c = 3*10**8
    cm1toGHz = 100*c / 10**9
    omega1 = c / (vis_um * 1000)                   #GHz, visible
    omega2 = float(freq_axis[i]) * cm1toGHz        #GHz, infrared
    omega = omega1 + omega2
    sfg_um = c / (omega*1000)
    
    n1_omega1 = 1     #air
    n1_omega2 = 1     #air
    n1_omega = 1      #air
    
    #H2O
    n2_omega1 = (np.interp(vis_um,xW_n,yW_n) + 1j*np.interp(vis_um,xW_k,yW_k))
    #X. Zhuang, P. B. Miranda, D. Kim, and Y. R. Shen, Phys. Rev. B 59, 12632 (1999)
    n = n2_omega1
    ni_omega1 = ((n**4 + 5*(n**2))/(4*(n**2) + 2))**0.5
    
    n2_omega2 = (np.interp(IR_um,xW_n,yW_n) + 1j*np.interp(IR_um,xW_k,yW_k))
    n = n2_omega2
    ni_omega2 = ((n**4 + 5*(n**2))/(4*(n**2) + 2))**0.5
    
    n2_omega = (np.interp(sfg_um,xW_n,yW_n) + 1j*np.interp(sfg_um,xW_k,yW_k))
    n = n2_omega
    ni_omega = ((n**4 + 5*(n**2))/(4*(n**2) + 2))**0.5
    
    vis_angle_value = float(sys.argv[4])     #incoming laser angles wrt surface normal
    ir_angle_value = float(sys.argv[5])
    #vis_angle_value = 67
    #ir_angle_value = 63

    beta1 = vis_angle_value*math.pi/180  # visible
    beta2 = ir_angle_value*math.pi/180  # IR
    beta = math.asin((omega1*math.sin(beta1) + omega2*math.sin(beta2))/omega)  # sfg
    
    Lxx_omega   = fresnel_factor('xx', n1_omega,  n2_omega,  ni_omega,  beta)
    Lxx_omega1  = fresnel_factor('xx', n1_omega1, n2_omega1, ni_omega1, beta1)
    Lxx_omega2  = fresnel_factor('xx', n1_omega2, n2_omega2, ni_omega2, beta2)
    Lyy_omega   = fresnel_factor('yy', n1_omega,  n2_omega,  ni_omega,  beta)
    Lyy_omega1  = fresnel_factor('yy', n1_omega1, n2_omega1, ni_omega1, beta1)
    Lyy_omega2  = fresnel_factor('yy', n1_omega2, n2_omega2, ni_omega2, beta2)
    Lzz_omega   = fresnel_factor('zz', n1_omega,  n2_omega,  ni_omega,  beta)
    Lzz_omega1  = fresnel_factor('zz', n1_omega1, n2_omega1, ni_omega1, beta1)
    Lzz_omega2  = fresnel_factor('zz', n1_omega2, n2_omega2, ni_omega2, beta2)
    
    ssp_amplitude.append(Lzz_omega2*Lyy_omega1*Lyy_omega* cmath.sin(beta2) *(xxz_real[i]+1j*xxz_imag[i]))
    ssp_intensity.append(abs(ssp_amplitude[i])**2)
    ssp_string = str(freq_axis[i]).rjust(8) + "     " + str(ssp_intensity[i]).rjust(16) + "\n"
    ssp_stream.write(ssp_string)
    pss_amplitude.append(Lyy_omega2*Lyy_omega1*Lzz_omega* cmath.sin(beta) *(zxx_real[i]+1j*zxx_imag[i]))
    pss_intensity.append(abs(pss_amplitude[i])**2)
    pss_string = str(freq_axis[i]).rjust(8) + "     " + str(pss_intensity[i]).rjust(16) + "\n"
    pss_stream.write(pss_string)
    sps_amplitude.append(Lyy_omega2*Lzz_omega1*Lyy_omega* cmath.sin(beta1) *(xzx_real[i]+1j*xzx_imag[i]))
    sps_intensity.append(abs(sps_amplitude[i])**2)
    sps_string = str(freq_axis[i]).rjust(8) + "     " + str(sps_intensity[i]).rjust(16) + "\n"
    sps_stream.write(sps_string)
    
    ppp_amplitude.append(- Lzz_omega2*Lxx_omega1*Lxx_omega * cmath.cos(beta)*cmath.cos(beta1)*cmath.sin(beta2) * (xxz_real[i]+1j*xxz_imag[i]) - Lxx_omega2*Lzz_omega1*Lxx_omega * cmath.cos(beta)*cmath.sin(beta1)*cmath.cos(beta2) * (xzx_real[i]+1j*xzx_imag[i]) + Lxx_omega2*Lxx_omega1*Lzz_omega * cmath.sin(beta)*cmath.cos(beta1)*cmath.cos(beta2) * (zxx_real[i]+1j*zxx_imag[i]) + Lzz_omega2*Lzz_omega1*Lzz_omega * cmath.sin(beta)*cmath.sin(beta1)*cmath.sin(beta2) * (zzz_real[i]+1j*zzz_imag[i]))
    ppp_intensity.append(abs(ppp_amplitude[i])**2)
    ppp_string = str(freq_axis[i]).rjust(8) + "     " + str(ppp_intensity[i]).rjust(16) + "\n"
    ppp_stream.write(ppp_string)

ssp_stream.close()
pss_stream.close()
sps_stream.close()
ppp_stream.close()
