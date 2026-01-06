from MITgcmutils import mds
import matplotlib.pyplot as plt
XC = mds.rdmds('XC'); YC = mds.rdmds('YC')
Eta = mds.rdmds('Eta', 77760)
plt.contourf(XC, YC, Eta, np.linspace(-0.02, 0.05,8), cmap='hot_r')
plt.colorbar(); plt.show()
