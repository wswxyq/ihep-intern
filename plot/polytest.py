#from tkinter import *
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

x = np.linspace(4200,4600,50)
y = (3.17395e+01)*x+(1.21695e-02)*x*x+(-4.02805e-06)*x*x*x

plt.figure()
plt.ylim(0,43500)
plt.plot(x,y)
plt.show()
