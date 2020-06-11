import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

with open('/Users/marcoturrini/Documents/EPFL/ENDLab/CNN_tests/LRtest/trainHistory.p','rb') as file_nt:
              opt = pickle.load(file_nt)
            
sns.set(color_codes=True)

I = range(1,len(opt)+1)

x = np.double(I)/100000
y = opt

fig, ax = plt.subplots()
ax.set_xlabel('Learning rate') 
ax.set_ylabel('val loss')
ax.plot(x,y)
# ax1 = sns.regplot(x=x, y=y,order=5,ci=50)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.show()



