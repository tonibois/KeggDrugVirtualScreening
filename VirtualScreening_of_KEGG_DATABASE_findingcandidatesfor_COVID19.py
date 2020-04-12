
# coding: utf-8

# # Study of ATC drug properties and structure from KEGG_DRUG database
# 
# ------------
# + Antonio Oliver Gelabert, April 2020*
# + ORCID    : http://orcid.org/0000-0001-8571-2733 ) *
# + Linkedin : https://www.linkedin.com/in/aoliverg/?locale=en_US
# ------------
# 
# ATC is european systemathic classification of accepted medicines by therapeutic usage. To see details of this classification better take a look at : 
# 
# + https://en.wikipedia.org/wiki/Anatomical_Therapeutic_Chemical_Classification_System
# + https://www.whocc.no/atc_ddd_index/
# 
# KEGG_DRUG database is a repository of compounds which structure is available from: 
# 
# + https://www.genome.jp/kegg/drug/
# 
# The individual information about compounds can be consulted at LigandBox:
# 
# + http://www.mypresto5.com/ligandbox/
# 
# Compound structures contains information that are not always to structure. So there are some ways to structure information as
# + Fingerprints (MACCS,FP2,...). http://openbabel.org/wiki/Tutorial:Fingerprints
# + Molecular descriptors (USR,ElectroShape,PED,...) 
# 
# For instance, PED descriptors are described in scientific literature
# 
# + https://www.nature.com/articles/srep43738
# 
# From the curated dataset combining information of this two links:
# 
# + http://www.mypresto5.com/ligandbox/cgi-bin/lbox_download.cgi?LANG=en  (MOL2 FILES OF KEGG_DRUG DATABASE)
# + https://www.genome.jp/kegg-bin/get_htext#A1 (KEGG ATC codes)
# 
# The information can be combined in the file 'Kegg_ATC_database_info.csv'. Lets open the file and make rapid visualization of the kind of information available:

# In[1]:


import pandas as pd
kegginfo=pd.read_csv('Kegg_ATC_database_info.csv')
kegginfo.head()


# Now open the structure file that have been previously converted from MOL2 file to SMILES using OpenBabel (http://openbabel.org/wiki/Main_Page) and it has been converted with annotated properties. Let's see the command:
# 
# obabel -imol2 KEGG_DRUG.mol2 -osmi -O KEGG_DRUG_PROPS.smi --append 'MW' 'logP' 'MR' 'TPSA' 'HBA1' 'HBD' 'HBA2' 'InChIKey' 'formula' 'cansmiNS' 'InChI'
# 
# Adding a suitable header and converting to CSV: 

# In[2]:


keggdata=pd.read_csv('ATC_KEGG_DRUG_DATABASE.csv')
keggdata.head()


# Now lets combine both files to one using merge and KEGG codes:

# In[3]:


keggcompl=pd.merge(kegginfo,keggdata,on='KEGG_code')
keggcompl


# And remove duplicates:

# In[4]:


keggclean=keggcompl.drop_duplicates()
keggclean


# In[5]:


# There are a total of :
len(keggclean)    #np.shape(keggclean)[0]
# Compounds in KEGG DRUG ATC DATABASE that area analyzed here


# Now, lets take a subset of columns:

# In[6]:


keggclean[['SMILES','KEGG_code','ATC_full_code']].to_csv('Kegg_atc.smi', sep=',', encoding='utf-8')


# Or a subset of rows by focusing on specific therapeutic activities, like the A01A category, which corresponds in ATC codification to the first one (i.e. stomathological preparations https://www.atccode.com/A01A)

# In[7]:


keggsubset=keggclean[keggclean['ATC_full_code'].str[:4]=='A01A']


# Next: 
# + Export the subset to CSV format
# + Reformat from CSV to SMILES
# + Take the first 20 rows 
# + Convert rows to SVG files for visualization using openbabel

# In[8]:


keggsubset[['SMILES','KEGG_code']].to_csv('subset.smi', sep=' ', encoding='utf-8')
get_ipython().system('gawk "{print $2,$3}" subset.smi | sed \'1d\'  > subsetcur.smi')
get_ipython().system('head -n 20 subsetcur.smi > subsetcur20.smi')
get_ipython().system('obabel -ismi subsetcur20.smi -osvg -O subset.svg')


# Display in a notebook environment the created SVG:

# In[9]:


from IPython.display import SVG, display
display(SVG('subset.svg'))


# Next, we can search for specific compound name, like: Ibuprofen or Danuravir

# In[10]:


keggmol=keggclean[keggclean['CompoundName'].str[:10]==' Ibuprofen']
keggmol


# In[11]:


keggmol=keggclean[keggclean['CompoundName'].str[:10]==' Darunavir']
keggmol


# Lets apply now this analysis to search for accepted drugs as antivirals of the **protease inhibitors (ATC code = J05AE)**, one of the approaches used to search for **COVID-19 potential treatment**:
# 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4774534/

# In[12]:


keggsubset2=keggclean[keggclean['ATC_full_code'].str[:5]=='J05AE']
keggsubset2[['SMILES','KEGG_code']].to_csv('subset2.smi', sep=' ', encoding='utf-8')
get_ipython().system('gawk "{print $2,$3}" subset2.smi | sed \'1d\'  > subsetcur2.smi')
get_ipython().system(' head -n 20 subsetcur2.smi > ProteaseinhibJ05AE.smi')
get_ipython().system('obabel -ismi ProteaseinhibJ05AE.smi -osvg -O ProteaseinhibJ05AE.svg')


# In[13]:


keggsubset2[['SMILES','KEGG_code']].head(20)


# In[14]:


display(SVG('ProteaseinhibJ05AE.svg'))


# Now lets study some statistical properties of that subgroup of compounds like:
#     
# + Molecular Weight (MW)
# + Total Polar Surface area (TPSA)
# + Octanol-Water Partition coefficient (logP)
# + Hydrogen Bond Donors (HBD)
# + Hydrogen Bond Acceptors (HBA)
# + Molar Refractivity (MR)
# 

# In[15]:


import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 3, constrained_layout=True,figsize=(15,15))

axs[0, 0].hist(keggsubset2['MW']);
axs[0, 0].set_title('MW_J05AE')
axs[0, 1].hist(keggsubset2['logP']);
axs[0, 1].set_title('logP_J05AE')
axs[0, 2].hist(keggsubset2['HBA1']);
axs[0, 2].set_title('HBA1_J05AE')
axs[1, 1].hist(keggsubset2['HBD']);
axs[1, 1].set_title('HBD_J05AE')
axs[1, 0].hist(keggsubset2['TPSA']);
axs[1, 0].set_title('TPSA_J05AE')
axs[1, 2].hist(keggsubset2['MR']);
axs[1, 2].set_title('MR_J05AE')

for ax in axs.flat:
    ax.set(xlabel='', ylabel='')

plt.tight_layout()
#plt.subplots_adjust(left= 0.1, bottom=0.05 , right=0.9, top=0.95, wspace=0.5, hspace=0.5)

plt.show()


# In[16]:


import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 3, constrained_layout=True,figsize=(15,15))

axs[0, 0].boxplot(keggsubset2['MW'],labels=['J05AE']);
axs[0, 0].set_title('MW')
axs[0, 1].boxplot(keggsubset2['logP'],labels=['J05AE']);
axs[0, 1].set_title('logP')
axs[0, 2].boxplot(keggsubset2['HBA1'],labels=['J05AE']);
axs[0, 2].set_title('HBA1')
axs[1, 1].boxplot(keggsubset2['HBD'],labels=['J05AE']);
axs[1, 1].set_title('HBD')
axs[1, 0].boxplot(keggsubset2['TPSA'],labels=['J05AE']);
axs[1, 0].set_title('TPSA')
axs[1, 2].boxplot(keggsubset2['MR'],labels=['J05AE']);
axs[1, 2].set_title('MR')

for ax in axs.flat:
    ax.set(xlabel='', ylabel='')

plt.tight_layout()
#plt.subplots_adjust(left= 0.1, bottom=0.05 , right=0.9, top=0.95, wspace=0.5, hspace=0.5)

plt.show()


# We can explore as well some statistical parameters for Molecular Weight of this subgroup:

# In[17]:


import numpy as np
from scipy import io
from scipy.stats import kurtosis,skew,itemfreq,describe,moment,mode,sem,hmean

data=[]
param='MW'
data=keggsubset2[param]
print("***********************************")
print("Statistical Summary for "+param+" :")
print("***********************************")
print("mean:",round(np.mean(data),2))
print("mode counts:",mode(data)[1][0])
print("mode :",round(mode(data)[0][0],2))
print("SEM (standard error of the mean):",round(sem(data),2))
print("std:",round(np.std(data),2))
print('skew:',round(skew(data),2))
print("Kurtosis:",round(kurtosis(data),2))
print("max:",round(np.max(data),2))
print("min:",round(np.min(data),2))
print("Percentile 25 :",round(np.percentile(data,25),2))
print("Percentile 75 :",round(np.percentile(data,75),2))
print("Percentile 10 :",round(np.percentile(data,10),2))
print("Percentile 90 :",round(np.percentile(data,90),2))
print("***********************************")
print("***********************************")


# In[18]:


import numpy as np
from scipy import io
from scipy.stats import kurtosis,skew,itemfreq,describe,moment,mode,sem,hmean

data=[]
param='logP'
data=keggsubset2[param]
print("***********************************")
print("Statistical Summary for "+param+" :")
print("***********************************")
print("mean:          ",round(np.mean(data),2))
print("mode counts:   ",mode(data)[1][0])
print("mode:          ",round(mode(data)[0][0],2))
print("SEM:           ",round(sem(data),2))
print("std:           ",round(np.std(data),2))
print('skew:          ',round(skew(data),2))
print("Kurtosis:      ",round(kurtosis(data),2))
print("max:           ",round(np.max(data),2))
print("min:           ",round(np.min(data),2))
print("Percentile 25 :",round(np.percentile(data,25),2))
print("Percentile 75 :",round(np.percentile(data,75),2))
print("Percentile 10 :",round(np.percentile(data,10),2))
print("Percentile 90 :",round(np.percentile(data,90),2))
print("***********************************")
print("***********************************")


# Following, a clustering based on descriptors as MW, logP, HBD, HBA, TPSA is performed.
# 
# The approach performed at :
# + https://builtin.com/data-science/unsupervised-learning-python (K-MEANS,TSNE,...)
# 
# SOM implementation here : 
# + https://peterwittek.com/somoclu-in-python.html

# In[19]:


# Importing Modules
from sklearn import datasets
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import SVG

kegginfo=pd.read_csv('Kegg_ATC_database_info.csv')
keggdata=pd.read_csv('ATC_KEGG_DRUG_DATABASE.csv')
keggcompl=pd.merge(kegginfo,keggdata,on='KEGG_code')
keggclean=keggcompl.drop_duplicates()

# Loading dataset
df = keggclean
df.data = keggclean[['MW','logP','TPSA','HBD','HBA1','HBA2']]
# Defining Model
model = TSNE(learning_rate=100)

# Fitting Model
transformed = model.fit_transform(df.data)

# Plotting 2d t-Sne
x=[]
y=[]
coord=[]
labels=[]
codes=keggclean['ATC_full_code']
labs=codes.str[:3].reset_index(drop=True)
#labs.reset_index(drop=True)
for i in range(0,len(labs)-1,1):
    x.append(transformed[i, 0])
    y.append(transformed[i, 1])
    coord.append((x, y))
    labels.append(labs[i])

plt.figure(figsize=(40, 40))

for i in range(0,len(labs)-1,1): #get (0, label)
   # x, y = coord[i] #2 dim
    plt.scatter(x[i], y[i])
    plt.annotate(labels[i],xy=(x[i], y[i]),xytext=(5, 2),textcoords='offset points',ha='right',va='bottom')
    
plt.savefig("KEGG_MAP2.png",dpi=500)
plt.savefig('KEGG_MAP2.svg')
plt.show()


# Finally, we can study the properties of all compounds of ATC database in KEGG_DRUG repository (lets observe some of them are duplicated because have more than single classification so they belong at different ATC subgroups). 

# In[20]:


import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 3, constrained_layout=True,figsize=(15,15))

axs[0, 0].hist(keggclean['MW']);
axs[0, 0].set_title('MW')
axs[0, 1].hist(keggclean['logP']);
axs[0, 1].set_title('logP')
axs[0, 2].hist(keggclean['HBA1']);
axs[0, 2].set_title('HBA1')
axs[1, 1].hist(keggclean['HBD']);
axs[1, 1].set_title('HBD')
axs[1, 0].hist(keggclean['TPSA']);
axs[1, 0].set_title('TPSA')
axs[1, 2].hist(keggclean['MR']);
axs[1, 2].set_title('MR')

for ax in axs.flat:
    ax.set(xlabel='', ylabel='')

plt.tight_layout()
#plt.subplots_adjust(left= 0.1, bottom=0.05 , right=0.9, top=0.95, wspace=0.5, hspace=0.5)

plt.show()


# In[21]:


import numpy as np
from scipy import io
from scipy.stats import kurtosis,skew,itemfreq,describe,moment,mode,sem,hmean

data=[]
param='logP'
data=keggclean[param]
print("***********************************")
print("Statistical Summary for "+param+" :")
print("***********************************")
print("mean:          ",round(np.mean(data),2))
print("mode counts:   ",mode(data)[1][0])
print("mode:          ",round(mode(data)[0][0],2))
print("SEM:           ",round(sem(data),2))
print("std:           ",round(np.std(data),2))
print('skew:          ',round(skew(data),2))
print("Kurtosis:      ",round(kurtosis(data),2))
print("max:           ",round(np.max(data),2))
print("min:           ",round(np.min(data),2))
print("Percentile 25 :",round(np.percentile(data,25),2))
print("Percentile 75 :",round(np.percentile(data,75),2))
print("Percentile 10 :",round(np.percentile(data,10),2))
print("Percentile 90 :",round(np.percentile(data,90),2))
print("***********************************")
print("***********************************")


# In[22]:


data=[]
param='MW'
data=keggclean[param]
print("***********************************")
print("Statistical Summary for "+param+" :")
print("***********************************")
print("mean:          ",round(np.mean(data),2))
print("mode counts:   ",mode(data)[1][0])
print("mode:          ",round(mode(data)[0][0],2))
print("SEM:           ",round(sem(data),2))
print("std:           ",round(np.std(data),2))
print('skew:          ',round(skew(data),2))
print("Kurtosis:      ",round(kurtosis(data),2))
print("max:           ",round(np.max(data),2))
print("min:           ",round(np.min(data),2))
print("Percentile 25 :",round(np.percentile(data,25),2))
print("Percentile 75 :",round(np.percentile(data,75),2))
print("Percentile 10 :",round(np.percentile(data,10),2))
print("Percentile 90 :",round(np.percentile(data,90),2))
print("***********************************")
print("***********************************")


# In[23]:


data=[]
param='TPSA'
data=keggclean[param]
print("***********************************")
print("Statistical Summary for "+param+" :")
print("***********************************")
print("mean:          ",round(np.mean(data),2))
print("mode counts:   ",mode(data)[1][0])
print("mode:          ",round(mode(data)[0][0],2))
print("SEM:           ",round(sem(data),2))
print("std:           ",round(np.std(data),2))
print('skew:          ',round(skew(data),2))
print("Kurtosis:      ",round(kurtosis(data),2))
print("max:           ",round(np.max(data),2))
print("min:           ",round(np.min(data),2))
print("Percentile 25 :",round(np.percentile(data,25),2))
print("Percentile 75 :",round(np.percentile(data,75),2))
print("Percentile 10 :",round(np.percentile(data,10),2))
print("Percentile 90 :",round(np.percentile(data,90),2))
print("***********************************")
print("***********************************")


# In[24]:


data=[]
param='HBA1'
data=keggclean[param]
print("***********************************")
print("Statistical Summary for "+param+" :")
print("***********************************")
print("mean:          ",round(np.mean(data),2))
print("mode counts:   ",mode(data)[1][0])
print("mode:          ",round(mode(data)[0][0],2))
print("SEM:           ",round(sem(data),2))
print("std:           ",round(np.std(data),2))
print('skew:          ',round(skew(data),2))
print("Kurtosis:      ",round(kurtosis(data),2))
print("max:           ",round(np.max(data),2))
print("min:           ",round(np.min(data),2))
print("Percentile 25 :",round(np.percentile(data,25),2))
print("Percentile 75 :",round(np.percentile(data,75),2))
print("Percentile 10 :",round(np.percentile(data,10),2))
print("Percentile 90 :",round(np.percentile(data,90),2))
print("***********************************")
print("***********************************")


# In[25]:


data=[]
param='HBD'
data=keggclean[param]
print("***********************************")
print("Statistical Summary for "+param+" :")
print("***********************************")
print("mean:          ",round(np.mean(data),2))
print("mode counts:   ",mode(data)[1][0])
print("mode:          ",round(mode(data)[0][0],2))
print("SEM:           ",round(sem(data),2))
print("std:           ",round(np.std(data),2))
print('skew:          ',round(skew(data),2))
print("Kurtosis:      ",round(kurtosis(data),2))
print("max:           ",round(np.max(data),2))
print("min:           ",round(np.min(data),2))
print("Percentile 25 :",round(np.percentile(data,25),2))
print("Percentile 75 :",round(np.percentile(data,75),2))
print("Percentile 10 :",round(np.percentile(data,10),2))
print("Percentile 90 :",round(np.percentile(data,90),2))
print("***********************************")
print("***********************************")


# Let's verify if Lipinski Rule of five () is satisfied by the average:
# 
# 1. No more than 5 hydrogen bond donors (the total number of nitrogen–hydrogen and oxygen–hydrogen bonds)
# 2. No more than 10 hydrogen bond acceptors (all nitrogen or oxygen atoms)
# 3. A molecular mass less than 500 daltons
# 4. An octanol-water partition coefficient (log P) that does not exceed 5
# 
# https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five

# In[26]:


import numpy as np
round(np.mean(keggclean),2)


# As we can see, condition 1 is satisfied (HBD=2.64 < 5), condition 2 is also satisfied looking at 2 (HBA2 < 10), condition 3 is also satisfied (MW = 377.54 < 500) and also condition 4 (logP = 2.41 < 5) 

# To finish, lets compute and analyze a more sophisticated set of molecular descriptors based on PED (pairwise Energy descriptors: https://www.nature.com/articles/srep43738 ). This descriptors are a compat set that account for rellevant electrostatic patterns in small compounds based on atom partial charges (qi) and spatial coordinates of the atom points (xi,yi,zi) and interactions between them by Electrostatic energy:
# 
# $E_{i,j}=\frac{q_i*q_j}{r_{i,j}}$
# 
# Where $r_{i,j}$ is the Euclidian distance between atom points i and j and qi,qj their respective partial charges. A good method to obtain partial charges is by computing them using MMFF94 method and Openbabel software. 
# 
# Then, using the structure provided by original MOL2 files of LIGANDBOX, the MMFF94 partialcharges are computed:

# In[27]:


get_ipython().system('obabel -imol2 KEGG_DRUG.mol2 -omol2 -O KEGG_DRUG_MMFF.mol2 --partialcharge MMFF94')


# Now lets compute descriptors using the methodology described at https://www.nature.com/articles/srep43738 ) and read the file.

# In[39]:


keggped=pd.read_csv('KEGG_PED_ATC.csv')
keggped.head()


# Try now to make a MAP to see the clusters.

# In[48]:


# Importing Modules
from sklearn import datasets
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import SVG

# Loading dataset
df = keggped
df.data = keggped[['u1','u2','u3','u4','u9','u10','u11','u12']]
# Defining Model
model = TSNE(learning_rate=100)

# Fitting Model
transformed = model.fit_transform(df.data)

# Plotting 2d t-Sne
x=[]
y=[]
coord=[]
labels=[]
codes=keggped['Label']
labs=codes.str[7:9].reset_index(drop=True)
#labs.reset_index(drop=True)
for i in range(0,len(labs)-1,1):
    x.append(transformed[i, 0])
    y.append(transformed[i, 1])
    coord.append((x, y))
    labels.append(labs[i])

plt.figure(figsize=(40, 40))

for i in range(0,len(labs)-1,1): #get (0, label)
   # x, y = coord[i] #2 dim
    plt.scatter(x[i], y[i])
    plt.annotate(labels[i],xy=(x[i], y[i]),xytext=(5, 2),textcoords='offset points',ha='right',va='bottom')
    
plt.savefig("KEGG_MAP_PED_color.png",dpi=500)
plt.savefig('KEGG_MAP_PED_color.svg')



'''
# Plotting 2d t-Sne
x=[]
y=[]
coord=[]
labels=[]

#labs.reset_index(drop=True)
for i in range(0,len(labs)-1,1):
    x.append(transformed[i, 0])
    y.append(transformed[i, 1])

plt.figure(figsize=(40, 40))
plt.scatter(x, y)
    
plt.savefig("KEGG_PED8_MAP.png",dpi=500)
plt.savefig('KEGG_PED8_MAP.svg')
plt.show()'''

plt.show()


# In[49]:


from scipy import io
from scipy.stats import kurtosis,skew,itemfreq,describe,moment,mode,sem,hmean
print("Average of descriptors :",round(np.mean(keggped),2))#,"+/-",round(np.std(keggped),2))


# In[50]:


import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 3, constrained_layout=True,figsize=(15,15))

axs[0, 0].hist(keggped['u1']);
axs[0, 0].set_title('u1')
axs[0, 1].hist(keggped['u2']);
axs[0, 1].set_title('u2')
axs[0, 2].hist(keggped['u3']);
axs[0, 2].set_title('u3')
axs[1, 1].hist(keggped['u10']);
axs[1, 1].set_title('u10')
axs[1, 0].hist(keggped['u11']);
axs[1, 0].set_title('u11')
axs[1, 2].hist(keggped['u12']);
axs[1, 2].set_title('u12')

for ax in axs.flat:
    ax.set(xlabel='', ylabel='')

plt.tight_layout()
#plt.subplots_adjust(left= 0.1, bottom=0.05 , right=0.9, top=0.95, wspace=0.5, hspace=0.5)

plt.show()


# Now, lets take some drug candidates to be effective against COVID19 and try to search the most similar drugs to them in 
# KEGG database. For instance, which compounds in KEGG_DATABASE are more similar to this set of compounds:
# 
# ### Drugs that are currently testing against COVID19
# + From source:https://www.linkedin.com/pulse/covid-19-treatments-clinical-trials-ana-gavald%25C3%25A1/?trackingId=tzMbX9MOQUK8qJMQDKuE3w%3D%3D
# + REMDESIVIR - Antiviral, prefered treatment against COVID-19 (https://en.wikipedia.org/wiki/Remdesivir ), (https://www.foxnews.com/science/remdesivir-what-to-know-about-potential-coronavirus-treatment ) ,( https://www.nature.com/articles/d41573-020-00016-0 )
# + CHLOROQUINE - Antimalarial (https://en.wikipedia.org/wiki/Chloroquine )
# + HYDROXYCHLOROQUINE - Derivative that improves CHLOROQUINE effectiveness (https://www.nature.com/articles/s41421-020-0156-0.pdf )
# + LOPINAVIR - Antiviral against VIH (https://en.wikipedia.org/wiki/Lopinavir )
# + DARUNAVIR - Antiviral against VIH (https://en.wikipedia.org/wiki/Darunavir )
# + FAVIPIRAVIR - Effective with Ebola and other viruses (https://www.livescience.com/flu-drug-could-treat-coronavirus.html ) (https://www.redaccionmedica.com/secciones/sanidad-hoy/coronavirus-tratamiento-china-anuncia-resultados-del-antiviral-favipiravir-4773 )
# + THALIDOMIDE - (https://en.wikipedia.org/wiki/Thalidomide )
# + FINGOLIMOD - Inmunomodulator effect. Multiple sclerosis (MS) threatment (https://en.wikipedia.org/wiki/Fingolimod )
# + GALIDESIVIR - (https://en.wikipedia.org/wiki/Galidesivir )

# Some of this set are not yet in this version of KEGG database. Lets search for compounds similar to
# HYDROXYCHLOROQUINE...Is this compound in the dataset?

# In[72]:


import pandas as pd
kegginfo=pd.read_csv('Kegg_ATC_database_info.csv')
keggdata=pd.read_csv('ATC_KEGG_DRUG_DATABASE.csv')
keggcompl=pd.merge(kegginfo,keggdata,on='KEGG_code')
keggclean=keggcompl.drop_duplicates()
keggmol=keggclean[keggclean['CompoundName'].str[:19]==' Hydroxychloroquine']
keggmol


# Then search in PED KEGG database the KEGG codes: D08050

# In[73]:


keggmol=keggped[keggped['Label'].str[:6]=='D08050']
keggmol


# Now lets search for similar compounds using Manhattan distance metrics:
# 

# In[74]:


S=[]
desf = open('Similarities_'+str(keggmol.iloc[0:1,12:13].values[0][0])+'.shd', 'w')
desf.write('label,S\n')

for j in range(0,np.shape(keggped)[0]-1,1):
#    print(str(j))
    MhD=0.0
    for i in range(0,np.shape(keggped)[1]-2,1):
        MhD=MhD+1/12*np.absolute(keggped.iloc[j:j+1,i:i+1].values[0][0]-keggmol.iloc[0:1,i:i+1].values[0][0])
    S=1.0/(1.0+MhD)
    desf.write(str(keggped.iloc[j:j+1,12:13].values[0][0])+','+str(round(S,4))+'\n')


# In[75]:


desf.close()
df=pd.read_csv('Similarities_'+str(keggmol.iloc[0:1,12:13].values[0][0])+'.shd') #_'+str(keggmol.iloc[0:1,12:13])+'.shd')

dfs=df.sort_values(by='S', ascending=False).reset_index(drop=True)
dfs.head(n=20)


# In[76]:


simfound=[]
for j in range(0,25,1): 
#    str(df.iloc[j:j+1].values[0][0])[:6]
    simfound.append(keggcompl[keggcompl['KEGG_code']==str(dfs.iloc[j:j+1].values[0][0])[:6]])
    
#sb=simfound.reset_index(drop=True)


# In[77]:



swr = open('Candidates.smi', 'w')
for j in range(0,25,1): 
    if(len(simfound[j])>0):
#        print(simfound[j].values[0][5]+' '+simfound[j].values[0][0]+'\n')
        swr.write(simfound[j].values[0][5]+' '+simfound[j].values[0][0]+'\n')
swr.close()


# In[78]:


get_ipython().system('obabel -ismi Candidates.smi -osvg -O Candidates.svg')


# In[79]:


from IPython.display import SVG
display(SVG('Candidates.svg'))


# In[80]:


keggmol=keggclean[keggclean['KEGG_code'].str[:10]=='D02114']
keggmol


# In[81]:


keggmol=keggclean[keggclean['KEGG_code'].str[:10]=='D02243']
keggmol


# In[66]:


keggmol=keggclean[keggclean['CompoundName'].str[:10]==' Lopinavir']
keggmol


# In[67]:


keggmol=keggped[keggped['Label'].str[:6]=='D02498']
keggmol


# In[68]:


S=[]
desf = open('Similarities_'+str(keggmol.iloc[0:1,12:13].values[0][0])+'.shd', 'w')
desf.write('label,S\n')

for j in range(0,np.shape(keggped)[0]-1,1):
#    print(str(j))
    MhD=0.0
    for i in range(0,np.shape(keggped)[1]-2,1):
        MhD=MhD+1/12*np.absolute(keggped.iloc[j:j+1,i:i+1].values[0][0]-keggmol.iloc[0:1,i:i+1].values[0][0])
    S=1.0/(1.0+MhD)
    desf.write(str(keggped.iloc[j:j+1,12:13].values[0][0])+','+str(round(S,4))+'\n')
    
desf.close()
df=pd.read_csv('Similarities_'+str(keggmol.iloc[0:1,12:13].values[0][0])+'.shd') #_'+str(keggmol.iloc[0:1,12:13])+'.shd')

dfs=df.sort_values(by='S', ascending=False).reset_index(drop=True)
dfs.head(n=20)


# In[70]:


simfound=[]
for j in range(0,25,1): 
#    str(df.iloc[j:j+1].values[0][0])[:6]
    simfound.append(keggcompl[keggcompl['KEGG_code']==str(dfs.iloc[j:j+1].values[0][0])[:6]])
    
swr = open('Candidates.smi', 'w')
for j in range(0,25,1): 
    if(len(simfound[j])>0):
#        print(simfound[j].values[0][5]+' '+simfound[j].values[0][0]+'\n')
        swr.write(simfound[j].values[0][5]+' '+simfound[j].values[0][0]+'\n')
swr.close()
get_ipython().system('obabel -ismi Candidates.smi -osvg -O Candidates.svg')
display(SVG('Candidates.svg'))


# The most similar one is D00427. Lets search for that compound in the database:

# In[71]:


keggmol=keggclean[keggclean['KEGG_code'].str[:10]=='D00427']
keggmol


# Ritonavir, a possible compound with similar effects than Lopinavir:
# https://www.nature.com/articles/s41591-020-0849-9
