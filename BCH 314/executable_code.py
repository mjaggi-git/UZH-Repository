# adds the paths to the repo folder
import sys
sys.path.insert(1, '../preprocessing/') # preprocessing code
sys.path.insert(1, '../experiments/') # folder with configuration files
sys.path.insert(1, '../data/') # data code
sys.path.insert(1, '../model/') # model code
sys.path.insert(1, '../evaluation_BIMODAL/') # model code
import os

# repo modules
from main_preprocessor import preprocess_data #for data preprocessing
import configparser # to automatically change the .ini file
from sample import Sampler

# other modules
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw # for molecule depiction
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Scaffolds import MurckoScaffold
import ChemTools
%load_ext autoreload
%autoreload 2

model = 'BIMODAL'
start = 'fixed' 
augmentation_level = 1
network_size = 512

# training process
fine_tuning = 'actives_format_1000' 
epochs = 10

# sampling process
T_sampling = 0.7 
n_sampling = 2000

# Uncomment the following line to execute the pretreatment again. Attention: this will overwrite the existing file!
preprocess_data(filename_in='actives_format_1000', filename_out='../data/actives_format_1000_BIMODAL_fixed_512_aug_1.csv', model_type=model, starting_point=start, augmentation=augmentation_level)

input_csv = '../data/actives_format_1000_BIMODAL_fixed_512_aug_1.csv'
df = pd.read_csv(input_csv, header=None)

# Filter SMILES that don't end with either A or E
filter_condition = df[0].astype(str).str.endswith(('E', 'A'))
filtered_df = df[filter_condition].copy()
print(f"Removed {len(df) - len(filtered_df)} SMILES not ending with E/A")

# Padding to molecular size of 151
padding = 'A' * 45
df[0] = df[0].astype(str).apply(lambda x: padding + x + padding)
output_csv = '../data/actives_format_1000_BIMODAL_fixed_512_aug_1.csv'
df.to_csv(output_csv, header=False, index=False)

filepath = '../data/'+ fine_tuning + '_' + model + '_' + start + '_' + str(network_size) + '_' + 'aug' + '_' + str(augmentation_level) + '.csv' #path to the pretreated file based on the settings
ft_set = pd.read_csv(filepath,header=None) 
ft_set # SMILES overview

import utils
exp_name = utils.make_config(model_type=model, net_size=network_size, epochs=epochs, starting_point=start, fine_tuning='actives_format_1000', n_sampling=n_sampling,T_sampling=T_sampling,augmentation_level=augmentation_level)

import fine_tuner
t = fine_tuner.FineTuner(experiment_name = exp_name)
#Uncomment the following line to start the fine-tuning. Attention, results will be overwritten!
t.fine_tuning(stor_dir='../evaluation_BIMODAL/')

### ANALYSIS PART ###

# Entropy loss graph
import matplotlib.pyplot as plt

# Load data
try:
    loss = pd.read_csv(f'../evaluation_BIMODAL/BIMODAL_fixed_512_FineTuning/statistic/stat_fold_1.csv', header=None)
    print(f"Loaded data from: ../evaluation_BIMODAL/BIMODAL_fixed_512_FineTuning/statistic/stat_fold_1.csv")
except FileNotFoundError:
    raise SystemExit("Error: CSV file not found. Check the path or 'exp_name'.")
    
# Create plot 
loss.columns = ["Epochs", "Cross Entropy Loss"]
fig, ax = plt.subplots(figsize=(6, 6))  # Explicit figure/axis creation
loss.plot.line(x="Epochs", y="Cross Entropy Loss", ax=ax)
ax.set_title("Training Loss Curve")

# Save plot file
plt.savefig('../evaluation_BIMODAL/results/BIMODAL_fixed_512_aug_1_cross_entropy_loss.png', bbox_inches='tight', dpi=300)
plt.show()
plt.close() 

from ChemTools import novelty_vs_epochs
folder = '../evaluation_BIMODAL/BIMODAL_fixed_512_FineTuning/' # generates the folder path automatically
stats = novelty_vs_epochs(folder,ft_path='../example/actives_format_1000.csv',export=True)

results_path = '../evaluation_BIMODAL/BIMODAL_fixed_512_FineTuning/analysis'
if os.path.exists(results_path) is False: # checks if the folder exists, and otherwise creates it
    os.mkdir(results_path)

stats.to_csv(results_path + '/novelty_results.csv')

ax = stats.plot.bar(figsize=(10, 8), rot=45, title="Novelty/Validity/Uniqueness Results")
ax.set_xlabel("Epochs")
ax.set_ylabel("Molecules")
plt.tight_layout()  # Prevent label overlap
plt.savefig('../evaluation_BIMODAL/results/BIMODAL_fixed_512_aug_1_novelty_validity_uniqueness.png', bbox_inches='tight', dpi=300)
plt.show()
plt.close() 

from ChemTools import scaffolds_vs_epochs
folder = '../evaluation_BIMODAL/BIMODAL_fixed_512_FineTuning/molecules_novel/'
scaffold_res = scaffolds_vs_epochs(folder)

# saves the results
scaffold_res.to_csv(os.path.join(results_path, 'scaffolds_results.csv'))

# creates and saves the plot
fig, ax = plt.subplots(figsize=(10, 8))
scaffold_res.plot.bar(ax=ax, subplots=True)

plt.tight_layout()  # prevents label overlapping
plt.savefig(
    os.path.join('../evaluation_BIMODAL/results', 'BIMODAL_fixed_512_aug_1_scaffold_analysis.png'),
    bbox_inches='tight',
    dpi=300
)
plt.show()
plt.close()

from ChemTools import morgan_vs_epochs
folder = '../evaluation_BIMODAL/BIMODAL_fixed_512_FineTuning/molecules_novel/'
sim_res = ChemTools.morgan_vs_epochs(folder,path_ft='../example/actives_format_1000.csv')

import seaborn as sns

# saves the results
sim_res.to_csv(os.path.join(results_path, 'similarity_results.csv'))

# creates and displays the plot
plt.figure(figsize=(10, 6))  # Added figure size for better proportions
ax = sns.lineplot(data=sim_res, ci=1000)

# Set plot title and labels
ax.set_title("Morgan Similarity Analysis", fontsize=14)
ax.set_xlabel("Epochs", fontsize=12)
ax.set_ylabel("Similarity", fontsize=12)

# Save the plot as PNG
plt.savefig(
    os.path.join('../evaluation_BIMODAL/results', 'BIMODAL_fixed_512_aug_1_morgan_similarity.png'),
    bbox_inches='tight',
    dpi=300
)

# Display and close
plt.show()
plt.close()

epoch = 10
file_path = '../evaluation_BIMODAL/BIMODAL_fixed_512_FineTuning/molecules_novel/' + 'molecule_fold_1_epochs_' + str(epoch-1) + '.csv' # generates the filename
designs = Chem.SmilesMolSupplier(file_path,smilesColumn=0,titleLine=False,delimiter=',')

from ChemTools import frequent_scaffolds
from rdkit.Chem import Draw
import os

freq_scaffolds = frequent_scaffolds(designs)

k = 6  # number of scaffolds to depict

# Create the grid image
img = Draw.MolsToGridImage(
    freq_scaffolds[:k],
    molsPerRow=2,
    subImgSize=(200, 200),
    legends=[x.GetProp("_Name") for x in freq_scaffolds[:k]]
)

# Save as PNG
output_path = '../evaluation_BIMODAL/results/BIMODAL_fixed_512_aug_1_top_scaffolds.png'
os.makedirs(os.path.dirname(output_path), exist_ok=True)  # Ensure directory exists
img.save(output_path)

# Display the image (optional)
img.show()
