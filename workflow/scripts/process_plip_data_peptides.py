# Script to parse PLIP data
#
# Need to extract the following information from multiple PLIP output files (XML)
# - Hydrophobic contacts (`hydrophobic_interactions`)
# - Pi stacking (`pi_stacks`)
# - Hydrogen bonds (`hydrogen_bonds`)
# - Salt bridges (`salt_bridges`)
# - Cation-pi interactions (`pi_cation_interactions`)
#
# Updated this script on 14/01/2022 to also get the aggregated data for
# individual contacts (e.g. hydrogen bonds, salt bridges etc.), basically
# merging the hydrophobic contacts script with this!
#
# Adapted for use with snakemake
# Snakemake input is a list of extractedXXX pdb folders, so just need to find
# report.xml in these folders and add to the file list

import xml.etree.ElementTree as ET
import pandas as pd
import os

def get_interaction(et_root, interaction_type):
    '''
    General function for searching XML file for interactions, used for
    hydrophilicity score calculation

    For peptides rather than small molecule ligands here, so need lig_res_no

    Not extracting distances as a result, just need the residue number etc.
    '''

    path_int = f'bindingsite/interactions/{interaction_type}/*'
    data = []

    for interaction in et_root.findall(path_int):
        res_no = interaction.find('resnr').text
        res_type = interaction.find('restype').text
        lig_res_no = interaction.find('resnr_lig').text
        lig_res_type = interaction.find('restype_lig').text

        combined = [res_type+res_no, lig_res_type+lig_res_no, float(lig_res_no),
                    interaction_type]
        data.append(combined)

    return data

file_list = []
hydrophobic_contacts = []
hydrophilic_contacts = []

for f in snakemake.input:
    file_list.append(os.path.join(f,'report.xml'))

print(file_list[0:5])

# extract contacts
for frame, file in enumerate(file_list):
    tree = ET.parse(file)
    root = tree.getroot()

    # extract hydrophobic contacts
    hydrophobic_interactions = get_interaction(
        root, 'hydrophobic_interactions')
    pi_stacks = get_interaction(root, 'pi_stacks')
    for i in hydrophobic_interactions:
        i.insert(0, frame)
        hydrophobic_contacts.append(i)
    for i in pi_stacks:
        # in reality, there are no pi_stacks
        i.insert(0, frame)
        hydrophobic_contacts.append(i)

    # extract hydrophilic contacts
    h_bonds = get_interaction(root, 'hydrogen_bonds')
    salt_bridges = get_interaction(root, 'salt_bridges')
    cat_pi = get_interaction(root, 'pi_cation_interactions')
    for i in h_bonds:
        i.insert(0, frame)
        hydrophilic_contacts.append(i)
    for i in salt_bridges:
        i.insert(0, frame)
        hydrophilic_contacts.append(i)
    for i in cat_pi:
        i.insert(0, frame)
        hydrophilic_contacts.append(i)

cols = ['Frame', 'Protein_ID', 'Peptide_ID', 'Peptide_number', 'Type']
df_hydrophobic = pd.DataFrame(hydrophobic_contacts, columns=cols)
df_hydrophilic = pd.DataFrame(hydrophilic_contacts, columns=cols)

print(f'Total hydrophobic contacts: {len(hydrophobic_contacts)}')
print(f'Total hydrophilic contacts: {len(hydrophilic_contacts)}')

# Make dataframes with different breakdowns
hydrophobic_counts = df_hydrophobic.groupby(['Peptide_number']).size()
hydrophobic_counts = hydrophobic_counts.reindex(list(range(1,22)), fill_value=0)
hydrophilic_counts = df_hydrophilic.groupby(['Peptide_number']).size()
hydrophilic_counts = hydrophilic_counts.reindex(list(range(1, 22)), fill_value=0)

combined_df = pd.concat([hydrophobic_counts, hydrophilic_counts], axis=1)
combined_df.columns = ['Hydrophobic', 'Hydrophilic']

hydrophobic_breakdown = df_hydrophobic.groupby(['Peptide_number','Peptide_ID', 'Protein_ID', 'Type'], as_index=False).size()
hydrophilic_breakdown = df_hydrophilic.groupby(['Peptide_number', 'Peptide_ID', 'Protein_ID', 'Type'], as_index=False).size()

df_breakdown = pd.concat([hydrophobic_breakdown, hydrophilic_breakdown])
df_breakdown.fillna(value=0, inplace=True)
df_breakdown.sort_values(by='Peptide_number', inplace=True)

hbonds = hydrophilic_breakdown[hydrophilic_breakdown['Type']=='hydrogen_bonds']
hydrophobic = hydrophobic_breakdown[hydrophobic_breakdown['Type']=='hydrophobic_interactions']

with pd.ExcelWriter(snakemake.output[0]) as writer:
    combined_df.to_excel(writer, sheet_name='scores')
    df_breakdown.to_excel(writer, sheet_name='breakdown')
    hbonds.to_excel(writer, sheet_name='hbonds_only')
    hydrophobic.to_excel(writer, sheet_name='hydrophobic_only')
