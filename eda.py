import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# Step 1: Load your data
df = pd.read_csv("amr_multi_organism_kenya.csv")
df.columns = [c.lstrip('#') for c in df.columns]  # cleans the '#' NCBI puts on the first column name

# Step 2: Basic checks - fill these in using what you remember from last time
print(df.shape)
print(df["Organism group"].value_counts())  # how many isolates per organism?
print(df["Isolation source"].value_counts().head(10))

# Step 3: Parse dates into years (same approach as your Klebsiella project)
df['Create date'] = pd.to_datetime(df['Create date'], errors='coerce')
df['Year'] = df['Create date'].dt.year

# Step 4: Compare isolates-per-year ACROSS organisms
# Hint: use groupby() with two columns - Year and Organism group - then unstack for a comparison table
year_org = df.groupby(['Year', 'Organism group']).size().unstack()
print(year_org)

# Step 5: Plot it as a grouped bar chart
year_org.plot(kind='bar', figsize=(10,6))
plt.title("Isolates per year, by organism - Kenya clinical isolates")
plt.ylabel("Number of isolates")
plt.tight_layout()
plt.savefig("chart_organism_trends.png", dpi=150)
def parse_genes(cell):
    if pd.isna(cell):
        return []
    genes = []
    for entry in cell.split(','):
        entry = entry.strip()
        if '=' in entry:
            gene = entry.split('=')[0].strip()  # same logic as your Klebsiella project
        else:
            gene = entry
        if gene:
            genes.append(gene)
    return genes

df['gene_list'] = df['AMR genotypes'].apply(parse_genes)

# For each organism, find its top 10 genes
for organism in df['Organism group'].unique():
    subset = df[df['Organism group'] == organism]
    all_genes = [g for row in subset['gene_list'] for g in row]
    gene_counts = Counter(all_genes)
    print(f"\nTop 10 genes in {organism}:")
    for gene, count in gene_counts.most_common(10):
        print(f"  {gene}: {count}")

# Genes to drug class mapping (acquired resistance genes only —
# intrinsic genes from Day 1 are deliberately excluded)
GENE_TO_CLASS = {
    'blaTEM': 'Beta-lactams',
    'blaCTX-M': 'Beta-lactams',
    'blaSHV': 'Beta-lactams',
    'aph(3': 'Aminoglycosides',   # covers aph(3'')-Ib
    'aph(6)': 'Aminoglycosides',
    'aac(3)': 'Aminoglycosides',
    'aac(6': 'Aminoglycosides',
    'sul1': 'Sulfonamides',
    'sul2': 'Sulfonamides',
    'dfrA': 'Trimethoprim',
    'tet(': 'Tetracyclines',
    'catA': 'Phenicols',
    'qnr': 'Fluoroquinolones',
    'gyrA': 'Fluoroquinolones',
    'gyrB': 'Fluoroquinolones',
    'mcr': 'Colistin',
}

INTRINSIC_GENES = {'mdsA', 'mdsB', 'acrF', 'blaEC', 'mdtM', 'emrD', 'fosA', 'oqxA', 'oqxB'}

def get_drug_classes(gene_list):
    classes = set()
    for gene in gene_list:
        if gene in INTRINSIC_GENES:
            continue
        for prefix, drug_class in GENE_TO_CLASS.items():
            if gene.startswith(prefix):
                classes.add(drug_class)
    return classes

df['drug_classes'] = df['gene_list'].apply(get_drug_classes)
df['n_classes'] = df['drug_classes'].apply(len)
df['MDR'] = df['n_classes'] >= 3

print(df['MDR'].value_counts())
print(df.groupby('Organism group')['MDR'].mean())

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
import pandas as pd

# Step 1: Select your features and target
features = df[['Organism group', 'Isolation source', 'Year']].copy()
target = df['MDR']

# Step 2: One-hot encode the categorical columns
# Hint: pd.get_dummies() converts text columns into 0/1 numeric columns automatically
features_encoded = pd.get_dummies(features, columns=['Organism group', 'Isolation source'])

# Step 3: Handle any missing Year values (some rows might have NaN)
features_encoded['Year'] = features_encoded['Year'].fillna(features_encoded['Year'].median())

# Step 4: Split into train and test sets
# Hint: train_test_split() takes your features, your target, and a test_size (as a decimal, e.g. 20%)
X_train, X_test, y_train, y_test = train_test_split(
    features_encoded, target, test_size=0.2, random_state=42
)

# Step 5: Train a Random Forest model
model = RandomForestClassifier(random_state=42)
model.fit(X_train, y_train)   # hint: the method that trains/fits a model on data

# Step 6: Predict on the test set and evaluate
predictions = model.predict(X_test)   # hint: the method that generates predictions
print(confusion_matrix(y_test, predictions))
print(classification_report(y_test, predictions))
importances = pd.Series(model.feature_importances_, index=X_train.columns)
print(importances.sort_values(ascending=False).head(10))

import networkx as nx
from itertools import combinations
from collections import Counter
import matplotlib.pyplot as plt

# Step 1: Only look at acquired genes (exclude intrinsic ones from Day 1)
def get_acquired_genes(gene_list):
    return [g for g in gene_list if g not in INTRINSIC_GENES]

df['acquired_genes'] = df['gene_list'].apply(get_acquired_genes)

# Step 2: Count how often every PAIR of genes appears together
pair_counts = Counter()
for genes in df['acquired_genes']:
    # Only count genes that actually co-occur (2 or more in the same isolate)
    unique_genes = sorted(set(genes))  # hint: sorting makes pairs consistent (A,B) not (B,A)
    for pair in combinations(unique_genes, 2):
        pair_counts[pair] += 1

# Step 3: Keep only strong pairs (appearing together often) to avoid a messy graph
MIN_COOCCURRENCE = 1000  # try adjusting this number based on your dataset size
strong_pairs = {pair: count for pair, count in pair_counts.items() if count >= MIN_COOCCURRENCE}

print(f"Number of strong gene pairs: {len(strong_pairs)}")

# Step 4: Build the network graph
G = nx.Graph()
for (gene1, gene2), count in strong_pairs.items():
    G.add_edge(gene1, gene2, weight=count)

# Step 5: Draw it
plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G, k=0.5, seed=42)
nx.draw(G, pos, with_labels=True, node_color='lightgreen', 
        node_size=1500, font_size=8, edge_color='gray')
plt.title("Resistance Gene Co-occurrence Network")
plt.savefig("gene_network.png", dpi=150)
df.to_csv("processed_data.csv", index=False)
