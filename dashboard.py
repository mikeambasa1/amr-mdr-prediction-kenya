import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

st.title("AMR Surveillance Dashboard — Kenya")
st.write("Explore antimicrobial resistance patterns across clinical isolates in Kenya")
INTRINSIC_GENES = {'mdsA', 'mdsB', 'acrF', 'blaEC', 'mdtM', 'emrD', 'fosA', 'oqxA', 'oqxB'}
# Load your data (reuse your cleaned dataframe logic from eda.py)
df = pd.read_csv("processed_data.csv")
import ast
df['gene_list'] = df['gene_list'].apply(ast.literal_eval)
df.columns = [c.lstrip('#') for c in df.columns]
df['Create date'] = pd.to_datetime(df['Create date'], errors='coerce')
df['Year'] = df['Create date'].dt.year

# --- Sidebar filters ---
st.sidebar.header("Filters")

# Hint: st.sidebar.multiselect(label, options_list, default_value) lets users pick multiple organisms
organisms = st.sidebar.multiselect("Select organism(s)", df['Organism group'].unique(), default=df['Organism group'].unique())

# Hint: st.sidebar.slider(label, min_value, max_value) lets users pick a year range
year_range = st.sidebar.slider("Select year range", int(df['Year'].min()), int(df['Year'].max()), (2022, 2026))

# --- Apply filters ---
filtered = df[
    (df['Organism group'].isin(organisms)) &
    (df['Year'] >= year_range[0]) &
    (df['Year'] <= year_range[1])
]

st.write(f"Showing {len(filtered)} isolates")

# --- Chart 1: isolates per year, filtered ---
st.subheader("Isolates per year")
year_counts = filtered['Year'].value_counts().sort_index()
# Hint: st.bar_chart() takes a pandas Series or DataFrame directly - no matplotlib needed
st.bar_chart(year_counts)

# --- Chart 2: top isolation sources, filtered ---
st.subheader("Top isolation sources")
source_counts = filtered['Isolation source'].value_counts().head(10)
st.bar_chart(source_counts)
# --- MDR rate for current selection ---
st.subheader("MDR Rate for Selection")
mdr_rate = filtered['MDR'].mean() * 100
st.metric("Percent MDR", f"{mdr_rate:.1f}%")

# --- Top resistance genes, filtered ---
st.subheader("Top resistance genes (excluding intrinsic genes)")
all_genes = [g for row in filtered['gene_list'] for g in row if g not in INTRINSIC_GENES]
gene_counts = Counter(all_genes)
top = pd.Series(dict(gene_counts.most_common(10)))
st.bar_chart(top)