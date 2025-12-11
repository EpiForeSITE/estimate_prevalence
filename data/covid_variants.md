# COVID-19 Variants Over Time


## Overview

This document processes phylogenetic data from a Newick file containing
COVID-19 variant sequences from North America. The goal is to extract
the date/time information from each variant and count the number of
distinct variants observed at each time point.

## Load Libraries

``` python
from Bio import Phylo
import pandas as pd
import matplotlib.pyplot as plt
import re
import io
```

## Read and Parse the Phylogenetic Tree

``` python
# Read the Newick file
with open('nextstrain_ncov_gisaid_north-america_6m_timetree.nwk', 'r') as f:
    content = f.read()

# Parse the tree
tree = Phylo.read(io.StringIO(content), 'newick')

# Get all terminal nodes (leaves) - these are the virus variants
terminals = tree.get_terminals()
print(f"Total number of variants in the tree: {len(terminals)}")
```

    Total number of variants in the tree: 2869

## Extract Variant Information

The variant names follow the format:
`hCoV-19/{country}/{identifier}/{year}`

We will extract the year from each variant name.

``` python
# Extract variant names and years
variants_data = []

for terminal in terminals:
    name = terminal.name
    if name:
        # Extract year from the variant name (format: .../YYYY at the end)
        match = re.search(r'/(\d{4})$', name)
        if match:
            year = int(match.group(1))
            # Filter to valid COVID-19 pandemic years (2019-2025)
            if 2019 <= year <= 2025:
                variants_data.append({
                    'variant_name': name,
                    'year': year
                })

print(f"Variants with valid year data: {len(variants_data)}")
```

    Variants with valid year data: 2869

## Create Dataset by Year

``` python
# Create a DataFrame
df = pd.DataFrame(variants_data)

# Count variants by year
variants_by_year = df.groupby('year').size().reset_index(name='variant_count')
variants_by_year.columns = ['year', 'variant_count']

print("Variants by Year:")
print(variants_by_year.to_string(index=False))
```

    Variants by Year:
     year  variant_count
     2020            105
     2021            118
     2022            132
     2023            123
     2024            116
     2025           2275

## Cumulative Variants Over Time

``` python
# Calculate cumulative count
variants_by_year['cumulative_count'] = variants_by_year['variant_count'].cumsum()

print("\nVariants by Year with Cumulative Count:")
print(variants_by_year.to_string(index=False))
```


    Variants by Year with Cumulative Count:
     year  variant_count  cumulative_count
     2020            105               105
     2021            118               223
     2022            132               355
     2023            123               478
     2024            116               594
     2025           2275              2869

## Visualization

``` python
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Bar chart of variants per year
axes[0].bar(variants_by_year['year'], variants_by_year['variant_count'], color='steelblue')
axes[0].set_xlabel('Year')
axes[0].set_ylabel('Number of Variants')
axes[0].set_title('COVID-19 Variants per Year')
axes[0].set_xticks(variants_by_year['year'])

# Line chart of cumulative variants
axes[1].plot(variants_by_year['year'], variants_by_year['cumulative_count'], 
             marker='o', linewidth=2, color='darkgreen')
axes[1].set_xlabel('Year')
axes[1].set_ylabel('Cumulative Number of Variants')
axes[1].set_title('Cumulative COVID-19 Variants Over Time')
axes[1].set_xticks(variants_by_year['year'])
axes[1].fill_between(variants_by_year['year'], variants_by_year['cumulative_count'], 
                     alpha=0.3, color='green')

plt.tight_layout()
plt.savefig('covid_variants_plot.png', dpi=150, bbox_inches='tight')
plt.show()
```

<div id="fig-variants-by-year">

![](covid_variants_files/figure-commonmark/fig-variants-by-year-output-1.png)

Figure 1: Number of COVID-19 variants observed by year

</div>

## Save Data to CSV

``` python
# Save the dataset to CSV
variants_by_year.to_csv('covid_variants.csv', index=False)
print(f"Data saved to covid_variants.csv")
print(f"\nCSV Contents:")
print(variants_by_year.to_string(index=False))
```

    Data saved to covid_variants.csv

    CSV Contents:
     year  variant_count  cumulative_count
     2020            105               105
     2021            118               223
     2022            132               355
     2023            123               478
     2024            116               594
     2025           2275              2869

## Summary Statistics

``` python
print(f"\n=== Summary Statistics ===")
print(f"Total variants analyzed: {len(df)}")
print(f"Year range: {df['year'].min()} - {df['year'].max()}")
print(f"Year with most variants: {variants_by_year.loc[variants_by_year['variant_count'].idxmax(), 'year']} ({variants_by_year['variant_count'].max()} variants)")
print(f"Average variants per year: {variants_by_year['variant_count'].mean():.1f}")
```


    === Summary Statistics ===
    Total variants analyzed: 2869
    Year range: 2020 - 2025
    Year with most variants: 2025 (2275 variants)
    Average variants per year: 478.2

## Raw Data Sample

``` python
print("\nSample of variant names (first 10):")
for i, row in df.head(10).iterrows():
    print(f"  - {row['variant_name']} ({row['year']})")
```


    Sample of variant names (first 10):
      - hCoV-19/USA/Cruise-CDC-3054899-001/2020 (2020)
      - hCoV-19/USA/NY1-PV08001/2020 (2020)
      - hCoV-19/USA/CA-CDC-02993471-001/2020 (2020)
      - hCoV-19/USA/TX-HMH-1788/2020 (2020)
      - hCoV-19/USA/MA-MGH-00614/2020 (2020)
      - hCoV-19/Mexico/CMX-INER-IBT-4/2020 (2020)
      - hCoV-19/USA/WA-S2/2020 (2020)
      - hCoV-19/USA/WA-S3/2020 (2020)
      - hCoV-19/USA/MD-HP09649-PIDXSVKTIH/2020 (2020)
      - hCoV-19/USA/WA-S174/2020 (2020)
