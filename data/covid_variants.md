# COVID-19 Variants Over Time


## Overview

This document processes phylogenetic data from a Newick file containing
COVID-19 variant sequences from North America. The goal is to extract
the date/time information from each variant and count the number of
distinct variants observed at each time point.

The Newick file is a timetree where branch lengths are in units of
years. By calculating the distance from the root to each terminal node,
we can determine the date of each variant sample.

## Load Libraries

``` python
from Bio import Phylo
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
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

## Convert Branch Lengths to Dates

The timetree has branch lengths in units of years. The root represents
the most recent common ancestor of SARS-CoV-2. Based on the tree
structure, we set the root date to December 2019 (the beginning of the
COVID-19 pandemic).

``` python
def decimal_year_to_date(decimal_year):
    """Convert a decimal year value to a date.
    
    For example, 2025.5 would be approximately June 2025 (mid-year).
    """
    year = int(decimal_year)
    remainder = decimal_year - year
    base_date = datetime(year, 1, 1)
    days_in_year = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
    result_date = base_date + timedelta(days=remainder * days_in_year)
    return result_date

# Root date for the timetree (start of COVID-19 pandemic)
# Based on the tree structure, the root corresponds to approximately December 2019
ROOT_DATE_DECIMAL = 2019.95  # Approximately mid-December 2019

# Extract variant data with dates calculated from branch lengths
variants_data = []

for terminal in terminals:
    name = terminal.name
    if name:
        # Calculate distance from root (sum of all branch lengths to this node)
        distance_from_root = tree.distance(terminal)
        
        # Convert to decimal year by adding distance to root date
        decimal_year = ROOT_DATE_DECIMAL + distance_from_root
        
        # Convert to actual date
        sample_date = decimal_year_to_date(decimal_year)
        
        variants_data.append({
            'variant_name': name,
            'decimal_year': round(decimal_year, 4),
            'date': sample_date.strftime('%Y-%m-%d'),
            'year': sample_date.year,
            'month': sample_date.month
        })

print(f"Total variants processed: {len(variants_data)}")
```

    Total variants processed: 2869

## Examine the Date Distribution

``` python
# Create a DataFrame
df = pd.DataFrame(variants_data)
df['date'] = pd.to_datetime(df['date'])

print("Date range of samples:")
print(f"  Earliest: {df['date'].min().strftime('%Y-%m-%d')}")
print(f"  Latest: {df['date'].max().strftime('%Y-%m-%d')}")

print("\nSample of variant data (first 10):")
print(df[['variant_name', 'decimal_year', 'date']].head(10).to_string(index=False))
```

    Date range of samples:
      Earliest: 2020-01-09
      Latest: 2025-09-04

    Sample of variant data (first 10):
                               variant_name  decimal_year       date
    hCoV-19/USA/Cruise-CDC-3054899-001/2020      2020.107 2020-02-09
               hCoV-19/USA/NY1-PV08001/2020      2020.129 2020-02-17
       hCoV-19/USA/CA-CDC-02993471-001/2020      2020.027 2020-01-10
               hCoV-19/USA/TX-HMH-1788/2020      2020.328 2020-04-30
              hCoV-19/USA/MA-MGH-00614/2020      2020.205 2020-03-16
         hCoV-19/Mexico/CMX-INER-IBT-4/2020      2020.208 2020-03-17
                     hCoV-19/USA/WA-S2/2020      2020.107 2020-02-09
                     hCoV-19/USA/WA-S3/2020      2020.126 2020-02-16
     hCoV-19/USA/MD-HP09649-PIDXSVKTIH/2020      2020.265 2020-04-06
                   hCoV-19/USA/WA-S174/2020      2020.180 2020-03-06

## Create Dataset by Month

For more granular analysis, we’ll count variants by month.

``` python
# Create year-month column for grouping
df['year_month'] = df['date'].dt.to_period('M')

# Count variants by month
variants_by_month = df.groupby('year_month').size().reset_index(name='variant_count')
variants_by_month['year_month'] = variants_by_month['year_month'].astype(str)

# Calculate cumulative count
variants_by_month['cumulative_count'] = variants_by_month['variant_count'].cumsum()

print("Variants by Month (showing first 20 rows):")
print(variants_by_month.head(20).to_string(index=False))
```

    Variants by Month (showing first 20 rows):
    year_month  variant_count  cumulative_count
       2020-01              2                 2
       2020-02              6                 8
       2020-03             16                24
       2020-04              8                32
       2020-05              8                40
       2020-06             13                53
       2020-07             11                64
       2020-08              6                70
       2020-09              8                78
       2020-10             12                90
       2020-11              9                99
       2020-12             10               109
       2021-01             10               119
       2021-02              8               127
       2021-03             10               137
       2021-04             10               147
       2021-05             12               159
       2021-06              4               163
       2021-07             10               173
       2021-08             11               184

## Create Dataset by Year

``` python
# Count variants by year
variants_by_year = df.groupby('year').size().reset_index(name='variant_count')
variants_by_year.columns = ['year', 'variant_count']

# Calculate cumulative count
variants_by_year['cumulative_count'] = variants_by_year['variant_count'].cumsum()

print("\nVariants by Year:")
print(variants_by_year.to_string(index=False))
```


    Variants by Year:
     year  variant_count  cumulative_count
     2020            109               109
     2021            118               227
     2022            136               363
     2023            120               483
     2024            114               597
     2025           2272              2869

## Visualization

``` python
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Bar chart of variants per year
axes[0, 0].bar(variants_by_year['year'], variants_by_year['variant_count'], color='steelblue')
axes[0, 0].set_xlabel('Year')
axes[0, 0].set_ylabel('Number of Variants')
axes[0, 0].set_title('COVID-19 Variants per Year')
axes[0, 0].set_xticks(variants_by_year['year'])

# Line chart of cumulative variants by year
axes[0, 1].plot(variants_by_year['year'], variants_by_year['cumulative_count'], 
             marker='o', linewidth=2, color='darkgreen')
axes[0, 1].set_xlabel('Year')
axes[0, 1].set_ylabel('Cumulative Number of Variants')
axes[0, 1].set_title('Cumulative COVID-19 Variants Over Time (by Year)')
axes[0, 1].set_xticks(variants_by_year['year'])
axes[0, 1].fill_between(variants_by_year['year'], variants_by_year['cumulative_count'], 
                     alpha=0.3, color='green')

# Monthly variant counts
month_dates = pd.to_datetime(variants_by_month['year_month'])
axes[1, 0].bar(range(len(month_dates)), variants_by_month['variant_count'], color='coral', width=0.8)
axes[1, 0].set_xlabel('Time')
axes[1, 0].set_ylabel('Number of Variants')
axes[1, 0].set_title('COVID-19 Variants per Month')
# Set x-axis labels for every 6 months
tick_positions = range(0, len(month_dates), 6)
tick_labels = [month_dates.iloc[i].strftime('%Y-%m') if i < len(month_dates) else '' for i in tick_positions]
axes[1, 0].set_xticks(tick_positions)
axes[1, 0].set_xticklabels(tick_labels, rotation=45, ha='right')

# Cumulative variants by month
axes[1, 1].plot(range(len(month_dates)), variants_by_month['cumulative_count'], 
                linewidth=2, color='purple')
axes[1, 1].set_xlabel('Time')
axes[1, 1].set_ylabel('Cumulative Number of Variants')
axes[1, 1].set_title('Cumulative COVID-19 Variants Over Time (by Month)')
axes[1, 1].set_xticks(tick_positions)
axes[1, 1].set_xticklabels(tick_labels, rotation=45, ha='right')
axes[1, 1].fill_between(range(len(month_dates)), variants_by_month['cumulative_count'], 
                        alpha=0.3, color='purple')

plt.tight_layout()
plt.savefig('covid_variants_plot.png', dpi=150, bbox_inches='tight')
plt.show()
```

<div id="fig-variants-over-time">

![](covid_variants_files/figure-commonmark/fig-variants-over-time-output-1.png)

Figure 1: Number of COVID-19 variants observed over time

</div>

## Save Data to CSV

We save both the monthly and yearly aggregated data.

``` python
# Save the monthly dataset to CSV (more granular - primary output)
variants_by_month.to_csv('covid_variants.csv', index=False)
print(f"Monthly data saved to covid_variants.csv")

# Also save yearly summary
variants_by_year.to_csv('covid_variants_yearly.csv', index=False)
print(f"Yearly data saved to covid_variants_yearly.csv")

print(f"\nCSV Contents (covid_variants.csv - first 20 rows):")
print(variants_by_month.head(20).to_string(index=False))
```

    Monthly data saved to covid_variants.csv
    Yearly data saved to covid_variants_yearly.csv

    CSV Contents (covid_variants.csv - first 20 rows):
    year_month  variant_count  cumulative_count
       2020-01              2                 2
       2020-02              6                 8
       2020-03             16                24
       2020-04              8                32
       2020-05              8                40
       2020-06             13                53
       2020-07             11                64
       2020-08              6                70
       2020-09              8                78
       2020-10             12                90
       2020-11              9                99
       2020-12             10               109
       2021-01             10               119
       2021-02              8               127
       2021-03             10               137
       2021-04             10               147
       2021-05             12               159
       2021-06              4               163
       2021-07             10               173
       2021-08             11               184

## Summary Statistics

``` python
print(f"\n=== Summary Statistics ===")
print(f"Total variants analyzed: {len(df)}")
print(f"Date range: {df['date'].min().strftime('%Y-%m-%d')} to {df['date'].max().strftime('%Y-%m-%d')}")
print(f"Number of months covered: {len(variants_by_month)}")
print(f"Month with most variants: {variants_by_month.loc[variants_by_month['variant_count'].idxmax(), 'year_month']} ({variants_by_month['variant_count'].max()} variants)")
print(f"Average variants per month: {variants_by_month['variant_count'].mean():.1f}")
```


    === Summary Statistics ===
    Total variants analyzed: 2869
    Date range: 2020-01-09 to 2025-09-04
    Number of months covered: 69
    Month with most variants: 2025-04 (574 variants)
    Average variants per month: 41.6

## Raw Data Sample

``` python
print("\nSample of variant data with dates (first 15):")
for i, row in df.head(15).iterrows():
    print(f"  - {row['variant_name']}")
    print(f"      Date: {row['date'].strftime('%Y-%m-%d')} (decimal year: {row['decimal_year']})")
```


    Sample of variant data with dates (first 15):
      - hCoV-19/USA/Cruise-CDC-3054899-001/2020
          Date: 2020-02-09 (decimal year: 2020.107)
      - hCoV-19/USA/NY1-PV08001/2020
          Date: 2020-02-17 (decimal year: 2020.129)
      - hCoV-19/USA/CA-CDC-02993471-001/2020
          Date: 2020-01-10 (decimal year: 2020.027)
      - hCoV-19/USA/TX-HMH-1788/2020
          Date: 2020-04-30 (decimal year: 2020.328)
      - hCoV-19/USA/MA-MGH-00614/2020
          Date: 2020-03-16 (decimal year: 2020.205)
      - hCoV-19/Mexico/CMX-INER-IBT-4/2020
          Date: 2020-03-17 (decimal year: 2020.208)
      - hCoV-19/USA/WA-S2/2020
          Date: 2020-02-09 (decimal year: 2020.107)
      - hCoV-19/USA/WA-S3/2020
          Date: 2020-02-16 (decimal year: 2020.126)
      - hCoV-19/USA/MD-HP09649-PIDXSVKTIH/2020
          Date: 2020-04-06 (decimal year: 2020.265)
      - hCoV-19/USA/WA-S174/2020
          Date: 2020-03-06 (decimal year: 2020.18)
      - hCoV-19/USA/IL-CDC-02983522-001/2020
          Date: 2020-01-09 (decimal year: 2020.022)
      - hCoV-19/USA/CA-SCCPHD-UC103/2020
          Date: 2020-02-16 (decimal year: 2020.126)
      - hCoV-19/Canada/BC-BCCDC-166226/2020
          Date: 2020-10-16 (decimal year: 2020.79)
      - hCoV-19/USA/IL-RUSH-00555/2020
          Date: 2020-03-22 (decimal year: 2020.224)
      - hCoV-19/Cuba/USAFSAM-S031/2020
          Date: 2020-03-12 (decimal year: 2020.194)
