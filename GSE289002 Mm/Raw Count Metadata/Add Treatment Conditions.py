import pandas as pd
import re


def parse_soft_file(soft_file_path):
    """Parse the SOFT file to extract mouse sample information"""
    samples = []
    current_sample = None

    with open(soft_file_path, 'r') as file:
        for line in file:
            line = line.strip()

            # Start of a new sample
            if line.startswith('^SAMPLE ='):
                if current_sample:
                    samples.append(current_sample)
                sample_id = line.split('=')[1].strip()
                current_sample = {'geo_accession': sample_id}

            # Sample attributes
            elif line.startswith('!Sample_') and current_sample:
                parts = line.split('=', 1)
                if len(parts) == 2:
                    key = parts[0].strip()[8:]  # Remove '!Sample_'
                    value = parts[1].strip()

                    # For characteristics, parse the field: value format
                    if key.startswith('characteristics_ch1'):
                        match = re.match(r'(.+?):\s*(.+)', value)
                        if match:
                            char_name, char_value = match.groups()
                            current_sample[char_name.strip()] = char_value.strip()
                    else:
                        current_sample[key] = value

    # Add the last sample
    if current_sample:
        samples.append(current_sample)

    # Filter for mouse samples only
    mouse_samples = [s for s in samples if s.get('organism_ch1') == 'Mus musculus']
    return mouse_samples


def create_metadata_df(mouse_samples):
    """Create a metadata dataframe from mouse samples"""
    metadata = []

    for sample in mouse_samples:
        record = {
            'geo_accession': sample.get('geo_accession'),
            'title': sample.get('title'),
            'treatment': sample.get('treatment'),
            'sex': sample.get('Sex'),
            'region': sample.get('region'),
            'strain': sample.get('strain'),
            'batch': sample.get('batch')
        }
        metadata.append(record)

    return pd.DataFrame(metadata)


def main():
    # File paths
    soft_file_path = '/GSE289002 Mm/GSE289002_family.soft'
    counts_file_path = '/GSE289002 Mm/GSE289002_mouse_raw_counts.csv'

    # Parse SOFT file and extract mouse samples
    mouse_samples = parse_soft_file(soft_file_path)
    print(f"Mouse samples found: {len(mouse_samples)}")

    # Create metadata dataframe
    metadata_df = create_metadata_df(mouse_samples)

    # Load counts file to examine structure
    counts_df = pd.read_csv(counts_file_path, nrows=5)
    print(f"Counts file columns: {len(counts_df.columns)}")
    print(f"First few columns: {counts_df.columns[:5].tolist()}")

    # Save metadata to file
    output_path = '/GSE289002 Mm/mouse_metadata.csv'
    metadata_df.to_csv(output_path, index=False)
    print(f"Metadata saved to: {output_path}")

    # Print summary
    print("\nSample metadata preview:")
    print(metadata_df.head())


if __name__ == "__main__":
    main()