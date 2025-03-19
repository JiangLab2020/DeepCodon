from typing import Dict, List, Tuple


def get_min_max_percentage(
    dna: str,
    codon_frequencies: Dict[str, Tuple[List[str], List[float]]],
    window_size: int = 18,
) -> List[float]:
    """
    Calculate the %MinMax metric for a DNA sequence.

    Args:
        dna (str): The DNA sequence.
        codon_frequencies (Dict[str, Tuple[List[str], List[float]]]): Codon frequency distribution per amino acid.
        window_size (int): Size of the window to calculate %MinMax.

    Returns:
        List[float]: List of %MinMax values for the sequence.

    Credit: https://github.com/chowington/minmax
    """
    # Get a dictionary mapping each codon to its respective amino acid
    codon2amino = {
        codon: amino
        for amino, (codons, frequencies) in codon_frequencies.items()
        for codon in codons
    }

    min_max_values = []
    codons = [dna[i : i + 3] for i in range(0, len(dna), 3)]  # Split DNA into codons

    # Iterate through the DNA sequence using the specified window size
    for i in range(len(codons) - window_size + 1):
        codon_window = codons[
            i : i + window_size
        ]  # List of the codons in the current window

        Actual = 0.0  # Average of the actual codon frequencies
        Max = 0.0  # Average of the min codon frequencies
        Min = 0.0  # Average of the max codon frequencies
        Avg = 0.0  # Average of the averages of all the frequencies associated with each amino acid

        # Sum the frequencies for codons in the current window
        for codon in codon_window:
            aminoacid = codon2amino[codon]
            frequencies = codon_frequencies[aminoacid][1]
            codon_index = codon_frequencies[aminoacid][0].index(codon)
            codon_frequency = codon_frequencies[aminoacid][1][codon_index]

            Actual += codon_frequency
            Max += max(frequencies)
            Min += min(frequencies)
            Avg += sum(frequencies) / len(frequencies)

        # Divide by the window size to get the averages
        Actual = Actual / window_size
        Max = Max / window_size
        Min = Min / window_size
        Avg = Avg / window_size

        # Calculate %MinMax
        percentMax = ((Actual - Avg) / (Max - Avg)) * 100
        percentMin = ((Avg - Actual) / (Avg - Min)) * 100

        # Append the appropriate %MinMax value
        if percentMax >= 0:
            min_max_values.append(percentMax)
        else:
            min_max_values.append(-percentMin)

    # Populate the last floor(window_size / 2) entries of min_max_values with None
    for i in range(int(window_size / 2)):
        min_max_values.append(None)

    return min_max_values
