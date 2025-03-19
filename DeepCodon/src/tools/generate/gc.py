def calculate_gc(sequence):
    count_g = 0
    count_c = 0
    for char in sequence:
        if char == "G":
            count_g += 1
        elif char == "C":
            count_c += 1
    return (count_g + count_c) / len(sequence)
