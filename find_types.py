def extract_columns(line):
    # Split the line by whitespace and return the first and second elements
    columns = line.strip().split()
    if len(columns) >= 2:
        return columns[0], columns[1]
    return None, None

def compare_files(file1_path, file2_path):
    differences = []  # List to store the differences found

    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        # Read lines from both files
        lines_file1 = file1.readlines()
        lines_file2 = file2.readlines()

        pairs_file1 = []
        pairs_file2 = []
        
        # Fill the pairs
        for line in lines_file1:
            column1, column2 = extract_columns(line)
            pairs_file1.append((column1, column2))
            
        for line in lines_file2:
            column1, column2 = extract_columns(line)
            pairs_file2.append((column1, column2))
        
        # Compare the pairs
        for i in range(len(pairs_file1)):
            for j in range(len(pairs_file2)):
                if pairs_file1[i][1] == pairs_file2[j][1]:
                    if pairs_file1[i][0] != pairs_file2[j][0]:
                        differences.append((i, pairs_file1[i][0], pairs_file1[i][1], pairs_file2[j][0], pairs_file2[j][1]))
                        break

    return differences

# Example usage
file1_path = "interface/DijetHistosFill.h"
file2_path = "foo.h"
differences = compare_files(file1_path, file2_path)
if differences:
    print("Differences found:")
    print(f"File 1: {file1_path} File 2: {file2_path}")
    for line_num, col1_file1, col2_file1, col1_file2, col2_file2 in differences:
        print(f"Line {line_num}: File 1: {col1_file1} {col2_file1}, File 2: {col1_file2} {col2_file2}")
else:
    print("No differences found.")
