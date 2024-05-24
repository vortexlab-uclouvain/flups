import os
import re

data_ref_folder = 'data_ref'
data_folder = 'data'
line_pattern = '9 '  # Pattern to search for in the line

# Get the list of files in the data folder
data_files = os.listdir(data_folder)

# Iterate over the data files
for file_name in data_files:
        # Skip files that do not have 'typeGreen=0' in their name
    if '44' not in file_name:
        continue
    
    print(file_name)
    new_file_path = os.path.join(data_ref_folder, file_name)
    ref_file_path = os.path.join(data_folder, file_name)

    # Read the reference file
    with open(ref_file_path, 'r') as ref_file:
        reference_lines = ref_file.readlines()
        print(reference_lines)

    # Read the new results file
    with open(new_file_path, 'r') as new_file:
        new_results_lines = new_file.readlines()
        print(new_results_lines)

    # Find the line to update in the new results
    new_number = None
    for i, line in enumerate(new_results_lines):
        if re.search(line_pattern, line):
            new_number = i
            print(new_number)
            break
        
    ref_number = None
    for i, line in enumerate(reference_lines):
        if re.search(line_pattern, line):
            ref_number = i
            print(ref_number)
            break

    # Update the specific line in the new results
    if new_number is not None:
        new_results_lines[new_number] = reference_lines[ref_number]
        print(f"New results: {new_results_lines}")

        # Write the updated new results back to the file
        with open(new_file_path, 'w') as new_file:
            new_file.writelines(new_results_lines)
    else:
        print(f"Line not found in {file_name}.")