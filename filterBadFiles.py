# Purpose: filter out the bad files found from log files, put them in a
#          badfiles_tmp.txt (so can remove from original file lists) and
#          move them to badfiles_tmp subdirectory (for later removal)

import os
import re

# Define the directory where your log files are located
log_directory = 'logs'

# Create a list to store the bad file names
bad_files = []

# Iterate through each log file
for filename in os.listdir(log_directory):
    if filename.startswith('log_') and filename.endswith('v34.txt'):
        with open(os.path.join(log_directory, filename), 'r') as log_file:
            lines = log_file.readlines()
            for line in lines:
                # Check if the line contains "/media/DATA" and "Error"
                if '/media/DATA' in line and (('Error' in line) or ('File: /media/DATA2/NANO_MC' in line)):
                    # Extract the file name from the line
                    file_match = re.search(r'(/media/DATA\S+\.(?:root|txt))', line)
                    if file_match:
                        bad_file = file_match.group(1)
                        bad_files.append(bad_file)

# Create the 'badfiles_tmp' directory if it doesn't exist
badfiles_tmp_directory = 'badfiles_tmp'
os.makedirs(badfiles_tmp_directory, exist_ok=True)

# Write the bad file names to a separate file
with open(os.path.join(badfiles_tmp_directory, 'badFiles_tmp.txt'), 'w') as bad_files_file:
    bad_files_file.write('\n'.join(bad_files))

print(f"{len(bad_files)} bad files found and listed in 'badFiles_tmp.txt'. They have been moved to the 'badfiles_tmp' directory.")

