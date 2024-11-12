# ___________________________________________________________________
#
# Script to change include file in main.cpp
#
# ___________________________________________________________________

import argparse

main_path = 'main.cpp'
include_name = 'include.hpp'

# ___________________________________________________________________
# Parse arguments

parser = argparse.ArgumentParser(description='Change include file in main.cpp')
parser.add_argument('-i', '--include', type=str, default=None, help='Name of the include file')
parser.add_argument('-m', '--main', type=str, default=None, help='Path of the main file')
args = parser.parse_args()

# get arguments
if args.include is not None:
    include_name = args.include

if args.main is not None:
    main_path = args.main

# ___________________________________________________________________
# Summary

print('include: {}'.format(include_name))
print('main: {}'.format(main_path))

# ___________________________________________________________________

# read the main file
with open(main_path, 'r') as file:
    main_file = file.readlines()

# find where the setup include is located
start_index = 0
while("Load a setup" not in main_file[start_index]):
    start_index += 1

print(' > start index: {}'.format(start_index))

# delete all include
for iline in range(start_index+1, len(main_file), 1):
    if '#include' in main_file[iline]:
        main_file[iline] = '\n'

main_file[start_index+1] = '#include "{}" \n'.format(include_name)

# ___________________________________________________________________
# Save the new main file

with open(main_path, 'w') as file:
    file.writelines(main_file)
