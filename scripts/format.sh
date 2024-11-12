# Bash script for clang-format
# Look for all .h and .cpp files in the src directory and format them
# with clang-format

# Get all .cpp and .hpp in src directory
FILES=$(find ./src -type f \( -name "*.cpp" -o -name "*.hpp" \))

# Format all files
for FILE in $FILES
do
  echo "Formatting $FILE"
  clang-format -i $FILE
done
