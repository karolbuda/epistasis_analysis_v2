#!/bin/bash

# Specify the directory where the "first_order_box_plot.svg" and "higher_order_box_plot.svg" files will be copied.
box_dest="/Users/karolbuda/Library/CloudStorage/OneDrive-UBC/PhD/Epistasis_Lit/figures/comms-rev/Box plot files"

# Specify the directory where the "traj_observed.svg" files will be copied.
traj_dest="/Users/karolbuda/Library/CloudStorage/OneDrive-UBC/PhD/Epistasis_Lit/figures/comms-rev/Landscape files"

# Search for "first_order_box_plot.svg" and "higher_order_box_plot.svg" files in subdirectories.
find . -type f -name "first_order_box_plot.svg" -o -name "higher_order_box_plot.svg" |
while read -r file; do
    # Get the subdirectory name and concatenate it with the filename.
    subdirectory="$(dirname "$file")"
    filename="$(basename "$file")"
    new_filename="${subdirectory}_$filename"

    # Copy the file to the specified destination directory.
    cp "$file" "$box_dest/$new_filename"
done

# Search for "traj_observed.svg" files in subdirectories.
find . -type f -name "traj_observed.svg" |
while read -r file; do
    # Get the subdirectory name and concatenate it with the filename.
    subdirectory="$(dirname "$file")"
    filename="$(basename "$file")"
    new_filename="${subdirectory}_$filename"

    # Copy the file to the specified destination directory.
    cp "$file" "$traj_dest/$new_filename"
done

echo "Files copied successfully."