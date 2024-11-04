#!/bin/bash
# List of folders to iterate through
folders=("expAE02" "expAE03" "expAE04" "expAE05" "ctrlAE")

# Loop through each folder
for folder in "${folders[@]}"; do
    # Input NetCDF files
    orog_file="${folder}_16/orog_AIS_UCM_Yelmo_${folder}.nc"
    lithk_file="${folder}_16/lithk_AIS_UCM_Yelmo_${folder}.nc"
    output_file="${folder}_16/base_AIS_UCM_Yelmo_${folder}.nc"

    # Step 1: Merge the orog and lithk files into one file
    # ncks -A $lithk_file $orog_file -o $merged_file
    cdo sub $orog_file $lithk_file $output_file
    ncrename -v orog,base $output_file


    # Step 3: Modify the standard_name attribute to "ice_base"
    ncatted -a standard_name,base,o,c,"ice_base" $output_file

    # Step 4: Modify the long_name attribute to "ice base by subtracting thickness from surface"
    ncatted -a long_name,base,o,c,"ice base by subtracting thickness from surface" $output_file

    # Step 5: Ensure the variable "base" is of type float (with 5 significant digits precision)
    ncks -O --ppc default=5 $output_file $output_file

      # Print a completion message
    echo "Processed folder '$folder' and created output file '$output_file'"
done
# Print final completion message
echo "All folders processed." 

