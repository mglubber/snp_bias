awk '{print $10" "$11"_"$12}' allsnps | sort | uniq -c | sort -nr > snps_mods_count 
