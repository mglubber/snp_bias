awk '{print $11"_"$12}' allsnps | sort | uniq -c | sort -nr > snps_type_count 
