while IFS=, read -r eps minsamp alpha; do 
   python ./clim_dist/clim_dist_analysis.py $eps $minsamp $alpha
done < clust_param_vals.csv
