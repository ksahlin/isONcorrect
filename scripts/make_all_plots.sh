

EVAL_BASE="/Users/kxs624/Documents/workspace/isONcorrect/"
FILE_BASE="/Users/kxs624/tmp/ISONCORRECT/RESULTS_2020_05_17/"


# Fig 1a, 3a,3b
#python $EVAL_BASE/evaluation/plots.py  $FILE_BASE/isONcorrect_dros_full.csv  $FILE_BASE/plots/ drosophila

# Fig 1b
#python $EVAL_BASE/evaluation/plots.py  $FILE_BASE/isONcorrect_sirv_full.csv  $FILE_BASE/plots/ sirv


# Fig 2a
# python $EVAL_BASE/evaluation_sim/plot_error_rates.py  $FILE_BASE/isONcorrect_sim_error_rate.csv  $FILE_BASE/plots/figure_2a.pdf 

# Fig 2b
#python $EVAL_BASE/evaluation_sirv/plots.py  $FILE_BASE/isONcorrect_sirv_subsampling.csv  $FILE_BASE/plots/ 

# Fig 4
#python $EVAL_BASE/evaluation_sim/plot_abundance_diff.py  $FILE_BASE/isONcorrect_sim_overcorrection.csv  $FILE_BASE/plots/figure_4.pdf 

# Fig Sx
# python $EVAL_BASE/evaluation_sim/plots.py  $FILE_BASE/isONcorrect_sim_full.csv  $FILE_BASE/plots/ sim


# Fig Sy
#python $EVAL_BASE/evaluation_sim/plot_abundance_diff.py  $FILE_BASE/comp_sim_overcorrection.csv  $FILE_BASE/plots/figure_Sy.pdf 
# Fig Sz
#python $EVAL_BASE/evaluation_sim/plot_error_rates.py  $FILE_BASE/comp_sim_error_rate.csv  $FILE_BASE/plots/figure_Sz.pdf 


# Fig S_dros
# python $EVAL_BASE/evaluation/plots.py  $FILE_BASE/comp_dros_full.csv  $FILE_BASE/plots_comp/ drosophila


# Fig sim 4%
# python $EVAL_BASE/evaluation_sim/plot_error_rates.py  $FILE_BASE/4_results.csv  $FILE_BASE/plots/sim/figure_4_percent.pdf 
# python $EVAL_BASE/evaluation_sim/plot_abundance_diff.py  $FILE_BASE/4_abundance.csv  $FILE_BASE/plots/sim/figure_4_percent_overcorr.pdf 

# Fig sim 12%
# python $EVAL_BASE/evaluation_sim/plot_error_rates.py  $FILE_BASE/12_results.csv  $FILE_BASE/plots/sim/figure_12_percent.pdf 
# python $EVAL_BASE/evaluation_sim/plot_abundance_diff.py  $FILE_BASE/12_abundance.csv  $FILE_BASE/plots/sim/figure_12_percent_overcorr.pdf 

# Fig SIRV iso_cov plots
python $EVAL_BASE/evaluation_sirv_iso_cov/plot.py $FILE_BASE/isONcorrect_sirv_iso_cov.csv $FILE_BASE/plots/sim/
python $EVAL_BASE/evaluation_sirv_iso_cov/plot_overcorrected_isoforms.py $FILE_BASE/overcorrected_isoforms.csv $FILE_BASE/plots/sim/
