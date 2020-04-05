

EVAL_BASE="/Users/kxs624/Documents/workspace/isONcorrect/"
FILE_BASE="/Users/kxs624/tmp/ISONCORRECT/RESULTS_2020_04_03/"


# Fig 1a, 3a,3b
python $EVAL_BASE/evaluation/plots.py  $FILE_BASE/isONcorrect_dros_full.csv  $FILE_BASE/plots/ drosophila

# Fig 1b
python $EVAL_BASE/evaluation/plots.py  $FILE_BASE/isONcorrect_sirv_full.csv  $FILE_BASE/plots/ sirv


# Fig 2a
#python $EVAL_BASE/evaluation_sim/plot_error_rates.py  $FILE_BASE/isONcorrect_sim_error_rate.csv  $FILE_BASE/plots/figure_2a.pdf 

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